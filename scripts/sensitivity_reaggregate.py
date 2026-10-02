#!/usr/bin/env python3
"""
sensitivity_reaggregate.py — scDRS 敏感性层（细胞级重聚合，零重打分）

为什么可以这样做
----------------
2026-10-01 实测（sgn_atlas，400 细胞子集 vs 全图谱，同 n_ctrl=1000 / 同 seed）：
    raw_score   max|diff| = 7.2e-09
    norm_score  max|diff| = 4.8e-02
即 **逐细胞打分与"本次纳入哪些细胞"无关**（raw_score 逐位一致；norm_score 的
残差来自 scDRS 内部的分布匹配步骤以当次纳入细胞的 score 分布为参照）。

因此所有**只改变细胞集合**的敏感性分析，都不必重新打分，直接对已落盘的
score_<trait>.csv.gz 重新聚合即可：
  - 留一数据集（leave-one-dataset-out, LODO）
  - 仅对照细胞（剔除 Cisplatin / Noise 等处理条件，排除"处理效应冒充遗传富集"）
  - 留一条件（处理组内部一致性）
  - 对照基因集数降采样（n_ctrl 1000 → 200，检验 MC p 的分辨率依赖）

⚠️ 边界（必须写入 Methods / Supplementary）
------------------------------------------
1. 本脚本**不做**任何改变基因集成员的操作。基因集规模（500/2000）与同源映射
   替换会改变打分输入，必须走 run_scdrs_pipeline.py 重跑。
2. 预处理（normalize / HVG / PCA / neighbors / scdrs.pp.preprocess）仍是在**全图谱**
   上完成的，本脚本不重做。故 LODO 检验的是"结论是否由某一数据集的细胞驱动"，
   而非"完全独立地在该数据集上重跑一遍流水线"。这一区别在论文中须写明。
3. n_ctrl 降采样取**前 N 个对照列**。对照基因集是按 random_seed=0 依次抽样的
   i.i.d. 随机集，取前 N 个与"以 n_ctrl=N 重跑"在分布上等价（顺序固定），
   但并非逐位相同；仅作分辨率敏感性检查，不替换主检验的 MC p。

统计量（与 compare_endpoints.py 严格一致，禁止改动）
----------------------------------------------------
    组内 norm_score 的 95 分位数 q95
    端点差 D = q95_A − q95_B
    零分布   = ctrl_q95_A[k] − ctrl_q95_B[l]（全枚举，上限 --n-null）
    mc_p    = (1 + #{|null| ≥ |D|}) / (1 + |null|)，双侧
    BH-FDR 跨（细胞类型 × 端点对）

用法
----
python scripts/sensitivity_reaggregate.py \
  --dir   analysis/scdrs_main/whole_cochlea_atlas \
  --atlas analysis/scdrs_main/whole_cochlea_atlas/preprocessed.h5ad \
  --out   analysis/sensitivity/reaggregate_whole \
  --group-col cell_type_analysis \
  --filter "all" --filter "dataset!=Miao2024" --filter "dataset!=Sun2023_full" \
  --filter "condition==Control" --filter "ctrl_sub=200"
"""
import argparse
import gzip
import os
import numpy as np
import pandas as pd

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import anndata as ad

Q = 0.95
MAX_NULL = 1_000_000


def bh_fdr(p):
    p = np.asarray(p, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(q, 0, 1)
    return out


def read_q95(path, labels, mask_cells, min_cells, ctrl_cols_keep=None):
    """读一个端点的 score csv → (obs q95, ctrl q95 DataFrame, 各类型细胞数)"""
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split(",")
    dtype = {c: np.float32 for c in header[1:]}
    d = pd.read_csv(path, index_col=0, dtype=dtype)
    ctrl_cols = [c for c in d.columns if c.startswith("ctrl_norm_score_")]
    if ctrl_cols_keep is not None and ctrl_cols_keep < len(ctrl_cols):
        # 对照集是 i.i.d. 抽样且顺序固定，取前 N 个 = n_ctrl 降采样
        ctrl_cols = sorted(ctrl_cols, key=lambda c: int(c.rsplit("_", 1)[1]))[:ctrl_cols_keep]
    use = ["norm_score"] + ctrl_cols
    d = d[use]
    ct = labels.reindex(d.index)
    # mask_cells 以图谱 obs_names 为索引，必须按 d.index 对齐后再取值；
    # 直接用 pd.Series(mask_cells, index=d.index) 会按位置套用，顺序不一致时静默错位。
    m = mask_cells.reindex(d.index).fillna(False).astype(bool)
    d = d[m.values & ct.notna().values]
    d["_ct"] = ct.reindex(d.index).values
    n_by_ct = d["_ct"].value_counts()
    keep = n_by_ct[n_by_ct >= min_cells].index
    d = d[d["_ct"].isin(keep)]
    g = d.groupby("_ct")[use].quantile(Q)
    if isinstance(g.index, pd.MultiIndex):
        g.index = g.index.droplevel(1)
    g.index.name = "cell_type"
    del d
    return g["norm_score"], g[ctrl_cols], n_by_ct


def parse_filter(spec, obs):
    """
    返回 (标签名, 细胞布尔掩码 Series, 对照列数上限 or None)
    支持语法：
      all                → 全细胞
      <col>==<val>       → 保留该值
      <col>!=<val>       → 剔除该值
      <col> in a|b|c     → 保留多值（| 分隔）
      ctrl_sub=N         → 对照列降采样到 N
      <expr>+ctrl_sub=N  → 组合（+ 连接）
    """
    parts = [s.strip() for s in spec.split("+") if s.strip()]
    mask = pd.Series(True, index=obs.index)
    n_ctrl = None
    label_parts = []
    for p in parts:
        if p == "all":
            label_parts.append("all")
            continue
        if p.startswith("ctrl_sub="):
            n_ctrl = int(p.split("=", 1)[1])
            label_parts.append(f"ctrl{n_ctrl}")
            continue
        if "!=" in p:
            col, val = [x.strip() for x in p.split("!=", 1)]
            mask &= obs[col].astype(str).ne(val)
            label_parts.append(f"{col}≠{val}")
        elif "==" in p:
            col, val = [x.strip() for x in p.split("==", 1)]
            mask &= obs[col].astype(str).eq(val)
            label_parts.append(f"{col}={val}")
        elif " in " in p:
            col, val = [x.strip() for x in p.split(" in ", 1)]
            vals = [v.strip() for v in val.split("|")]
            mask &= obs[col].astype(str).isin(vals)
            label_parts.append(f"{col}∈{'+'.join(vals)}")
        else:
            raise SystemExit(f"无法解析过滤表达式: {p}")
    name = " & ".join(label_parts) if label_parts else "all"
    return name, mask, n_ctrl


def run_config(scores, labels, mask, traits, min_cells, n_ctrl, n_null, seed):
    obs_q95, ctrl_q95, n_by_ct = {}, {}, None
    for t in traits:
        o, c, n = read_q95(scores[t], labels, mask, min_cells, n_ctrl)
        obs_q95[t], ctrl_q95[t] = o, c
        if n_by_ct is None:
            n_by_ct = n
    cts = sorted(set.intersection(*[set(obs_q95[t].index) for t in traits]))
    recs = []
    rng = np.random.default_rng(seed)
    for ct in cts:
        n = int(n_by_ct.get(ct, 0))
        row = {"cell_type": ct, "n_cell": n}
        for t in traits:
            s = float(obs_q95[t][ct])
            c = ctrl_q95[t].loc[ct].to_numpy(dtype=float)
            row[f"q95_{t}"] = s
            row[f"mcz_{t}"] = (s - c.mean()) / c.std()
        for i in range(len(traits)):
            for j in range(i + 1, len(traits)):
                A, B = traits[i], traits[j]
                D = float(obs_q95[A][ct] - obs_q95[B][ct])
                cA = ctrl_q95[A].loc[ct].to_numpy(dtype=float)
                cB = ctrl_q95[B].loc[ct].to_numpy(dtype=float)
                if len(cA) * len(cB) <= n_null:
                    null = (cA[:, None] - cB[None, :]).ravel()
                else:
                    k = int(np.sqrt(n_null))
                    null = (rng.choice(cA, k)[:, None] - rng.choice(cB, k)[None, :]).ravel()
                p = (1 + np.count_nonzero(np.abs(null) >= abs(D))) / (1 + null.size)
                recs.append(dict(cell_type=ct, n_cell=n, ep_A=A, ep_B=B,
                                 q95_A=float(obs_q95[A][ct]), q95_B=float(obs_q95[B][ct]),
                                 diff=D, direction=("A>B" if D > 0 else "A<B"),
                                 mc_p=float(p), n_null=int(null.size)))
    pw = pd.DataFrame(recs)
    if len(pw):
        pw["qval_BH"] = bh_fdr(pw["mc_p"].values)
        pw = pw.sort_values("mc_p").reset_index(drop=True)
    return pw, cts, n_by_ct


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True, help="含 score_<trait>.csv.gz 的目录")
    ap.add_argument("--atlas", required=True, help="preprocessed.h5ad（取 obs）")
    ap.add_argument("--out", required=True)
    ap.add_argument("--group-col", default="cell_type_analysis")
    ap.add_argument("--filters", nargs="+", required=True)
    ap.add_argument("--traits", default=None)
    ap.add_argument("--min-cells", type=int, default=30)
    ap.add_argument("--n-null", type=int, default=MAX_NULL)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--focus", default=None,
                    help="重点关注条目，如 'IPhC_IBC|SEN|CON'；汇总时置顶")
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)

    scores = {}
    for fn in sorted(os.listdir(args.dir)):
        if fn.startswith("score_") and fn.endswith(".csv.gz"):
            scores[fn[len("score_"):-len(".csv.gz")]] = os.path.join(args.dir, fn)
    traits = args.traits.split(",") if args.traits else sorted(scores)
    print(f"端点: {traits}")

    a = ad.read_h5ad(args.atlas, backed="r")
    obs = a.obs
    col = args.group_col
    if col not in obs.columns:
        for _alt in ("cell_type_transferred", "cell_type"):
            if _alt in obs.columns:
                print(f"使用备选分组列 '{_alt}'")
                col = _alt
                break
        else:
            raise SystemExit(f"图谱无 '{col}' 列；候选: {[c for c in obs.columns if 'cell_type' in c]}")
    labels = obs[col].astype(str)
    print(f"图谱: {a.n_obs:,} 细胞，分组列 '{col}'，{labels.nunique()} 类")

    all_pw, summary_rows = {}, []
    for spec in args.filters:
        name, mask, n_ctrl = parse_filter(spec, obs)
        print(f"\n===== 配置 [{name}] 细胞 {int(mask.sum()):,} / {len(mask):,}；"
              f"对照列 {n_ctrl or '全部'} =====", flush=True)
        pw, cts, n_by_ct = run_config(scores, labels, mask, traits,
                                      args.min_cells, n_ctrl, args.n_null, args.seed)
        tag = spec.replace(" ", "").replace("==", "-").replace("!=", "-not-")
        tag = tag.replace("+", "_").replace("|", "-").replace("/", "-")
        pw.to_csv(os.path.join(args.out, f"pairwise_{tag}.csv"), index=False)
        all_pw[name] = pw
        sig = pw[pw["qval_BH"] < 0.05] if len(pw) else pw
        print(f"  进入检验细胞类型 {len(cts)}；BH-FDR<0.05: {len(sig)} 项", flush=True)
        for _, r in (sig.iterrows() if len(sig) else []):
            summary_rows.append(dict(config=name, cell_type=r["cell_type"],
                                     pair=f"{r['ep_A']}>{r['ep_B']}", diff=r["diff"],
                                     mc_p=r["mc_p"], qval_BH=r["qval_BH"]))
        # 关注条目（无论是否显著都记录）
        if args.focus and len(pw):
            ct_, a_, b_ = args.focus.split("|")
            hit = pw[(pw["cell_type"] == ct_) &
                     (((pw["ep_A"] == a_) & (pw["ep_B"] == b_)) |
                      ((pw["ep_A"] == b_) & (pw["ep_B"] == a_)))]
            if len(hit):
                r = hit.iloc[0]
                summary_rows.append(dict(config=name, cell_type=ct_, pair=f"{a_}vs{b_}[FOCUS]",
                                         diff=r["diff"] if r["ep_A"] == a_ else -r["diff"],
                                         mc_p=r["mc_p"], qval_BH=r["qval_BH"]))

    foc = pd.DataFrame(summary_rows)
    foc.to_csv(os.path.join(args.out, "focus_summary.csv"), index=False)

    # ---------- 汇总 ----------
    L = ["# 敏感性层（细胞级重聚合，零重打分）", "",
         f"- 打分目录: `{args.dir}`",
         f"- 图谱: `{os.path.basename(args.atlas)}`（{a.n_obs:,} 细胞），分组列 `{col}`",
         f"- 端点: {', '.join(traits)}；统计量: 组内 norm_score 的 {Q:.0%} 分位数，"
         "零分布 = 对照基因集同位置分位数之差（双侧 MC p），跨（细胞类型 × 端点对）BH-FDR",
         "", "## 方法学边界（写进 Methods 时不得删）", "",
         "- 逐细胞打分与纳入细胞集合无关（实测 raw_score 逐位一致，Δ=7.2e-09），",
         "  故仅改变细胞集合的敏感性分析可直接重聚合已落盘分数。",
         "- 但预处理（归一化 / HVG / PCA / neighbors / scdrs.pp.preprocess）仍在**全图谱**",
         "  上完成，本层**未**在各子集上独立重跑流水线。LODO 回答的是",
         "  “结论是否由某一数据集的细胞驱动”，不是“独立重跑”。",
         "- 本层不改变基因集成员；基因集规模与同源映射的敏感性须重跑打分。", ""]
    if args.focus and len(foc):
        ct_, a_, b_ = args.focus.split("|")
        f2 = foc[foc["pair"].str.contains("FOCUS")]
        L += ["## 主结论稳健性（关注条目）", "",
              f"- 关注: `{ct_}`，对比 {a_} vs {b_}", "",
              f2.round(6).to_markdown(index=False) if len(f2) else "（无）", "",
              f"- **判据**：所有配置下 `diff` 同号且 `qval_BH` 仍 <0.05 → 主结论稳健；",
              "  任一配置翻转符号或失去显著性 → 须在正文明确标注该稳健性边界。", ""]
    L += ["## 各配置下 BH-FDR<0.05 的全部条目", ""]
    sig_all = (foc[~foc["pair"].str.contains("FOCUS")] if len(foc)
               else pd.DataFrame(columns=["config", "cell_type", "pair", "diff", "mc_p", "qval_BH"]))
    L += [sig_all.round(6).to_markdown(index=False) if len(sig_all) else "（无）", "",
          "---", "", "逐配置明细见 `pairwise_*.csv`；关注条目汇总见 `focus_summary.csv`。"]
    with open(os.path.join(args.out, "sensitivity_summary.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(L) + "\n")
    print(f"\n[写出] {args.out}/sensitivity_summary.md")
    print("DONE")


if __name__ == "__main__":
    main()
