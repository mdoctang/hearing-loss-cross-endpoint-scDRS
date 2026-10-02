#!/usr/bin/env python3
"""
jackknife_endpoints.py — 跨端点富集差异的基因块 jackknife 稳健性检查

定位（见 study_design.md §6.4 第二分支，2026-09-30 结果揭晓前范围修订）
----------------------------------------------------------------------
本脚本**不产出 p 值**，只回答一个问题：
    "跨端点富集强度之差，在多大程度上依赖基因集里具体包含哪些基因？"

与主检验的分工：
  - 主检验（compare_endpoints.py）：MC 置换，零分布来自 1000 个匹配随机基因集
    → 回答"该差异相对随机基因集有多意外"（推断性）。
  - 本脚本：基因块 jackknife → 回答"该差异对基因集组成有多敏感"（稳健性）。
  二者互补，**口径不得混用**：本脚本用 n_ctrl=200 降分辨率换取可行性，
  其输出只能作为"低分辨率稳健性检查"呈报。

方法
----
1. 将 .gs 中每个 trait 的 top-N 基因等分为 B 块（默认 20）。
2. 逐块剔除后重算 scDRS 细胞分数（weight_opt="vs"，n_ctrl=200）。
3. 对每个细胞类型、每个端点，取组内 norm_score 的 95 分位数（与 scDRS 官方一致）。
4. 端点对之差 D_b = q95_A,b − q95_B,b。
5. Jackknife：SE = sqrt((B-1)/B · Σ(D_b − D̄)²)，并给出删除-±1 块区间。

成本（实测基准：8,916 细胞 / n_ctrl=1000 / 单 trait ≈ 5 分钟）
----------------------------------------------------------------
跑数 = (B+1) × n_traits。本脚本按预登记收窄为**仅主图谱**（50,757 细胞），
B=20、n_ctrl=200，预计 5–6 小时。**不要**在 whole_cochlea_atlas（101,082
细胞）上直接跑全量，会到数十小时量级。

用法
----
python scripts/jackknife_endpoints.py \
  --atlas analysis/scrna_integrated/sv_atlas_labeled.h5ad \
  --cache analysis/scdrs_main/sv_atlas/preprocessed.h5ad \
  --gs analysis/hl_endpoints.gs \
  --out analysis/jackknife_sv \
  --traits SEN,CON,NOISE --n-blocks 20 --n-ctrl 200
"""
import argparse
import os

# BLAS 线程限制（与 run_scdrs_pipeline.py 同因：满线程在 Windows 下概率性段错误）
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
os.environ.setdefault("TQDM_DISABLE", "1")

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import scdrs

# NumPy 2.x 兼容垫片（与 run_scdrs_pipeline.py 一致；见 analysis/ENVIRONMENT.md）
if not hasattr(np, "float_"):
    np.float_ = np.float64

SCDRS_DATA = os.path.join(os.path.dirname(scdrs.__file__), "data")
HOMOLOG = os.path.join(SCDRS_DATA, "mouse_human_homologs.txt")
Q = 0.95


def log(m):
    print(m, flush=True)


def load_homolog():
    h = pd.read_csv(HOMOLOG, sep="\t")
    h.columns = [c.upper() for c in h.columns]
    h = h.dropna().drop_duplicates(subset=["MOUSE_GENE_SYM"])
    return dict(zip(h["HUMAN_GENE_SYM"], h["MOUSE_GENE_SYM"]))


def parse_gs(gs_path):
    out = {}
    with open(gs_path, encoding="utf-8") as fh:
        fh.readline()
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            items = []
            for tok in parts[1].split(","):
                if tok:
                    g, w = tok.rsplit(":", 1)
                    items.append((g, float(w)))
            out[parts[0]] = items
    return out


def q95_by_group(df_score, labels):
    d = df_score[["norm_score"]].copy()
    d["_ct"] = labels.reindex(d.index).values
    q = d.groupby("_ct")["norm_score"].quantile(Q)
    return q


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--atlas", required=True, help="图谱 h5ad（用于取细胞类型标签）")
    ap.add_argument("--cache", required=True, help="已预处理的 h5ad（含 SCDRS_PARAM）")
    ap.add_argument("--gs", required=True, help="scDRS .gs（人类符号）")
    ap.add_argument("--out", required=True)
    ap.add_argument("--traits", default=None, help="逗号分隔；默认全部")
    ap.add_argument("--n-blocks", type=int, default=20)
    ap.add_argument("--n-ctrl", type=int, default=200,
                    help="jackknife 专用；主检验仍为 1000，口径不得混用")
    ap.add_argument("--group-col", default=None)
    ap.add_argument("--min-cells", type=int, default=30)
    ap.add_argument(
        "--keep-cells", default=None,
        help="逗号分隔的细胞类型；只对这些细胞打分（其余细胞不进入统计）。\n"
             "依据 2026-10-01 实测：scDRS 逐细胞 raw_score 与'本次纳入哪些细胞'无关\n"
             "（400 细胞子集 vs 全图谱，同参数下 raw_score 逐位一致 Δ=7.2e-09），\n"
             "而组内 q95 只用到组内细胞，故可把 101,082 细胞的全图谱收窄到目标细胞类型，\n"
             "使 whole_cochlea 上的基因块 jackknife 从'数十小时'降到可行量级。\n"
             "代价：norm_score 依赖 scDRS 内部的分布匹配（以当次纳入细胞的 score 分布\n"
             "为参照），子集与全图谱会有 ~0.05 的系统性平移；但全量组与删除-1 块组\n"
             "在同一子集内计算，该平移为共模，在差值 D 上基本抵消。本层不产出 p 值，\n"
             "只报告方向一致性与 jackknife SE，故此近似可接受；须在 Methods 写明。",
    )
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    log(f"# 基因块 jackknife（稳健性层，n_ctrl={args.n_ctrl}，B={args.n_blocks}）")

    # ---------- 1. 预处理对象 ----------
    if not os.path.exists(args.cache):
        raise SystemExit(f"找不到预处理缓存 {args.cache}；先跑 run_scdrs_pipeline.py（会写出该文件）")
    adata = sc.read_h5ad(args.cache)
    log(f"预处理对象: {adata.n_obs} 细胞 × {adata.n_vars} 基因")
    # 还原被打包的 pandas.Series：run_scdrs_pipeline.py 写缓存时为绕开
    # anndata>=0.11 的 IORegistryError，把 COV_GENE_MEAN 存成 dict 落盘；
    # 而 _compute_raw_score 需要 Series 的 .loc 索引，不还原会抛
    # AttributeError: 'dict' object has no attribute 'loc'。
    _sp = adata.uns.get("SCDRS_PARAM", {})
    for _k, _v in list(_sp.items()):
        if isinstance(_v, dict) and _v.get("__pd_series__"):
            _sp[_k] = pd.Series(_v["values"], index=_v["index"])
            log(f"uns/SCDRS_PARAM[{_k}] 已从打包结构还原为 Series")
    if "SCDRS_PARAM" not in adata.uns or "GENE_STATS" not in adata.uns["SCDRS_PARAM"]:
        raise SystemExit("缓存内缺 SCDRS_PARAM.GENE_STATS，请用 run_scdrs_pipeline.py 生成缓存。")

    # ---------- 2. 标签 ----------
    a0 = ad.read_h5ad(args.atlas, backed="r")
    col = args.group_col
    if col is None:
        col = "cell_type_transferred" if "cell_type_transferred" in a0.obs.columns else "cell_type"
    labels = a0.obs[col].astype(str).reindex(adata.obs_names)
    if labels.isna().any():
        raise SystemExit(f"标签无法对齐（列 '{col}'）；检查 cache 与 atlas 是否同一图谱。")
    # 只保留目标细胞类型（见 --keep-cells 说明；使大图谱上的 jackknife 可行）
    if args.keep_cells:
        want = set(s.strip() for s in args.keep_cells.split(","))
        sel = labels[labels.isin(want)].index
        miss = want - set(labels.unique())
        if miss:
            log(f"⚠️ 以下细胞类型在图谱中不存在，已忽略: {sorted(miss)}")
        if len(sel) == 0:
            raise SystemExit("--keep-cells 未匹配到任何细胞")
        adata = adata[sel].copy()
        labels = labels.reindex(adata.obs_names)
        log(f"--keep-cells: {len(want)} 类 → 保留 {adata.n_obs:,} 细胞（原全图谱子集化）")

    n_ct = labels.value_counts()
    keep_ct = n_ct[n_ct >= args.min_cells].index
    log(f"分组列 '{col}'；{len(n_ct)} 类，其中 {len(keep_ct)} 类 ≥{args.min_cells} 细胞进入统计")

    # ---------- 3. 基因集 ----------
    h2m = load_homolog()
    gs_all = parse_gs(args.gs)
    traits = args.traits.split(",") if args.traits else list(gs_all.keys())
    log(f"gs traits: {list(gs_all.keys())}；运行: {traits}")

    gene_map = {}
    for t in traits:
        gm, wm = [], []
        for g, w in gs_all[t]:
            mg = h2m.get(g.upper())
            if mg is not None and mg in adata.var_names:
                gm.append(mg)
                wm.append(w)
        gene_map[t] = (gm, wm)
        log(f"  [{t}] 图谱内可用基因 {len(gm)}")
        if len(gm) < 2 * args.n_blocks:
            raise SystemExit(f"[{t}] 可用基因 {len(gm)} 太少，无法切 {args.n_blocks} 块")

    # ---------- 4. 打分：全量 + 逐块剔除 ----------
    def score(gm, wm):
        ds = scdrs.method.score_cell(adata, gm, gene_weight=wm, weight_opt="vs",
                                     n_ctrl=args.n_ctrl, random_seed=0)
        return ds

    qobs, qblocks = {}, {t: [] for t in traits}
    for t in traits:
        gm, wm = gene_map[t]
        log(f"[{t}] 全量打分（{len(gm)} 基因, n_ctrl={args.n_ctrl}）…")
        qobs[t] = q95_by_group(score(gm, wm), labels)

        idx = np.arange(len(gm))
        blocks = np.array_split(idx, args.n_blocks)
        for b, blk in enumerate(blocks):
            dep = np.setdiff1d(idx, blk, assume_unique=False)
            gmb = [gm[i] for i in dep]
            wmb = [wm[i] for i in dep]
            qblocks[t].append(q95_by_group(score(gmb, wmb), labels))
            log(f"  块 {b + 1}/{args.n_blocks}（剔除 {len(blk)} 基因，余 {len(gmb)}）完成")

    # ---------- 5. Jackknife ----------
    # 注意：A、B 两个 trait 的基因集不同，其"第 b 块"并非同一批基因，
    # 强行按块索引配对会引入无关噪声。故**分别扰动两侧基因集**：
    #   SE_A：只剔除 A 的基因块（B 保持全量）→ A 侧基因组成带来的不确定性
    #   SE_B：只剔除 B 的基因块（A 保持全量）→ B 侧基因组成带来的不确定性
    #   SE   = sqrt(SE_A² + SE_B²)
    def jack_se(v):
        v = np.asarray(v, dtype=float)
        nb = len(v)
        if nb < 3:
            return np.nan
        return float(np.sqrt((nb - 1) / nb * np.sum((v - v.mean()) ** 2)))

    recs = []
    for i in range(len(traits)):
        for j in range(i + 1, len(traits)):
            A, B = traits[i], traits[j]
            for ct in keep_ct:
                if ct not in qobs[A].index or ct not in qobs[B].index:
                    continue
                qA_full, qB_full = float(qobs[A][ct]), float(qobs[B][ct])
                D_full = qA_full - qB_full
                DA = [float(qb.get(ct, np.nan)) - qB_full for qb in qblocks[A]]
                DB = [qA_full - float(qb.get(ct, np.nan)) for qb in qblocks[B]]
                DA = [v for v in DA if np.isfinite(v)]
                DB = [v for v in DB if np.isfinite(v)]
                se_a, se_b = jack_se(DA), jack_se(DB)
                if not np.isfinite(se_a) or not np.isfinite(se_b):
                    continue
                se = float(np.sqrt(se_a ** 2 + se_b ** 2))
                allv = DA + DB
                recs.append(dict(cell_type=ct, n_cell=int(n_ct.get(ct, 0)),
                                 ep_A=A, ep_B=B, diff_full=D_full,
                                 jack_se_A=se_a, jack_se_B=se_b, jack_se=se,
                                 ci_lo=D_full - se, ci_hi=D_full + se,
                                 range_lo=float(np.min(allv)), range_hi=float(np.max(allv)),
                                 n_blocks=len(DA),
                                 sign_consistent=bool(np.all(np.sign(allv) == np.sign(D_full))),
                                 n_ctrl_jackknife=args.n_ctrl))
    out = pd.DataFrame(recs).sort_values(["ep_A", "ep_B", "cell_type"])
    out.to_csv(f"{args.out}/jackknife_results.csv", index=False)
    log(f"[写出] {args.out}/jackknife_results.csv")

    # ---------- 6. 汇总 ----------
    L = ["# 基因块 jackknife（稳健性层）", "",
         f"- 图谱: `{os.path.basename(args.atlas)}`",
         f"- 基因集: `{os.path.basename(args.gs)}`；端点: {', '.join(traits)}",
         f"- 参数: B={args.n_blocks}，n_ctrl={args.n_ctrl}（**低分辨率**，勿与主检验的 MC p 混用）",
         "- 本表**不产出 p 值**，只报告差异对基因集组成的敏感性。", "",
         out.round(4).to_markdown(index=False) if len(out) else "（无结果）", "", "---", "",
         "## 读法",
         "",
         "- `diff_full`：全基因集下的端点间 q95 之差（观测值）。",
         "- `jack_se_A` / `jack_se_B`：分别只扰动 A 侧 / B 侧基因块所得的标准误",
         "  （A、B 基因集不同，其\"第 b 块\"不是同一批基因，故不按块索引配对）；",
         "  `jack_se = sqrt(se_A² + se_B²)`。",
         "- `sign_consistent`：**关键判据**。若为 True，说明剔除任意单个基因块都不会翻转",
         "  差异方向——即主结论不依赖某几个基因。若为 False，需在正文中明确该差异不稳健。",
         "- `ci_lo`/`ci_hi` = diff_full ± 1×jack_se，仅作区间展示，非正式置信区间；",
         "  `range_lo`/`range_hi` 为全部删除-1 块估计的极值范围。"]
    with open(f"{args.out}/jackknife_summary.md", "w", encoding="utf-8") as fh:
        fh.write("\n".join(L) + "\n")
    log(f"[写出] {args.out}/jackknife_summary.md")
    if len(out):
        log(f"符号一致率: {out['sign_consistent'].mean():.1%}（{int(out['sign_consistent'].sum())}/{len(out)}）")
    log("DONE")


if __name__ == "__main__":
    main()
