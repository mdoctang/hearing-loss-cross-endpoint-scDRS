#!/usr/bin/env python3
"""
compare_endpoints.py — 跨端点细胞类型富集谱对比（CAND-C v2 核心增量）

设计依据：study_design.md §6 第 4 条
  "各端点在同一整合图谱上的富集谱（细胞类型 × 端点矩阵）→ 配对对比检验"

⚠️ 统计量选择（重要，v1 → v2 的关键修正）
------------------------------------------------
scDRS 的官方分组统计量（scdrs/method.py:791-801）**不是**组内 norm_score 的均值，
而是**组内 norm_score 的 95 分位数**，且其零分布由**同一细胞类型下 1000 个对照
基因集的同位置 95 分位数**构造：

    score_q95      = quantile(group_norm_score, 0.95)
    v_ctrl_q95     = quantile(group_ctrl_norm_score, 0.95, axis=0)     # 长度 n_ctrl
    assoc_mcp      = ((v_ctrl_q95 >= score_q95).sum() + 1) / (n_ctrl + 1)
    assoc_mcz      = (score_q95 - v_ctrl_q95.mean()) / v_ctrl_q95.std()

**若改用"逐细胞均值的配对比 Z 检验"会得到大量假阳性**（实测：两个完全随机的
MOCK 基因集之间即报出 3 项 FDR<0.05）。原因是单个随机基因集在特定细胞类型中
仍有系统性表达偏倚，只有对照基因集校准才能扣除它。因此本脚本 v2 严格沿用
scDRS 的 95 分位数 + 对照校准口径，并把它**扩展到两端点之差**：

    D          = q95_A(group) − q95_B(group)                     ← 观测差异
    null(1000×1000) = ctrl_q95_A[k] − ctrl_q95_B[l]              ← 零分布
    mc_p       = (1 + #{ |null| >= |D| }) / (1 + 1e6)            ← 双侧

统计含义：在扣除"匹配随机基因集"的表达偏倚后，端点 A 与端点 B 在该细胞类型
上的富集强度是否存在差异 —— 即**病因-细胞类型匹配**的直接检验。

输入
----
--dir       一次 scDRS 运行的输出目录（内含 score_<trait>.csv.gz）
--atlas     该次运行所用图谱 h5ad（只读 obs 取细胞类型标签）
--group-col 细胞类型列名（默认 cell_type_analysis）

输出
----
compare_endpoint_matrix.csv  细胞类型 × 端点 富集矩阵（scDRS 官方 assoc_mcz/mcp）
compare_pairwise_tests.csv   每个（细胞类型, 端点对）的校准后差检验 + BH-FDR
compare_summary.md           人类可读汇总
"""
import argparse
import gzip
import os
import numpy as np
import pandas as pd

# BLAS 线程限制（与 run_scdrs_pipeline.py 同因：满线程在 Windows 下概率性段错误）
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import anndata as ad

Q = 0.95          # 与 scDRS 官方一致
MAX_NULL = 1_000_000   # 1000×1000，超出则抽样以控内存


def bh_fdr(p):
    """Benjamini-Hochberg FDR，返回与输入同序的 q 值。"""
    p = np.asarray(p, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(q, 0, 1)
    return out


def group_q95(path, labels, min_cells):
    """
    对单个端点：逐细胞类型计算 norm_score 与全部对照列的 95 分位数。
    返回 (obs: Series[cell_type], ctrl: DataFrame[cell_type × n_ctrl])
    —— 一次只驻留一个端点的矩阵，避免多端点同时入内存。
    """
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split(",")
    # 只对数据列设 float32（不能整体 dtype=float32，否则首列细胞名也会被转数值）
    dtype = {c: np.float32 for c in header[1:]}
    d = pd.read_csv(path, index_col=0, dtype=dtype)
    ctrl_cols = [c for c in d.columns if c.startswith("ctrl_norm_score_")]
    d["_ct"] = labels.reindex(d.index).values
    if d["_ct"].isna().any():
        raise SystemExit("图谱 obs_names 与 score 文件索引无法对齐。")
    d = d[d["_ct"].notna()]
    n_by_ct = d["_ct"].value_counts()
    keep = n_by_ct[n_by_ct >= min_cells].index
    d = d[d["_ct"].isin(keep)]
    q = d.groupby("_ct")[["norm_score"] + ctrl_cols].quantile(Q)
    # pandas 对单键 groupby + 标量 quantile 可能返回单层索引，视版本而定
    if isinstance(q.index, pd.MultiIndex):
        q.index = q.index.droplevel(1)
    q.index.name = "cell_type"
    del d
    return q["norm_score"], q[ctrl_cols], n_by_ct


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True)
    ap.add_argument("--atlas", required=True)
    ap.add_argument("--group-col", default="cell_type_analysis")
    ap.add_argument("--min-cells", type=int, default=30,
                    help="细胞数下限；小于此的细胞类型不检验（scDRS 默认 30）")
    ap.add_argument("--n-null", type=int, default=MAX_NULL,
                    help="零分布抽样上限（默认 1e6 = 1000×1000 全枚举）")
    ap.add_argument("--out-prefix", default=None)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)
    pref = args.out_prefix or os.path.join(args.dir, "compare")

    # ---------- 1. 端点清单 ----------
    scores = {}
    for fn in sorted(os.listdir(args.dir)):
        if fn.startswith("score_") and fn.endswith(".csv.gz"):
            scores[fn[len("score_"):-len(".csv.gz")]] = os.path.join(args.dir, fn)
    if len(scores) < 2:
        raise SystemExit(f"需要 ≥2 个端点，实见 {list(scores)}")
    traits = list(scores)
    print(f"端点: {traits}")

    # ---------- 2. 标签 ----------
    a = ad.read_h5ad(args.atlas, backed="r")
    col = args.group_col
    if col not in a.obs.columns:
        # run_scdrs_pipeline.py 在 preprocess 时构建 cell_type_analysis 并写入
        # <dir>/preprocessed.h5ad；原始 atlas 不一定有该列（如 sgn_atlas 只有
        # cell_type）。优先从预处理缓存读取——它就是打分时实际使用的标签列。
        _pp = os.path.join(args.dir, "preprocessed.h5ad")
        if os.path.exists(_pp):
            a = ad.read_h5ad(_pp, backed="r")
            print(f"原始图谱无 '{col}' 列，改用预处理缓存: {_pp}")
    if col not in a.obs.columns:
        # 仍缺失时按 pipeline 的标签构建逻辑回退（cell_type_transferred 优先）
        for _alt in ("cell_type_transferred", "cell_type"):
            if _alt in a.obs.columns:
                print(f"使用备选分组列 '{_alt}'（与 pipeline 标签构建逻辑一致）")
                col = _alt
                break
        else:
            cands = [c for c in a.obs.columns if "cell_type" in c]
            raise SystemExit(f"图谱无 '{col}' 列；候选: {cands}")
    labels = a.obs[col].astype(str)
    print(f"图谱: {a.n_obs} 细胞，分组列 '{col}'，{labels.nunique()} 类")

    # ---------- 3. 逐端点求 q95 ----------
    obs_q95, ctrl_q95, n_by_ct = {}, {}, None
    for t in traits:
        o, c, n = group_q95(scores[t], labels, args.min_cells)
        obs_q95[t], ctrl_q95[t] = o, c
        if n_by_ct is None:
            n_by_ct = n
        print(f"  {t}: 有效细胞类型 {len(o)}，对照列 {c.shape[1]}")

    # 对齐细胞类型集合
    cts = sorted(set.intersection(*[set(obs_q95[t].index) for t in traits]))
    print(f"进入检验的细胞类型: {len(cts)}")

    # ---------- 4. 富集矩阵（scDRS 官方口径，逐端点复算） ----------
    rows = []
    for ct in cts:
        rec = {"cell_type": ct, "n_cell": int(n_by_ct.get(ct, 0))}
        for t in traits:
            s = float(obs_q95[t][ct])
            c = ctrl_q95[t].loc[ct].to_numpy(dtype=float)
            rec[f"q95_{t}"] = s
            rec[f"mcp_{t}"] = (int((c >= s).sum()) + 1) / (len(c) + 1)
            rec[f"mcz_{t}"] = (s - c.mean()) / c.std()
        rows.append(rec)
    mat = pd.DataFrame(rows).set_index("cell_type")
    mat.to_csv(f"{pref}_endpoint_matrix.csv")
    print(f"[写出] {pref}_endpoint_matrix.csv")

    # ---------- 5. 端点对差异检验（校准后） ----------
    recs = []
    for ct in cts:
        n = int(n_by_ct.get(ct, 0))
        for i in range(len(traits)):
            for j in range(i + 1, len(traits)):
                A, B = traits[i], traits[j]
                D = float(obs_q95[A][ct] - obs_q95[B][ct])
                cA = ctrl_q95[A].loc[ct].to_numpy(dtype=float)
                cB = ctrl_q95[B].loc[ct].to_numpy(dtype=float)
                # 全枚举 nA×nB 差值；过大则抽样
                if len(cA) * len(cB) <= args.n_null:
                    null = (cA[:, None] - cB[None, :]).ravel()
                else:
                    k = int(np.sqrt(args.n_null))
                    null = (rng.choice(cA, k)[:, None] - rng.choice(cB, k)[None, :]).ravel()
                p = (1 + np.count_nonzero(np.abs(null) >= abs(D))) / (1 + null.size)
                recs.append(dict(cell_type=ct, n_cell=n, ep_A=A, ep_B=B,
                                 q95_A=float(obs_q95[A][ct]), q95_B=float(obs_q95[B][ct]),
                                 diff=D, direction=("A>B" if D > 0 else "A<B"),
                                 mc_p=float(p), n_null=int(null.size)))
    pw = pd.DataFrame(recs).sort_values("mc_p").reset_index(drop=True)
    pw["qval_BH"] = bh_fdr(pw["mc_p"].values)
    pw.to_csv(f"{pref}_pairwise_tests.csv", index=False)
    print(f"[写出] {pref}_pairwise_tests.csv")

    # ---------- 6. 汇总 ----------
    sig = pw[pw["qval_BH"] < 0.05]
    L = ["# 跨端点细胞类型富集对比", "",
         f"- 图谱: `{os.path.basename(args.atlas)}`（{a.n_obs:,} 细胞）",
         f"- 分组列: `{col}`；进入检验的细胞类型: {len(cts)}（下限 {args.min_cells} 细胞）",
         f"- 端点: {', '.join(traits)}",
         f"- 统计量: 组内 norm_score 的 {Q:.0%} 分位数，零分布 = 1000 对照基因集同位置分位数之差（双侧 MC p）",
         "", "## 富集矩阵（scDRS 官方口径：mcz = Mc-Z，mcp = MC p）", ""]
    cols = ["n_cell"] + [f"{k}_{t}" for t in traits for k in ("q95", "mcz", "mcp")]
    L += [mat[cols].round(4).to_markdown(), "", "## 配对比检验", ""]
    if len(sig) == 0:
        L += ["**BH-FDR<0.05 的端点间差异：0 项。**", "",
              "> 模拟（随机基因集）运行下这正是**预期结果**，说明检验已正确校准；",
              "> 若为真实端点，则说明未检出病因特异性富集。"]
    else:
        L += [f"**BH-FDR<0.05 的端点间差异：{len(sig)} 项**", "",
              sig.round(6).to_markdown(index=False)]
    L += ["", "## 全部检验（按 mc_p 升序）", "", pw.round(6).to_markdown(index=False), "",
          "---", "",
          "## TODO（正式分析前需补齐）",
          "",
          "- 预登记要求 `block jackknife over genes` 估计富集系数差的标准误；需按基因块",
          "  重算 scDRS 分数（依赖 preprocessed.h5ad + .gs），作为独立脚本补齐。",
          "- 先验三组（血管纹 / 感觉上皮 / SGN）的 Bonferroni 校正尚未施加（见 §6.8）。",
          "- `hetero_mcz`（组内异质性）直接取 scDRS `group_*` 表，本脚本不重算。"]
    with open(f"{pref}_summary.md", "w", encoding="utf-8") as fh:
        fh.write("\n".join(L) + "\n")
    print(f"[写出] {pref}_summary.md")
    print("DONE")


if __name__ == "__main__":
    main()
