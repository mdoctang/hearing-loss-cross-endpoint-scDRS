#!/usr/bin/env python3
"""
build_scrna_atlas.py — 耳蜗跨病因轴 scRNA 整合图谱构建（CAND-C v2 scRNA 层）

输入（project/datasets/scrna/）:
  sun2023_cra004814_full_atlas.h5ad   衰老全耳蜗（27 类注释，Ensembl+symbol）
  sun2023_cra004814_stria.h5ad        衰老血管纹子集（5 类）
  milon2021_gse168041_lateralwall.h5ad 噪声侧壁（5 类）
  milon2021_gse168041_sgn.h5ad        噪声 SGN（5 类）
  korrapati2019_naive_sv.h5ad         naive SV 对照（部分注释）
  GSE165662/GSM5047933_filtered_gene_bc_matrices_h5.h5  顺铂 SV（无注释）
  GSM8618786/7/8/9_..._filtered_feature_bc_matrix.h5    顺铂/对照全耳蜗管 snRNA（无注释）

输出（project/analysis/scrna_integrated/）:
  sv_atlas.h5ad            血管纹/侧壁图谱（Sun stria + Milon LW + Korrapati + Taukulis）
  whole_cochlea_atlas.h5ad 全耳蜗（Sun 全图谱 + Miao GSE281324）
  sgn_atlas.h5ad           SGN（Milon）
  integration_log.md       逐步计数与细胞类型分布

统一约定:
  - var_names 使用小鼠基因符号（去重）
  - obs 统一字段: dataset, axis, condition, cell_type, cell_type_original, sample, n_genes, n_counts, pct_mt
  - QC: 细胞 ≥200 基因；基因在 ≥3 细胞中检出；pct_mt<20%
  - .X 保留原始计数（稀疏 int32→float32）
"""
import os
import sys
import glob
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy import sparse

BASE = "C:/Users/docta/WorkBuddy/纯生信/project"
SRC = f"{BASE}/datasets/scrna"
OUT = f"{BASE}/analysis/scrna_integrated"
os.makedirs(OUT, exist_ok=True)

LOG = []


def log(msg):
    print(msg, flush=True)
    LOG.append(msg)


def to_symbols(a):
    """Ensembl 索引 + gene_symbol 列 → 符号索引（去重）"""
    if "gene_symbol" in a.var.columns:
        a.var["ensembl"] = a.var_names
        a.var_names = a.var["gene_symbol"].astype(str).values
    a.var_names_make_unique()
    return a


def basic_qc(a, name):
    n0, g0 = a.n_obs, a.n_vars
    sc.pp.filter_cells(a, min_genes=200)
    sc.pp.filter_genes(a, min_cells=3)
    # 线粒体
    mt = np.asarray(a.var_names.str.startswith(("mt-", "MT-")))
    if mt.sum() > 0:
        a.var["mt"] = mt
        sc.pp.calculate_qc_metrics(a, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
        a = a[a.obs["pct_counts_mt"] < 20].copy()
    else:
        a.obs["pct_mt"] = np.nan
    log(f"  [{name}] QC: cells {n0}->{a.n_obs}, genes {g0}->{a.n_vars}")
    return a


def ensure_counts(a):
    if a.X.dtype != np.float32:
        a.X = a.X.astype(np.float32)
    if sparse.issparse(a.X):
        a.X = a.X.tocsr()
    a.var_names_make_unique()
    return a


def load_h5ad(path, name):
    a = sc.read_h5ad(path)
    a.var["ensembl"] = a.var_names
    a = to_symbols(a)
    a.var_names_make_unique()
    log(f"  [{name}] loaded {a.n_obs} cells x {a.n_vars} genes; obs cols: {list(a.obs.columns)}")
    return a


def annotate(a, dataset, axis, condition, cell_type, cell_type_original=None, sample=None):
    a.obs["dataset"] = dataset
    a.obs["axis"] = axis
    a.obs["condition"] = condition
    a.obs["cell_type"] = cell_type
    a.obs["cell_type_original"] = cell_type_original if cell_type_original is not None else cell_type
    a.obs["sample"] = sample if sample is not None else a.obs.get("sample", pd.Series("NA", index=a.obs_names))
    return a


SV_MAP = {
    # Sun stria
    "Intermediate cell": "Intermediate_cell",
    "Basal cell": "Basal_cell",
    "Marginal cell": "Marginal_cell",
    "Capillary endothelial cell": "Capillary_endothelial",
    "PVM-like melanocyte": "Melanocyte_PVM",
    # Milon lateral wall
    "Intermediate Cells": "Intermediate_cell",
    "Marginal Cells": "Marginal_cell",
    "Basal Cells": "Basal_cell",
    "Fibrocytes": "Fibrocyte",
    "Spindle/Root Cells": "Spindle_root_cell",
    # Korrapati
    "Intermediate": "Intermediate_cell",
    "Marginal": "Marginal_cell",
    "Basal": "Basal_cell",
    "Fibrocyte": "Fibrocyte",
    "Spindle/Root": "Spindle_root_cell",
    "Unknown": "Unannotated",
    "Erythrocyte": "Other",
    "Erythroblast": "Other",
    "B-Cell": "Immune",
    "Macrophage": "Immune",
    "Neutrophil": "Immune",
}


def main():
    log("# scRNA 整合日志")
    log("")

    # ---------- 1. Sun 血管纹子集（衰老轴） ----------
    log("## 1. Sun 2023 stria (aging)")
    sun_sv = load_h5ad(f"{SRC}/sun2023_cra004814_stria.h5ad", "Sun_stria")
    sun_sv = annotate(
        sun_sv, "Sun2023_stria", "aging",
        condition=sun_sv.obs["age"].astype(str),
        cell_type=sun_sv.obs["CellType"].map(SV_MAP).fillna("Other"),
        cell_type_original=sun_sv.obs["CellType"].astype(str),
        sample=sun_sv.obs["sample"].astype(str),
    )
    sun_sv = basic_qc(sun_sv, "Sun_stria")

    # ---------- 2. Milon 侧壁（噪声轴） ----------
    log("## 2. Milon 2021 lateral wall (noise)")
    milon_sv = load_h5ad(f"{SRC}/milon2021_gse168041_lateralwall.h5ad", "Milon_LW")
    milon_sv = annotate(
        milon_sv, "Milon2021_LW", "noise",
        condition=milon_sv.obs["condition"].astype(str),
        cell_type=milon_sv.obs["celltype"].map(SV_MAP).fillna("Other"),
        cell_type_original=milon_sv.obs["celltype"].astype(str),
        sample=milon_sv.obs["replicate"].astype(str) if "replicate" in milon_sv.obs else None,
    )
    milon_sv = basic_qc(milon_sv, "Milon_LW")

    # ---------- 3. Korrapati naive SV（顺铂轴对照） ----------
    log("## 3. Korrapati 2019 naive SV (cisplatin control)")
    kor = load_h5ad(f"{SRC}/korrapati2019_naive_sv.h5ad", "Korrapati")
    kor = annotate(
        kor, "Korrapati2019", "ototoxic",
        condition="Control",
        cell_type=kor.obs["cell_type"].map(SV_MAP).fillna("Other"),
        cell_type_original=kor.obs["cell_type"].astype(str),
        sample=kor.obs["replicate"].astype(str) if "replicate" in kor.obs else None,
    )
    kor = basic_qc(kor, "Korrapati")

    # ---------- 4. Taukulis 顺铂 SV（无注释） ----------
    log("## 4. Taukulis 2021 cisplatin SV (unannotated)")
    tau = sc.read_10x_h5(f"{SRC}/GSE165662/GSM5047933_filtered_gene_bc_matrices_h5.h5")
    tau = ensure_counts(tau)
    tau = annotate(tau, "Taukulis2021", "ototoxic", condition="Cisplatin",
                   cell_type="Unannotated", cell_type_original="NA", sample="GSM5047933")
    tau = basic_qc(tau, "Taukulis")

    # ---------- 5. SV 合并 ----------
    log("## 5. Merge SV atlas")
    sv = ad.concat([sun_sv, milon_sv, kor, tau], join="outer", label="batch", index_unique="-")
    sv.obs_names_make_unique()
    log(f"  merged SV: {sv.n_obs} cells x {sv.n_vars} genes")
    sc.pp.filter_genes(sv, min_cells=3)
    log(f"  after gene filter: {sv.n_vars} genes")
    sv.write(f"{OUT}/sv_atlas.h5ad")
    log(f"  saved sv_atlas.h5ad")
    log("")
    log("  SV cell_type x dataset:")
    log("  " + pd.crosstab(sv.obs["cell_type"], sv.obs["dataset"]).to_string().replace("\n", "\n  "))
    log("")
    log("  SV condition x dataset:")
    log("  " + pd.crosstab(sv.obs["condition"], sv.obs["dataset"]).to_string().replace("\n", "\n  "))

    # ---------- 6. Sun 全图谱 + Miao GSE281324（全耳蜗） ----------
    log("")
    log("## 6. Whole-cochlea atlas")
    sun_full = load_h5ad(f"{SRC}/sun2023_cra004814_full_atlas.h5ad", "Sun_full")
    sun_full = annotate(
        sun_full, "Sun2023_full", "aging",
        condition=sun_full.obs["age"].astype(str),
        cell_type=sun_full.obs["CellType"].astype(str),
        cell_type_original=sun_full.obs["CellType"].astype(str),
        sample=sun_full.obs["sample"].astype(str),
    )
    sun_full = basic_qc(sun_full, "Sun_full")

    miao = []
    for gsm, cond in [("GSM8618786", "Control"), ("GSM8618787", "Control"),
                      ("GSM8618788", "Cisplatin"), ("GSM8618789", "Cisplatin")]:
        f = glob.glob(f"{SRC}/{gsm}_*.h5")[0]
        m = sc.read_10x_h5(f)
        m = ensure_counts(m)
        m = annotate(m, "Miao2024", "ototoxic", condition=cond,
                     cell_type="Unannotated", cell_type_original="NA", sample=gsm)
        m = basic_qc(m, gsm)
        miao.append(m)
    miao = ad.concat(miao, join="outer", index_unique="-")
    miao.obs_names_make_unique()
    log(f"  Miao merged: {miao.n_obs} cells x {miao.n_vars} genes")

    wc = ad.concat([sun_full, miao], join="outer", label="batch", index_unique="-")
    wc.obs_names_make_unique()
    log(f"  merged whole-cochlea: {wc.n_obs} cells x {wc.n_vars} genes")
    sc.pp.filter_genes(wc, min_cells=3)
    wc.write(f"{OUT}/whole_cochlea_atlas.h5ad")
    log("  saved whole_cochlea_atlas.h5ad")
    log("  " + pd.crosstab(wc.obs["condition"], wc.obs["dataset"]).to_string().replace("\n", "\n  "))

    # ---------- 7. SGN（Milon） ----------
    log("")
    log("## 7. SGN atlas")
    sgn = load_h5ad(f"{SRC}/milon2021_gse168041_sgn.h5ad", "Milon_SGN")
    sgn = annotate(
        sgn, "Milon2021_SGN", "noise",
        condition=sgn.obs["condition"].astype(str),
        cell_type=sgn.obs["celltype"].astype(str),
        cell_type_original=sgn.obs["celltype"].astype(str),
        sample=sgn.obs["replicate"].astype(str) if "replicate" in sgn.obs else None,
    )
    sgn = basic_qc(sgn, "Milon_SGN")
    sgn.write(f"{OUT}/sgn_atlas.h5ad")
    log("  saved sgn_atlas.h5ad")
    log("  " + pd.crosstab(sgn.obs["cell_type"], sgn.obs["condition"]).to_string().replace("\n", "\n  "))

    with open(f"{OUT}/integration_log.md", "w", encoding="utf-8") as fh:
        fh.write("\n".join(LOG) + "\n")
    log("")
    log("DONE")


if __name__ == "__main__":
    main()
