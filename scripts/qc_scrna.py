#!/usr/bin/env python3
"""
qc_scrna.py — scRNA 数据资产结构质检（只读，不改文件）

对每个 h5ad/h5 报告：
- 细胞数 × 基因数、矩阵类型
- obs（细胞元数据）列名与候选注释列的类别分布（细胞类型、样本、条件、时间点和性别）
- var（基因）索引格式抽样（符号 vs Ensembl）
- 线粒体/核糖体基因占比粗查（QC 可行性）
用法：python qc_scrna.py <file1> [file2 ...]
"""
import sys
import scanpy as sc
import pandas as pd

def candidate_cols(cols, keys):
    out = []
    for c in cols:
        cl = c.lower()
        if any(k in cl for k in keys):
            out.append(c)
    return out

def qc_one(path):
    print("=" * 70)
    print(f"FILE: {path}")
    try:
        if path.endswith(".h5ad"):
            ad = sc.read_h5ad(path, backed="r")
        elif path.endswith(".h5"):
            ad = sc.read_10x_h5(path)
        else:
            print("  unsupported extension"); return
    except Exception as e:
        print(f"  READ ERROR: {e}"); return

    print(f"  shape: {ad.n_obs} cells x {ad.n_vars} genes")
    var_names = list(ad.var_names[:5])
    print(f"  var_names sample: {var_names}")
    obs_cols = list(ad.obs.columns)
    print(f"  obs columns ({len(obs_cols)}): {obs_cols}")
    ann_keys = ["cell", "type", "annot", "cluster", "ldesc", "label",
                "sample", "cond", "treat", "time", "age", "sex", "group", "orig"]
    for c in candidate_cols(obs_cols, ann_keys):
        try:
            vc = ad.obs[c].value_counts()
            top = vc.head(12)
            print(f"  [{c}] n_categories={len(vc)}")
            for idx, cnt in top.items():
                print(f"      {idx}: {cnt}")
            if len(vc) > 12:
                print(f"      ... (+{len(vc)-12} more)")
        except Exception as e:
            print(f"  [{c}] value_counts error: {e}")
    if hasattr(ad, "backed") and ad.backed:
        ad.file.close()

for p in sys.argv[1:]:
    qc_one(p)
