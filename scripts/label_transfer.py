#!/usr/bin/env python3
"""
label_transfer.py — 未注释细胞的小鼠耳蜗细胞类型迁移（kNN + 置信度门控）

方法（写进 Methods 的表述）:
  1) 参考集与查询集合并 → 共同高变基因（HVG, 2000）→ 标准化 + log1p → PCA(30)
  2) 以参考集已注释细胞训练 kNN（k=15，欧氏距离）
  3) 查询集细胞取近邻多数标签；置信度 = 近邻中该标签占比
  4) 置信度 ≥ 阈值（默认 0.6）才写入 cell_type_transferred，否则保留原标签（Unannotated）

输出: 就地更新 h5ad 的两个新列（cell_type_transferred, transfer_conf）并保存副本 *_labeled.h5ad
用法:
  python label_transfer.py <atlas.h5ad> <ref_datasets逗号分隔> <query_datasets逗号分隔> <out.h5ad>
"""
import sys
import numpy as np
import scanpy as sc
import pandas as pd
from sklearn.neighbors import KNeighborsClassifier

def main():
    atlas_path, ref_names, qry_names, out_path = sys.argv[1:5]
    conf_thr = float(sys.argv[5]) if len(sys.argv) > 5 else 0.6
    ref_names = ref_names.split(",")
    qry_names = qry_names.split(",")

    ad = sc.read_h5ad(atlas_path)
    print(f"loaded {ad.n_obs} x {ad.n_vars}", flush=True)
    ref_mask = ad.obs["dataset"].isin(ref_names).values
    qry_mask = ad.obs["dataset"].isin(qry_names).values
    print(f"ref cells {ref_mask.sum()}, query cells {qry_mask.sum()}", flush=True)

    sc.pp.normalize_total(ad, target_sum=1e4)
    sc.pp.log1p(ad)
    sc.pp.highly_variable_genes(ad, n_top_genes=2000, batch_key="dataset")
    hvg = ad.var["highly_variable"].values
    print(f"HVG: {hvg.sum()}", flush=True)

    # 关键：仅在 HVG 子集上稠密化 + 标准化 + PCA，避免全矩阵 densify 爆内存
    adh = ad[:, hvg].copy()
    adh.X = adh.X.astype(np.float32)
    sc.pp.scale(adh, max_value=10, zero_center=True)
    sc.tl.pca(adh, n_comps=30, svd_solver="arpack", random_state=0)
    print("PCA done", flush=True)

    ad.obsm["X_pca"] = adh.obsm["X_pca"]
    X = ad.obsm["X_pca"]
    idx_ref = np.where(ref_mask)[0]
    idx_qry = np.where(qry_mask)[0]
    y_ref = ad.obs["cell_type"].values[idx_ref]

    knn = KNeighborsClassifier(n_neighbors=15, weights="uniform", n_jobs=-1)
    knn.fit(X[idx_ref], y_ref)
    proba = knn.predict_proba(X[idx_qry])
    classes = knn.classes_
    best = proba.argmax(1)
    conf = proba[np.arange(len(idx_qry)), best]
    pred = classes[best]

    final = ad.obs["cell_type"].values.copy()
    assign_mask = conf >= conf_thr
    qry_labels = ad.obs["cell_type"].values[idx_qry].copy()
    qry_labels[assign_mask] = pred[assign_mask]
    final[idx_qry] = qry_labels

    ad.obs["cell_type_transferred"] = final
    ad.obs["transfer_conf"] = np.nan
    ad.obs.iloc[idx_qry, ad.obs.columns.get_loc("transfer_conf")] = conf

    print("转移后标签分布（查询集）:", flush=True)
    print(pd.Series(final[idx_qry]).value_counts().to_string(), flush=True)
    print(f"置信度 >= {conf_thr} 被赋标签: {assign_mask.sum()}/{len(idx_qry)}", flush=True)

    ad.write(out_path)
    print(f"saved {out_path}", flush=True)

if __name__ == "__main__":
    main()
