#!/usr/bin/env python3
"""
run_scdrs_pipeline.py — scDRS 富集主分析（CAND-C v2 主要终点）

输入:
  --atlas  整合图谱 h5ad（.X 为原始计数；obs 含 cell_type）
  --gs     scDRS .gs 文件（人类基因；脚本自动转小鼠同源符号）
  --out    输出前缀目录
可选:
  --traits 只跑指定 trait（逗号分隔；默认全部）
  --n-ctrl 对照基因集数（默认 1000）
  --group-cols 分组列（默认 cell_type,dataset,condition）
  --adj-prop 传入 scdrs.pp.preprocess 的 adj_prop（默认 cell_type）

流程（scDRS 官方四步）:
  1) 归一化 + log1p；计算 n_genes/n_counts 协变量
  2) scdrs.pp.preprocess
  3) 每个 trait: scdrs.method.score_cell（weight_opt='vs', n_ctrl=1000）
  4) scdrs.method.downstream_group_analysis（细胞类型级 FDR）
输出:
  <out>/score_<trait>.csv.gz      每细胞疾病分数（含 pval）
  <out>/group_<trait>_<col>.csv   分组富集表
  <out>/pipeline_log.md
"""
import argparse
import os
import sys

# ---------------------------------------------------------------------------
# 必须在 import numpy/scipy 之前设置，否则不生效。
#
# 1) BLAS 线程数：本机 32 逻辑核，OpenBLAS 默认开满线程时在 Windows 上会随机
#    段错误（EXIT=139，无 Python traceback，日志停在 score_cell 中途）。
#    实测限制为 1 线程后完全稳定。如需提速可显式覆盖，例如
#    OPENBLAS_NUM_THREADS=8 python run_scdrs_pipeline.py ...
# 2) TQDM_DISABLE：关闭 scDRS 内部进度条，否则 1000 个对照基因集会刷出
#    上万行，把 run.log 撑爆并冲掉真正的 traceback。
# ---------------------------------------------------------------------------
for _v in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
):
    os.environ.setdefault(_v, "1")
os.environ.setdefault("TQDM_DISABLE", "1")

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import scdrs

# anndata>=0.11 默认拒绝写 nullable string 列，需显式表态。
# 置 False = 写 anndata<0.11 兼容的非 nullable 格式，最大化下游工具互通性。
try:
    ad.settings.allow_write_nullable_strings = False
except Exception:  # 旧版 anndata 无此设置，忽略
    pass

# ---------------------------------------------------------------------------
# NumPy 2.x 兼容垫片
# scDRS 1.0.2 的 scdrs.method.gearys_c 仍使用已被 NumPy 2.0 移除的 np.float_
# （method.py:1055 与 :1062），会在 downstream_group_analysis 阶段抛
# AttributeError。此处恢复别名，避免修改 site-packages、保证仓库可复现。
# ---------------------------------------------------------------------------
if not hasattr(np, "float_"):
    np.float_ = np.float64

# ---------------------------------------------------------------------------
# NumPy 2 / SciPy 兼容垫片 #2
# scdrs.method.distribution_match 把 rankdata(..., method="ordinal") 的结果直接
# 当数组下标用（method.py:971）；scipy 的 rankdata 返回 float64，而 NumPy 2 起
# 禁止浮点索引 → IndexError: arrays used as indices must be of integer type。
# ordinal 秩本身就是 1..n 的唯一整数，转 int64 语义完全等价。
# 该函数内部是 `from scipy.stats import rankdata`（调用时才导入），
# 因此补丁 scipy.stats.rankdata 即可生效，无需改 site-packages。
# ---------------------------------------------------------------------------
import scipy.stats as _sst

if not getattr(_sst.rankdata, "_scdrs_int_shim", False):
    _rankdata_orig = _sst.rankdata

    def _rankdata_shim(a, method="average", **kw):
        r = _rankdata_orig(a, method=method, **kw)
        return r.astype(np.int64) if method == "ordinal" else r

    _rankdata_shim._scdrs_int_shim = True
    _sst.rankdata = _rankdata_shim

LOG = []
def log(m):
    print(m, flush=True)
    LOG.append(m)

SCDRS_DATA = os.path.join(os.path.dirname(scdrs.__file__), "data")
HOMOLOG = os.path.join(SCDRS_DATA, "mouse_human_homologs.txt")
# 敏感性映射：MGI HOM_MouseHumanSequence（由 build_ortholog_map.py 构建）
# 见 study_design.md §6.3 —— 主映射仍是 scDRS 内置表，MGI 表只作稳健性检查。
MGI_MAP = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       "datasets", "ortholog", "mouse_to_human_symbol.tsv")


def load_homolog(source="builtin"):
    """
    返回 {HUMAN_GENE_SYM: MOUSE_GENE_SYM}。

    builtin  : scDRS 自带 mouse_human_homologs.txt（16,466 对官方鼠人同源）——**主映射**
    mgi      : MGI HOM_MouseHumanSequence 全量（24,584 对）——敏感性映射
    mgi_1to1 : MGI 中 is_1to1=True 的严格一对一子集（18,784 对）——敏感性映射
    """
    if source == "builtin":
        h = pd.read_csv(HOMOLOG, sep="\t")
        h.columns = [c.upper() for c in h.columns]
        h = h.dropna().drop_duplicates(subset=["MOUSE_GENE_SYM"])
        return dict(zip(h["HUMAN_GENE_SYM"], h["MOUSE_GENE_SYM"]))
    if not os.path.exists(MGI_MAP):
        raise SystemExit(f"找不到 MGI 映射表 {MGI_MAP}；先跑 scripts/build_ortholog_map.py")
    m = pd.read_csv(MGI_MAP, sep="\t")
    if source == "mgi_1to1":
        m = m[m["is_1to1"] == True]  # noqa: E712 —— 列内是 Python bool，需恒等比较
    m = m.dropna().drop_duplicates(subset=["mouse_symbol"])
    return dict(zip(m["human_symbol"].astype(str).str.upper(),
                    m["mouse_symbol"].astype(str)))


def parse_gs(gs_path):
    """读 .gs → {trait: [(gene, weight), ...]}"""
    out = {}
    with open(gs_path, encoding="utf-8") as fh:
        header = fh.readline()
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            trait, geneset = parts[0], parts[1]
            items = []
            for tok in geneset.split(","):
                if not tok:
                    continue
                g, w = tok.rsplit(":", 1)
                items.append((g, float(w)))
            out[trait] = items
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--atlas", required=True)
    ap.add_argument("--gs", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--traits", default=None)
    ap.add_argument("--n-ctrl", type=int, default=1000)
    ap.add_argument("--group-cols", default="cell_type_analysis,dataset,condition")
    ap.add_argument("--adj-prop", default="cell_type")
    ap.add_argument(
        "--weight-opt", default="vs",
        choices=["uniform", "vs", "inv_std", "od"],
        help="scDRS 权重方案；'vs' = 使用传入的 gene_weight（MAGMA zstat）",
    )
    ap.add_argument("--cache-adata", default=None,
                    help="预处理缓存路径（默认 <out>/preprocessed.h5ad）")
    ap.add_argument("--reuse", action="store_true",
                    help="若预处理缓存存在则直接复用（跳过归一化/PCA/preprocess，省时）")
    ap.add_argument("--no-cache", action="store_true", help="不写预处理缓存")
    ap.add_argument("--homolog", default="builtin",
                    choices=["builtin", "mgi", "mgi_1to1"],
                    help="跨物种同源映射表（默认 builtin = scDRS 内置，主分析用；"
                         "mgi / mgi_1to1 为 study_design.md §6.3 预登记的敏感性映射）")
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    log(f"# scDRS pipeline: {os.path.basename(args.atlas)}")

    cache_path = None if args.no_cache else (
        args.cache_adata or os.path.join(args.out, "preprocessed.h5ad")
    )

    # ---------- 1. 加载 + 归一化 + 协变量 ----------
    if args.reuse and cache_path and os.path.exists(cache_path):
        adata = sc.read_h5ad(cache_path)
        log(f"复用预处理缓存: {cache_path}（{adata.n_obs} 细胞 × {adata.n_vars} 基因）")
        # 还原被打包的 pandas.Series（见写缓存处的说明）
        _sp = adata.uns.get("SCDRS_PARAM", {})
        for _k, _v in _sp.items():
            if isinstance(_v, dict) and _v.get("__pd_series__"):
                _sp[_k] = pd.Series(_v["values"], index=_v["index"])
                log(f"uns/SCDRS_PARAM[{_k}] 已从打包结构还原为 Series")
    else:
        adata = sc.read_h5ad(args.atlas)
        log(f"atlas: {adata.n_obs} cells x {adata.n_vars} genes")
        # 优先使用 label transfer 后的标签（若存在）
        if "cell_type_transferred" in adata.obs.columns:
            adata.obs["cell_type_analysis"] = adata.obs["cell_type_transferred"].astype(str)
            log("使用 cell_type_transferred 作为分析标签")
        else:
            adata.obs["cell_type_analysis"] = adata.obs["cell_type"].astype(str)
        log(f"分析标签分布: {adata.obs['cell_type_analysis'].value_counts().to_dict()}")

        if "n_counts" not in adata.obs or adata.obs["n_counts"].isna().any():
            adata.obs["n_counts"] = np.asarray(adata.X.sum(1)).ravel()
        if "n_genes" not in adata.obs or adata.obs["n_genes"].isna().any():
            adata.obs["n_genes"] = np.asarray((adata.X > 0).sum(1)).ravel()

        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat")
        # 仅在 HVG 上做 PCA，避免全基因空间稠密化（大图谱内存安全）
        # scanpy>=1.11: mask_var 取代已弃用的 use_highly_variable
        sc.pp.pca(adata, n_comps=20, mask_var="highly_variable", svd_solver="arpack", random_state=0)
        sc.pp.neighbors(adata, n_neighbors=15, n_pcs=20, random_state=0)

    cov = pd.DataFrame(index=adata.obs_names)
    cov["const"] = 1.0
    cov["n_genes"] = adata.obs["n_genes"].astype(float).values
    cov["n_counts"] = adata.obs["n_counts"].astype(float).values
    if args.adj_prop == "cell_type":
        args.adj_prop = "cell_type_analysis"
    adj = args.adj_prop if args.adj_prop and args.adj_prop in adata.obs.columns else None
    if "SCDRS_PARAM" in adata.uns and "GENE_STATS" in adata.uns["SCDRS_PARAM"]:
        log("缓存中已有 SCDRS_PARAM，跳过 preprocess")
    else:
        scdrs.pp.preprocess(adata, cov=cov, adj_prop=adj)
        log(f"preprocess 完成（adj_prop={adj}）")
        if cache_path:
            # 原子写入：先写 .tmp 再替换，避免中途失败留下半截缓存被 --reuse 误读
            # anndata >= 0.11 拒绝把 pandas.Series 直接写入 uns：
            # sparse X 路径下 scdrs.pp.preprocess 会把 COV_GENE_MEAN 存为 Series，
            # 导致 IORegistryError: No method registered for writing Series into Group
            # （whole_cochlea_atlas 实测命中；dense 路径的 sgn_atlas 无此键不受影响）。
            # 解决：写缓存前打包成 dict（ndarray + index），写完立即还原，
            # 保证内存中的 adata 仍可被 scdrs 打分（COV_GENE_MEAN.loc 需要 Series）。
            _sp = adata.uns.get("SCDRS_PARAM", {})
            _packed = []
            for _k, _v in _sp.items():
                if isinstance(_v, pd.Series):
                    _sp[_k] = {"__pd_series__": True,
                               "values": np.asarray(_v.values),
                               "index": np.asarray(_v.index, dtype=object)}
                    _packed.append(_k)
            tmp = cache_path + ".tmp"
            adata.write_h5ad(tmp, compression="gzip")
            for _k in _packed:
                _d = _sp[_k]
                _sp[_k] = pd.Series(_d["values"], index=_d["index"])
            os.replace(tmp, cache_path)
            log(f"预处理缓存已写: {cache_path}（打包 Series 键: {_packed or '无'}）")

    # ---------- 2. 基因集：人类 → 小鼠同源 ----------
    human2mouse = load_homolog(args.homolog)
    log(f"同源映射表: {args.homolog}（{len(human2mouse)} 对）")
    gs_all = parse_gs(args.gs)
    traits = args.traits.split(",") if args.traits else list(gs_all.keys())
    log(f"gs traits: {list(gs_all.keys())}；运行: {traits}")

    # ---------- 3. 逐 trait 打分 ----------
    score_tables = {}
    for trait in traits:
        items = gs_all[trait]
        genes_m, weights_m = [], []
        for g, w in items:
            mg = human2mouse.get(g.upper())
            if mg is not None:
                genes_m.append(mg)
                weights_m.append(w)
        present = [g for g in genes_m if g in adata.var_names]
        w_map = dict(zip(genes_m, weights_m))
        genes_final = present
        weights_final = [w_map[g] for g in genes_final]
        log(f"[{trait}] 人类基因 {len(items)} → 小鼠同源 {len(genes_m)} → 图谱内 {len(genes_final)}")
        if len(genes_final) < 50:
            log(f"[{trait}] 有效基因 <50，跳过"); continue
        # 断点续跑：score csv 已存在则复用，只补缺失的 group 表
        # （2026-10-01 实测：长任务可能被超时/OOM 中断，50k/101k 细胞打分
        #  单次 15-25 分钟，重打浪费且再次暴露于中断风险）
        score_path = f"{args.out}/score_{trait}.csv.gz"
        df_score = None
        if os.path.exists(score_path):
            # ⚠️ 只判断 os.path.exists 不够（2026-10-01 事故）：后台任务被回收时会在
            # 磁盘上留下**只写了一半**的 csv.gz（实测 75 MB vs 完整 477 MB），
            # pd.read_csv 能读出一个行数严重不足的 DataFrame，且下游 group analysis
            # 与 compare 都不会报错，只会静默产出错误的分位数。
            # ⇒ 复用前必须校验：① 能读开；② 行数 == adata.n_obs；③ 索引与 obs_names 有交集。
            try:
                _df = pd.read_csv(score_path, index_col=0)
                _n = _df.shape[0]
                _ov = len(set(_df.index) & set(adata.obs_names))
                if _n == adata.n_obs and _ov > 0 and _df.shape[1] > 5:
                    df_score = _df
                    log(f"[{trait}] 复用已有打分结果: {score_path}（{_df.shape}）")
                else:
                    log(f"[{trait}] 已有打分文件不完整（行数 {_n} vs 细胞数 {adata.n_obs}、"
                        f"索引交集 {_ov}、列数 {_df.shape[1]}）→ 删除并重新打分")
                    os.remove(score_path)
            except Exception as e:
                log(f"[{trait}] 已有打分文件无法解析（{type(e).__name__}: {e}）→ 删除并重新打分")
                try:
                    os.remove(score_path)
                except OSError:
                    pass
        if df_score is None:
            df_score = scdrs.method.score_cell(
                adata, genes_final, gene_weight=weights_final,
                weight_opt=args.weight_opt, n_ctrl=args.n_ctrl,
                return_ctrl_norm_score=True, random_seed=0,
            )
            df_score.to_csv(score_path, compression="gzip")
            log(f"[{trait}] 打分完成: {df_score.shape}; FDR<0.1 细胞占比 {np.mean(df_score['pval'] < 0.1):.4f}")
        score_tables[trait] = df_score

        for gcol in args.group_cols.split(","):
            if gcol not in adata.obs.columns:
                continue
            gpath = f"{args.out}/group_{trait}_{gcol}.csv"
            if os.path.exists(gpath):
                log(f"[{trait}] group [{gcol}] 已存在，跳过")
                continue
            res = scdrs.method.downstream_group_analysis(adata, df_score, group_cols=[gcol])
            if gcol in res:
                res[gcol].to_csv(gpath)
                log(f"[{trait}] group analysis [{gcol}] 完成")

    with open(f"{args.out}/pipeline_log.md", "w", encoding="utf-8") as fh:
        fh.write("\n".join(LOG) + "\n")
    log("DONE")


if __name__ == "__main__":
    main()
