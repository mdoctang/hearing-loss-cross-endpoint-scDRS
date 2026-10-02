#!/usr/bin/env python3
"""
prep_gene_z_from_fuma.py — FUMA magma.genes.out × N 端点 → scDRS zscore_file 矩阵

输入: 多个 FUMA MAGMA 基因级结果文件（magma.genes.out，含 GENE 与 ZSTAT 列）
输出: zscore_file tsv（首列表头 'GENE'，其后每列一个端点 trait 名）
      可直接喂给 `scdrs munge-gs --zscore_file <out> --weight zscore --n-max 1000`

⚠ 物种桥接（本脚本的关键职责）
--------------------------------
FUMA/MAGMA 的基因符号是 **人类 HGNC**（如 A1BG、A2M），
而三个耳蜗单细胞图谱的 var 索引是 **小鼠 MGI**（如 A1bg、A2m、0610005C13Rik）。
直接匹配命中数 ≈ 10/18906（等于零），scDRS 会拿不到基因。
故本脚本默认把人类符号经 **MGI 官方同源表** 映射为小鼠符号后再输出。

映射规则：
  - 默认只用 is_1to1 == True 的严格 1:1 同源对（避免一对多歧义）
  - A1BG(人) → A1bg(鼠)，大小写按同源表原样（与图谱 var 索引一致的精确匹配）

用法:
  python prep_gene_z_from_fuma.py out.tsv EP1=path1.genes.out EP2=path2.genes.out ...
  python prep_gene_z_from_fuma.py out.tsv --no-ortholog EP1=...   # 关闭映射（人类符号直出）
  python prep_gene_z_from_fuma.py out.tsv --ortholog-table <tsv> EP1=...
"""
import os
import sys
import pandas as pd

# 默认同源表（相对本脚本定位，保证任意工作目录下可用）
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_DEFAULT_ORTHOLOG = os.path.normpath(
    os.path.join(_SCRIPT_DIR, "..", "datasets", "ortholog", "mouse_to_human_symbol.tsv")
)


def load_ortholog(path, one2one_only=True):
    """人类符号 → 小鼠符号。返回 dict；一对多时取排序后首个（可复现）。"""
    o = pd.read_csv(path, sep="\t")
    req = {"mouse_symbol", "human_symbol"}
    if not req <= set(o.columns):
        raise SystemExit(f"{path}: 缺列 {req - set(o.columns)}，实见 {list(o.columns)}")
    if one2one_only and "is_1to1" in o.columns:
        o = o[o["is_1to1"].astype(str).str.lower() == "true"]
    h2m = {}
    for h, g in o.groupby("human_symbol"):
        h2m[str(h)] = sorted(g["mouse_symbol"].astype(str))[0]
    return h2m


def main():
    args = sys.argv[1:]
    if not args:
        raise SystemExit(__doc__)

    out_path = args[0]
    rest = args[1:]

    use_ortholog = "--no-ortholog" not in rest
    ortholog_path = _DEFAULT_ORTHOLOG
    if "--ortholog-table" in rest:
        i = rest.index("--ortholog-table")
        ortholog_path = rest[i + 1]
        del rest[i:i + 2]
    rest = [a for a in rest if a != "--no-ortholog"]

    pairs = [a.split("=", 1) for a in rest]
    if not pairs:
        raise SystemExit("未提供任何 trait=path 输入")

    h2m = None
    if use_ortholog:
        if not os.path.exists(ortholog_path):
            raise SystemExit(f"同源表不存在: {ortholog_path}（可用 --no-ortholog 关闭映射）")
        h2m = load_ortholog(ortholog_path)
        print(f"[ortholog] {ortholog_path}")
        print(f"[ortholog] human->mouse 1:1 对数: {len(h2m)}")

    merged = None
    stats = {}
    for trait, path in pairs:
        df = pd.read_csv(path, sep=r"\s+")
        cols = {c.upper(): c for c in df.columns}
        gcol = cols.get("GENE")
        zcol = cols.get("ZSTAT")
        if gcol is None or zcol is None:
            raise SystemExit(f"{path}: 缺 GENE/ZSTAT 列，实见列: {list(df.columns)}")
        # MAGMA 优先用 SYMBOL 列做键（GENE 是 Ensembl ID）
        scol = cols.get("SYMBOL")
        keycol = scol if scol is not None else gcol
        s = df.set_index(keycol)[zcol]
        s.index = s.index.astype(str)
        s = s[~s.index.duplicated(keep="first")].rename(trait)

        n_before = len(s)
        if h2m is not None:
            # 人类符号 → 小鼠符号；无同源映射的基因丢弃
            mapped = s.index.map(lambda g: h2m.get(g))
            keep = mapped.notna()
            s = pd.Series(s.values[keep], index=mapped[keep], name=trait)
            s = s[~s.index.duplicated(keep="first")]
        stats[trait] = (n_before, len(s))
        merged = s.to_frame() if merged is None else merged.join(s, how="outer")

    merged.index.name = "GENE"          # scDRS munge-gs 要求首列表头为 GENE
    merged.to_csv(out_path, sep="\t", float_format="%.6g")

    print(f"traits: {list(merged.columns)}")
    print(f"genes(output): {merged.shape[0]}")
    for t, (b, a) in stats.items():
        print(f"  {t}: {b} -> {a} 基因（映射保留 {100*a/b:.1f}%）")
    print(f"per-trait non-null: {merged.notna().sum().to_dict()}")
    print(f"output: {out_path}")
    if use_ortholog:
        print("提示：输出键为【小鼠 MGI 符号】，与图谱 var 索引一致，可直接喂 scdrs munge-gs。")


if __name__ == "__main__":
    main()
