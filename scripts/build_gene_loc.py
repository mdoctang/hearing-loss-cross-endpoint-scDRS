#!/usr/bin/env python3
"""
build_gene_loc.py — 从 GENCODE v19 (GRCh37) GTF 构建 MAGMA gene-loc 文件

输出格式（制表符分隔，无表头，与 MAGMA NCBI37.3.hgnc.gene.loc 惯例一致）:
    GENE_SYMBOL  CHR  START  END  STRAND

说明:
- 官方 NCBI37.3.gene.loc（Entrez ID）因 ctg.cncr.nl 被沙箱网络封锁不可得，
  改用 GENCODE v19（GRCh37 坐标）构建等价的 HGNC 符号版，
  与 MAGMA.Celltyping 使用的 NCBI37.3.hgnc.gene.loc 惯例一致（符号作为基因 ID）。
- 仅保留 protein_coding 基因（与官方文件口径一致）。
- 基因边界 = 该基因所有转录本的最小 start / 最大 end。
- 基因符号重复时（GENCODE 中罕见，如 PAR 区重复）保留首个并记录。
"""
import gzip
import re
import sys
from collections import OrderedDict

GTF_PATH = sys.argv[1] if len(sys.argv) > 1 else "project/datasets/gene_annot/gencode.v19.annotation.gtf.gz"
OUT_PATH = sys.argv[2] if len(sys.argv) > 2 else "project/datasets/gene_annot/gencode.v19.hgnc.gene.loc"

genes = OrderedDict()  # symbol -> [chr, min_start, max_end, strand, gene_type]
dup_skipped = 0
n_total = 0

with gzip.open(GTF_PATH, "rt", encoding="utf-8") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 9 or parts[2] != "gene":
            continue
        chrom = parts[0].replace("chr", "")
        if chrom in ("M", "MT", "X", "Y"):  # MAGMA 主分析只用常染色体
            continue
        if not chrom.isdigit():
            continue
        start, end, strand = int(parts[3]), int(parts[4]), parts[6]
        attrs = parts[8]
        m_type = re.search(r'gene_type "([^"]+)"', attrs)
        m_name = re.search(r'gene_name "([^"]+)"', attrs)
        if not m_name:
            continue
        gtype = m_type.group(1) if m_type else ""
        if gtype != "protein_coding":
            continue
        symbol = m_name.group(1)
        n_total += 1
        if symbol in genes:
            # 同一符号出现在多个位置（罕见）：合并为覆盖区间并标记
            g = genes[symbol]
            if g[0] == chrom:
                g[1] = min(g[1], start)
                g[2] = max(g[2], end)
            else:
                dup_skipped += 1  # 跨染色体重复（如 PAR_X/Y）：保留首个
            continue
        genes[symbol] = [chrom, start, end, strand]

with open(OUT_PATH, "w", encoding="utf-8") as out:
    for symbol, (chrom, start, end, strand) in genes.items():
        out.write(f"{symbol}\t{chrom}\t{start}\t{end}\t{strand}\n")

print(f"protein_coding genes parsed: {n_total}")
print(f"unique symbols written:      {len(genes)}")
print(f"cross-chr duplicates skipped: {dup_skipped}")
print(f"output: {OUT_PATH}")
