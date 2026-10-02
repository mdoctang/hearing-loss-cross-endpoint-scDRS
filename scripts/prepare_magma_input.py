#!/usr/bin/env python3
"""
prepare_magma_input.py — FinnGen R13 summary stats → MAGMA 输入文件

对每个端点输出两个文件（制表符分隔）:
  {EP}.snploc : SNP  CHR  POS           （无表头，供 --snp-loc）
  {EP}.pval   : SNP  P    N             （带表头，供 --pval，N=总样本量）

QC 规则（与预登记统计计划一致）:
1. 仅保留常染色体 1-22（参考面板为常染色体）。
2. 仅保留带 rsID 的变异（rsids 列取第一个 rs* 条目）；无 rsID 的indel/结构变异剔除并计数。
3. pval 缺失/非有限值/超出 (0,1] 的剔除。
4. 同一 rsID 重复时保留 pval 最小者。
5. 输出按染色体、位置排序。

用法:
  python prepare_magma_input.py <finngen.gz> <EP> <N_total> <outdir>
例:
  python prepare_magma_input.py finngen_R13_H8_HL_SEN_NAS.gz H8_HL_SEN_NAS 409661 out
"""
import gzip
import sys
import csv

def main():
    src, ep, n_total, outdir = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4]
    snploc_path = f"{outdir}/{ep}.snploc"
    pval_path = f"{outdir}/{ep}.pval"

    n_in = n_kept = n_no_rsid = n_bad_p = n_dup = n_nonauto = 0
    best = {}  # rsid -> (p, chrom, pos)

    with gzip.open(src, "rt", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        cols = reader.fieldnames
        # FinnGen R13 列名: #chrom,pos,ref,alt,rsids,nearest_genes,pval,mlogp,beta,sebeta,af_alt,...
        c_chr = "#chrom" if "#chrom" in cols else "chrom"
        for row in reader:
            n_in += 1
            chrom = row[c_chr].replace("chr", "").replace("23", "X")
            if not chrom.isdigit():
                n_nonauto += 1
                continue
            rsids = (row.get("rsids") or "").strip()
            rsid = ""
            for tok in rsids.replace(";", ",").split(","):
                tok = tok.strip()
                if tok.startswith("rs"):
                    rsid = tok
                    break
            if not rsid:
                n_no_rsid += 1
                continue
            try:
                p = float(row["pval"])
            except (ValueError, TypeError):
                n_bad_p += 1
                continue
            if not (0.0 < p <= 1.0):
                n_bad_p += 1
                continue
            pos = int(row["pos"])
            if rsid in best:
                n_dup += 1
                if p < best[rsid][0]:
                    best[rsid] = (p, int(chrom), pos)
            else:
                best[rsid] = (p, int(chrom), pos)

    rows = sorted(best.items(), key=lambda kv: (kv[1][1], kv[1][2]))
    with open(snploc_path, "w", encoding="utf-8") as f1, open(pval_path, "w", encoding="utf-8") as f2:
        f2.write("SNP\tP\tN\n")
        for rsid, (p, chrom, pos) in rows:
            f1.write(f"{rsid}\t{chrom}\t{pos}\n")
            f2.write(f"{rsid}\t{p:.6e}\t{n_total}\n")
            n_kept += 1

    print(f"[{ep}] input rows: {n_in}")
    print(f"[{ep}] kept rsID SNPs: {n_kept}")
    print(f"[{ep}] dropped: no_rsid={n_no_rsid}, bad_p={n_bad_p}, dup={n_dup}, non_autosome={n_nonauto}")
    print(f"[{ep}] N_total={n_total}")
    print(f"[{ep}] outputs: {snploc_path} , {pval_path}")

if __name__ == "__main__":
    main()
