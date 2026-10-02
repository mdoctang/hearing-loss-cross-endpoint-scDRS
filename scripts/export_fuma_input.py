#!/usr/bin/env python3
"""
export_fuma_input.py — FinnGen R13 summary stats → FUMA SNP2GENE 上传格式

输出列（制表符分隔，gzip）:
  SNP  A1  A2  FRQ  BETA  SE  P  N  CHR  BP

约定:
- A1 = 效应等位基因（FinnGen 的 alt，beta 为其效应）
- A2 = 非效应等位基因（FinnGen 的 ref）
- FRQ = af_alt（效应等位基因频率，全样本）
- 仅保留带 rsID 的常染色体变异（与 MAGMA 预处理口径一致，FUMA 以 rsID 映射 dbSNP）
- 同一 rsID 重复时保留 pval 最小者

用法:
  python export_fuma_input.py <finngen.gz> <EP> <N_total> <outdir>
"""
import gzip
import sys
import csv

def main():
    src, ep, n_total, outdir = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4]
    out_path = f"{outdir}/{ep}.fuma.tsv.gz"

    n_in = n_kept = n_no_rsid = n_bad = n_dup = n_nonauto = 0
    best = {}  # rsid -> (p, chrom, pos, a1, a2, frq, beta, se)

    with gzip.open(src, "rt", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        cols = reader.fieldnames
        c_chr = "#chrom" if "#chrom" in cols else "chrom"
        for row in reader:
            n_in += 1
            chrom = row[c_chr].replace("chr", "")
            if chrom == "23":
                chrom = "X"
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
                beta = float(row["beta"])
                se = float(row["sebeta"])
                frq = float(row["af_alt"])
            except (ValueError, TypeError):
                n_bad += 1
                continue
            if not (0.0 < p <= 1.0):
                n_bad += 1
                continue
            a2 = row["ref"].upper()
            a1 = row["alt"].upper()
            pos = int(row["pos"])
            rec = (p, int(chrom), pos, a1, a2, frq, beta, se)
            if rsid in best:
                n_dup += 1
                if p < best[rsid][0]:
                    best[rsid] = rec
            else:
                best[rsid] = rec

    rows = sorted(best.items(), key=lambda kv: (kv[1][1], kv[1][2]))
    with gzip.open(out_path, "wt", encoding="utf-8", newline="\n") as out:
        out.write("SNP\tA1\tA2\tFRQ\tBETA\tSE\tP\tN\tCHR\tBP\n")
        for rsid, (p, chrom, pos, a1, a2, frq, beta, se) in rows:
            out.write(f"{rsid}\t{a1}\t{a2}\t{frq:.6g}\t{beta:.6g}\t{se:.6g}\t{p:.6e}\t{n_total}\t{chrom}\t{pos}\n")
            n_kept += 1

    print(f"[{ep}] input rows: {n_in}")
    print(f"[{ep}] kept rsID SNPs: {n_kept}")
    print(f"[{ep}] dropped: no_rsid={n_no_rsid}, bad={n_bad}, dup={n_dup}, non_autosome={n_nonauto}")
    print(f"[{ep}] output: {out_path}")

if __name__ == "__main__":
    main()
