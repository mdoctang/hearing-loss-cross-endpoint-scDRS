#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
只读 QC：核查三个 FinnGen R13 端点的显著性 SNP 计数，用于 FUMA ERROR:005 故障定位。

不做任何写入/移动/删除，仅读取 .gz 文件并打印统计。

用法:
  python qc_significant_snps.py
"""

import gzip
import sys
import os
from collections import Counter

BASE = r"C:/Users/docta/WorkBuddy/纯生信/project"
RAW_DIR = os.path.join(BASE, "datasets", "finngen_r13")
FUMA_DIR = os.path.join(BASE, "analysis", "fuma_input_grch38")

ENDPOINTS = [
    ("H8_HL_SEN_NAS", "感音神经性听力损失"),
    ("H8_HL_CON_NAS", "传导性听力损失"),
    ("H8_NOISEINNER", "噪声性听力损失"),
]

RAW_FMT = "finngen_R13_{ep}.gz"
FUMA_FMT = "{ep}.grch38_3col.tsv.gz"

THRESHOLDS = [(5e-8, "5e-8"), (5e-6, "5e-6"), (1e-5, "1e-5")]

# 记录 P<5e-8 的 SNP 的最小 P 前 N 个
TOP_N = 20


def scan_raw(path):
    """流式扫描原始 FinnGen 13 列文件。返回统计字典。"""
    tot = 0
    cnt = {lab: 0 for _, lab in THRESHOLDS}
    maf_out = 0          # af_alt 在 [0.01,0.99] 之外
    maf_near = 0         # af_alt 在 [0.005,0.995] 之外（若也用 0.005 阈值可见差异）
    empty_rsid = 0       # rsids 列为空
    chrom_counter = Counter()
    top = []             # 最小 P 的 TOP_N
    sig_total = 0

    opener = gzip.open
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                # 表头行，跳过（第一列是 #chrom）
                if line.lstrip("#").startswith("chrom"):
                    continue
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 13:
                continue
            tot += 1
            chrom, pos, ref, alt, rsids, genes, pval_s, mlogp, beta, sebeta, af_s, afc, afk = parts[:13]
            try:
                pv = float(pval_s)
            except ValueError:
                continue
            for thr, lab in THRESHOLDS:
                if pv < thr:
                    cnt[lab] += 1
            if pv < 5e-8:
                sig_total += 1
                chrom_counter[chrom] += 1
                try:
                    af = float(af_s)
                except ValueError:
                    af = None
                if af is not None:
                    if not (0.01 <= af <= 0.99):
                        maf_out += 1
                    if not (0.005 <= af <= 0.995):
                        maf_near += 1
                if rsids.strip() in ("", ".", "NA"):
                    empty_rsid += 1
                top.append((pv, chrom, pos, rsids.strip(), af_s))

    top.sort(key=lambda x: x[0])
    return {
        "total": tot,
        "counts": cnt,
        "sig_total": sig_total,
        "maf_out": maf_out,
        "maf_near": maf_near,
        "empty_rsid": empty_rsid,
        "chrom": chrom_counter,
        "top": top[:TOP_N],
    }


def scan_fuma3col(path):
    """流式扫描 3 列 FUMA 上传文件 (chromosome, base_pair_location, p_value)。"""
    tot = 0
    cnt = {lab: 0 for _, lab in THRESHOLDS}
    chrom_counter = Counter()
    for line in _iter_gz(path):
        if line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 3:
            continue
        tot += 1
        try:
            pv = float(parts[2])
        except ValueError:
            continue
        for thr, lab in THRESHOLDS:
            if pv < thr:
                cnt[lab] += 1
        if pv < 5e-8:
            chrom_counter[parts[0]] += 1
    return {"total": tot, "counts": cnt, "chrom": chrom_counter}


def _iter_gz(path):
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            yield line


def main():
    summary = {}
    for ep, zh in ENDPOINTS:
        rp = os.path.join(RAW_DIR, RAW_FMT.format(ep=ep))
        fp = os.path.join(FUMA_DIR, FUMA_FMT.format(ep=ep))

        print("=" * 78)
        print(f"端点 {ep}  ({zh})")
        print("=" * 78)

        if not os.path.exists(rp):
            print(f"[缺失] 原始文件不存在: {rp}")
            continue

        print(f"[A] 原始 FinnGen 13 列文件: {rp}")
        r = scan_raw(rp)
        print(f"    A1 总数据行数(不含表头)          : {r['total']}")
        for _, lab in THRESHOLDS:
            print(f"    A2 P < {lab:<5} SNP 个数           : {r['counts'][lab]}")
        print(f"    A5 P<5e-8 中 af_alt 在 [0.01,0.99] 之外 : {r['maf_out']}"
              f"  (占比 {_pct(r['maf_out'], r['sig_total'])})")
        print(f"       P<5e-8 中 af_alt 在 [0.005,0.995] 之外(参考): {r['maf_near']}"
              f"  (占比 {_pct(r['maf_near'], r['sig_total'])})")
        print(f"    A6 P<5e-8 按染色体分布           : {dict(sorted(r['chrom'].items(), key=_chrom_key))}")
        print(f"    A7 P<5e-8 中 rsids 为空的数量     : {r['empty_rsid']}"
              f"  (占比 {_pct(r['empty_rsid'], r['sig_total'])})")
        print(f"    A8 P 值最小的前 {TOP_N} 个 SNP (chrom pos rsids pval af_alt):")
        for i, (pv, c, p, rs, af) in enumerate(r["top"], 1):
            print(f"       {i:>2}. {c:<3} {p:<12} {rs:<22} {pv:<12.4e} {af}")

        if os.path.exists(fp):
            print(f"[B] FUMA 上传 3 列文件: {fp}")
            f3 = scan_fuma3col(fp)
            print(f"    B9 总数据行数                    : {f3['total']}")
            for _, lab in THRESHOLDS:
                print(f"    B9 P < {lab:<5} SNP 个数           : {f3['counts'][lab]}")
            print(f"    B9.1 P<5e-8 按染色体分布          : {dict(sorted(f3['chrom'].items(), key=_chrom_key))}")
            for _, lab in THRESHOLDS:
                diff = f3["counts"][lab] - r["counts"][lab]
                print(f"    B10 P < {lab:<5} 差异(3列 - 原始)     : {diff:+d}")
        else:
            print(f"[B] [缺失] 3 列文件不存在: {fp}")

        summary[ep] = {"raw_sig": r["counts"]["5e-8"], "fuma_sig": None}
        print()

    print("=" * 78)
    print("[C] 跨端点比较: P < 5e-8 SNP 个数")
    print("=" * 78)
    for ep, zh in ENDPOINTS:
        s = summary.get(ep)
        print(f"    {ep:<16} {zh:<14} : {s['raw_sig'] if s else 'N/A'}")


def _pct(n, d):
    if not d:
        return "n/a"
    return f"{100.0 * n / d:.2f}%"


def _chrom_key(item):
    c = item[0]
    c2 = c.lower().replace("chr", "")
    if c2.isdigit():
        return (0, int(c2))
    return (1, c2)


if __name__ == "__main__":
    main()
