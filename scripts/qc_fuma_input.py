#!/usr/bin/env python3
"""
qc_fuma_input.py — FUMA SNP2GENE 上传文件体检

针对 FUMA 官方 ERROR:001（输入格式不正确）的已知成因逐项排查
（依据：FUMA GWAS users 组 troubleshooting 清单 + 官方 prepare_input_files 文档）：
  1. 表头列名与列数是否与数据行一致
  2. 是否含引号字符
  3. rsID 是否形如 rs<数字>
  4. 染色体是否在 1-23 / X / Y / MT 内
  5. **BP 是否整数、是否出现科学计数法**
  6. **P / BETA / SE / FRQ 是否科学计数法**（FUMA 对科学计数法敏感）
  7. 缺失值 / 非有限值 / FRQ 越界 / P<=0
  8. 等位基因是否为合法字符
  9. 重复 rsID
输出：逐项计数 + 每个问题的样例行（前 5 条）
"""
import gzip
import re
import sys
from collections import Counter

SCIRE = re.compile(r"^[+-]?(\d+\.?\d*|\.\d+)[eE][+-]?\d+$")
RSID = re.compile(r"^rs\d+$")
INT = re.compile(r"^-?\d+$")
OK_CHR = {str(i) for i in range(1, 23)} | {"X", "Y", "MT", "23", "24", "25", "26", "XY"}
ALLELE = re.compile(r"^[ACGTNacgtn]+$")


def scan(path, max_rows=None):
    probs = Counter()
    samples = {}
    n = 0
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        ncol = len(header)
        print(f"  表头({ncol}): {header}")
        for line in fh:
            n += 1
            if max_rows and n > max_rows:
                break
            if '"' in line or "'" in line:
                probs["含引号"] += 1
                samples.setdefault("含引号", line)
            f = line.rstrip("\n").split("\t")
            if len(f) != ncol:
                probs["列数不一致"] += 1
                samples.setdefault("列数不一致", line)
                continue
            d = dict(zip(header, f))
            snp = d.get("SNP", "")
            if not RSID.match(snp):
                probs["rsID格式非法"] += 1
                samples.setdefault("rsID格式非法", line)
            chr_ = d.get("CHR", "")
            if chr_.upper() not in OK_CHR:
                probs["CHR非法"] += 1
                samples.setdefault("CHR非法", line)
            bp = d.get("BP", "")
            if SCIRE.match(bp):
                probs["BP科学计数法"] += 1
                samples.setdefault("BP科学计数法", line)
            elif not INT.match(bp):
                probs["BP非整数"] += 1
                samples.setdefault("BP非整数", line)
            for col in ("P", "BETA", "SE", "FRQ", "N"):
                v = d.get(col, "")
                if v == "" or v.upper() in ("NA", "NAN", "."):
                    probs[f"{col}缺失"] += 1
                    samples.setdefault(f"{col}缺失", line)
                elif SCIRE.match(v):
                    probs[f"{col}科学计数法"] += 1
                    samples.setdefault(f"{col}科学计数法", line)
            try:
                p = float(d.get("P", "nan"))
                if not (0 < p <= 1):
                    probs["P越界(<=0或>1)"] += 1
                    samples.setdefault("P越界(<=0或>1)", line)
            except ValueError:
                probs["P不可解析"] += 1
                samples.setdefault("P不可解析", line)
            try:
                q = float(d.get("FRQ", "nan"))
                if not (0 <= q <= 1):
                    probs["FRQ越界"] += 1
                    samples.setdefault("FRQ越界", line)
            except ValueError:
                pass
            for col in ("A1", "A2"):
                if not ALLELE.match(d.get(col, "")):
                    probs[f"{col}非法字符"] += 1
                    samples.setdefault(f"{col}非法字符", line)
    return n, probs, samples


def main():
    rc = 0
    for path in sys.argv[1:]:
        print(f"\n=== {path} ===")
        n, probs, samples = scan(path)
        print(f"  扫描行数: {n:,}")
        if not probs:
            print("  ✓ 未发现已知 ERROR:001 成因")
        else:
            print("  发现问题：")
            for k, v in probs.most_common():
                print(f"    - {k}: {v:,} 行")
                print(f"      样例: {samples[k][:180].rstrip()}")
                rc = 1
    return rc


if __name__ == "__main__":
    sys.exit(main())
