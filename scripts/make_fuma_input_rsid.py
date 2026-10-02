#!/usr/bin/env python3
"""
make_fuma_input_rsid.py — 生成 FUMA SNP2GENE「纯 rsID 模式」上传文件（免勾选 GRCh38）

为什么需要它
------------
FUMA 官方 prepare_input_files 对 **GRCh38 且含 rsID** 的 sumstat 给出的解法是：
    "if your input GWAS sumstat file contains rsID, you can submit to FUMA by first
     removing the chromosome and position before submitting to FUMA.
     ... FUMA will use the provided rsID to look up the chromosomes and positions
     using dbSNP version 146"
    "IMPORTANT: make sure to completely remove the chromosome and position (in GRCh38)
     if present in your gwas sumstat."

这条路径的优势（相对 GRCh38 3col/6col 模式）：
  * 不依赖提交页面的 GRCh38 勾选框 → 消除「勾没勾」这一最大不确定因素；
  * 文件只要 4 列，不触发「必须恰好 3 或 6 列 / 多一列即 ERROR:001」的严格表头校验；
  * 列名 SNP / A1 / A2 / P 全部命中 FUMA 的自动列名识别表（snp→rsID, A1→effect,
    A2→non-effect, P→p），因此 **Section 1 的列名框可保持留空**；
  * FUMA 以 dbSNP146（GRCh37）回查坐标，输出天然是 GRCh37，与本项目
    gene.loc（GENCODE v19 = GRCh37）与 1000G EUR 面板（GRCh37）一致。

代价
----
源文件 `export_fuma_input.py` 已把「无 rsID / 非 1-22 / 重复 rsID」的行剔除，
因此本文件不含那 4.19% 无 rsID 的变异（约占原始行数 4.19%）。这些变异在 FUMA
参考面板中通常也无法匹配，属于常规取舍。

输入 → 输出
-----------
输入: analysis/fuma_input/<EP>.fuma.tsv.gz   （10 列: SNP A1 A2 FRQ BETA SE P N CHR BP）
输出: analysis/fuma_input_rsid/<EP>.fuma_rsid.tsv.gz （4 列: SNP A1 A2 P）

用法
----
python scripts/make_fuma_input_rsid.py analysis/fuma_input analysis/fuma_input_rsid
"""
import gzip
import os
import sys

EPISODES = ["H8_HL_SEN_NAS", "H8_HL_CON_NAS", "H8_NOISEINNER"]


def main():
    src_dir, out_dir = sys.argv[1], sys.argv[2]
    os.makedirs(out_dir, exist_ok=True)

    for ep in EPISODES:
        src = os.path.join(src_dir, f"{ep}.fuma.tsv.gz")
        if not os.path.exists(src):
            print(f"x 缺少源文件 {src}")
            continue
        out = os.path.join(out_dir, f"{ep}.fuma_rsid.tsv.gz")
        tmp = out + ".tmp"                     # 原子写入

        n_in = n_out = n_bad = 0
        # 读写两侧显式声明换行：Win 上 gzip.open("wt") 默认会把 "\n" 变成 "\r\n"，
        # 导致表头末列实际是 "P\r"，精确匹配的表头校验会失败。
        with gzip.open(src, "rt", encoding="utf-8", newline="") as fi, \
             gzip.open(tmp, "wt", encoding="utf-8", newline="\n", compresslevel=6) as fo:
            header = fi.readline().rstrip("\r\n").split("\t")
            idx = {c: i for i, c in enumerate(header)}
            need = {"SNP", "A1", "A2", "P"}
            if not need <= set(idx):
                raise SystemExit(f"{src}: 缺列，实见 {header}")

            fo.write("SNP\tA1\tA2\tP\n")
            for line in fi:
                f = line.rstrip("\r\n").split("\t")
                n_in += 1
                snp, a1, a2, p = f[idx["SNP"]], f[idx["A1"]], f[idx["A2"]], f[idx["P"]]
                if not snp.startswith("rs"):
                    n_bad += 1
                    continue
                fo.write(f"{snp}\t{a1}\t{a2}\t{p}\n")
                n_out += 1

        os.replace(tmp, out)
        print(f"v {ep}: 读入 {n_in:,} -> 写出 {n_out:,}（剔除 {n_bad:,}）；"
              f"{out}  {os.path.getsize(out)/1e6:.1f} MB")

    print("\n提交要点：默认模式（**不要**勾选 GRCh38）；Section 1 列名框留空；"
          "Section 2 sample size 填整数 N（SEN 480520 / CON 435524 / NOISE 474575）。")


if __name__ == "__main__":
    main()
