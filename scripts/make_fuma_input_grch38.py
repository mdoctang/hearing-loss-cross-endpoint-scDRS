#!/usr/bin/env python3
"""
make_fuma_input_grch38.py — 按 FUMA 官方 GRCh38 规定重制 SNP2GENE 上传文件

背景（为什么必须重做）
----------------------
实测本项目 GWAS 坐标**为 GRCh38**：多条染色体的最大 BP 超过 GRCh37 的染色体长度
（如 chr17 最大 83,232,370 > GRCh37 的 81,195,210；chr18/chr20/chr3/chr5 同理），
且各染色体末端坐标均紧贴 GRCh38 长度。

FUMA 官方 prepare_input_files 对 GRCh38 的硬性规定：
  - 若 sumstat 含 rsID：提交前**必须完全移除 chromosome 与 position**，让 FUMA 用
    dbSNP146 回查坐标（GRCh37 模式）；**或**
  - 勾选页面上的 GRCh38 选项，此时输入文件必须**恰好 3 列**：
        chromosome  base_pair_location  p_value
    或**恰好 6 列**：
        chromosome  base_pair_location  effect_allele  other_allele  beta  p_value
    （列名与顺序都必须完全一致；多一列即 ERROR:001）

我们此前的文件是 10 列、使用 GRCh37 命名（SNP/A1/A2/FRQ/BETA/SE/P/N/CHR/BP），
在 GRCh38 模式下表头不合规 → FUMA 报 ERROR:001（三个任务全部同一错误，与该诊断吻合）。

本脚本产出两种合规版本
----------------------
1) <trait>.grch38_3col.tsv.gz   —— 3 列，配「勾选 GRCh38」提交（**首选**，体积最小）
       chromosome  base_pair_location  p_value
2) <trait>.grch38_6col.tsv.gz   —— 6 列，配「勾选 GRCh38」提交（需保留等位基因/BETA 时用）
       chromosome  base_pair_location  effect_allele  other_allele  beta  p_value

样本量（N）不再作为列出现，改为在提交页 Section 2 的 sample size 处**填数值**：
   H8_HL_SEN_NAS = 480520
   H8_HL_CON_NAS = 435524
   H8_NOISEINNER = 474575

用法
----
python scripts/make_fuma_input_grch38.py analysis/fuma_input analysis/fuma_input_grch38
"""
import gzip
import os
import sys

# trait -> (文件名, 总样本量 N)
TRAITS = {
    "H8_HL_SEN_NAS": 480520,
    "H8_HL_CON_NAS": 435524,
    "H8_NOISEINNER": 474575,
}

# 源文件列名 -> 语义
SRC = {"SNP": "SNP", "A1": "EA", "A2": "NEA", "BETA": "BETA", "P": "P", "CHR": "CHR", "BP": "BP"}
OK_CHR = {str(i) for i in range(1, 23)}


def main():
    src_dir, out_dir = sys.argv[1], sys.argv[2]
    # 可选第 3 参数：逗号分隔的 trait 子集（用于断点重跑，避免重做已完成的大文件）
    only = None
    if len(sys.argv) > 3:
        only = {t.strip() for t in sys.argv[3].split(",") if t.strip()}
    os.makedirs(out_dir, exist_ok=True)

    for trait, n_total in TRAITS.items():
        if only and trait not in only:
            continue
        src = os.path.join(src_dir, f"{trait}.fuma.tsv.gz")
        if not os.path.exists(src):
            print(f"✗ 缺少源文件 {src}")
            continue

        p3 = os.path.join(out_dir, f"{trait}.grch38_3col.tsv.gz")
        p6 = os.path.join(out_dir, f"{trait}.grch38_6col.tsv.gz")
        t3, t6 = p3 + ".tmp", p6 + ".tmp"        # 原子写入：中断只留 .tmp，不会留下截断的成品
        n_in = n_written = n_badchr = 0

        # 关键：读写两侧都显式处理换行。
        # Win 上 gzip.open(..., "wt") 默认 newline=None → "\n" 被转成 os.linesep("\r\n")，
        # 会让**表头最后一列实际变成 "p_value\r"**，而 FUMA 的表头校验是精确匹配 → ERROR:001。
        # 故：读侧 newline=""（不翻译）+ 统一 rstrip("\r\n")；写侧强制 newline="\n"。
        with gzip.open(src, "rt", encoding="utf-8", newline="") as fi, \
             gzip.open(t3, "wt", encoding="utf-8", newline="\n", compresslevel=6) as f3, \
             gzip.open(t6, "wt", encoding="utf-8", newline="\n", compresslevel=6) as f6:
            header = fi.readline().rstrip("\r\n").split("\t")
            idx = {SRC[c]: i for i, c in enumerate(header) if c in SRC}
            need = {"CHR", "BP", "P", "EA", "NEA", "BETA"}
            if not need <= set(idx):
                raise SystemExit(f"{src}: 缺列，实见 {header}")

            f3.write("chromosome\tbase_pair_location\tp_value\n")
            f6.write("chromosome\tbase_pair_location\teffect_allele\tother_allele\tbeta\tp_value\n")

            for line in fi:
                f = line.rstrip("\r\n").split("\t")
                n_in += 1
                chr_ = f[idx["CHR"]]
                if chr_ not in OK_CHR:          # FUMA 只接受 1-23/X；其余会被丢弃
                    n_badchr += 1
                    continue
                bp, p = f[idx["BP"]], f[idx["P"]]
                f3.write(f"{chr_}\t{bp}\t{p}\n")
                f6.write(f"{chr_}\t{bp}\t{f[idx['EA']]}\t{f[idx['NEA']]}\t{f[idx['BETA']]}\t{p}\n")
                n_written += 1

        os.replace(t3, p3)      # 只有完整写完才落成正式文件名
        os.replace(t6, p6)

        print(f"✓ {trait}: 读入 {n_in:,} → 写出 {n_written:,}（非 1-22 剔除 {n_badchr:,}）；N={n_total:,}")
        for p in (p3, p6):
            print(f"    {p}  {os.path.getsize(p)/1e6:.1f} MB")

    print("\n提示：提交时在页面勾选 GRCh38，并在 Section 2 的 sample size 处填上面给出的 N 值（数值，不是列名）。")


if __name__ == "__main__":
    main()
