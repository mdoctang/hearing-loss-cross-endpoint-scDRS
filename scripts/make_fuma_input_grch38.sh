#!/usr/bin/env bash
# make_fuma_input_grch38.sh — FUMA GRCh38 合规上传文件生成（shell 版）
#
# 为什么另有一个 shell 版：
#   Python 版（make_fuma_input_grch38.py）逻辑正确，但本机在跑 19.5M 行的流式
#   gzip 读写时会**概率性 SIGILL 崩溃**（bash 报 "Illegal instruction"，无 traceback；
#   与 ENVIRONMENT.md §3.3 的 BLAS 段错误同属本机抖动类问题，但此处与 BLAS 无关）。
#   awk+gzip 管道全程无 Python 运行时、内存占用恒定，实测稳定，作为主用通道；
#   Python 版保留作可读的参考实现与列名映射文档。
#
# 输出（FUMA 勾选 GRCh38 时要求**恰好 3 列或恰好 6 列**，列名与顺序不可变）：
#   <trait>.grch38_3col.tsv.gz : chromosome  base_pair_location  p_value
#   <trait>.grch38_6col.tsv.gz : chromosome  base_pair_location  effect_allele  other_allele  beta  p_value
# 换行一律 LF（CRLF 会让表头末列变成 "p_value\r"，触发 FUMA ERROR:001）。
#
# 源文件列序（10 列）: SNP A1 A2 FRQ BETA SE P N CHR BP
#                    $1  $2 $3 $4  $5   $6 $7 $8 $9  $10
#
# 用法:
#   bash scripts/make_fuma_input_grch38.sh [SRC_DIR] [OUT_DIR] [TRAIT ...]
#   bash scripts/make_fuma_input_grch38.sh analysis/fuma_input analysis/fuma_input_grch38 H8_HL_CON_NAS
set -uo pipefail

SRC_DIR="${1:-analysis/fuma_input}"
OUT_DIR="${2:-analysis/fuma_input_grch38}"
shift 2 2>/dev/null || true
if [ "$#" -gt 0 ]; then TRAITS=("$@"); else TRAITS=(H8_HL_SEN_NAS H8_HL_CON_NAS H8_NOISEINNER); fi

mkdir -p "$OUT_DIR"

# awk 程序：BEGIN 打表头（OFS=\t），NR==1 跳过源表头，CHR 限 1-22 纯整数
read -r -d '' AWK3 <<'AWK' || true
BEGIN{OFS="\t"; print "chromosome","base_pair_location","p_value"}
NR==1{next}
{ c=$9+0; if ($9 ~ /^[0-9]+$/ && c>=1 && c<=22) print $9,$10,$7; n++ }
END{ printf "  src_rows=%d\n", n > "/dev/stderr" }
AWK

read -r -d '' AWK6 <<'AWK' || true
BEGIN{OFS="\t"; print "chromosome","base_pair_location","effect_allele","other_allele","beta","p_value"}
NR==1{next}
{ c=$9+0; if ($9 ~ /^[0-9]+$/ && c>=1 && c<=22) print $9,$10,$2,$3,$5,$7; n++ }
END{ printf "  src_rows=%d\n", n > "/dev/stderr" }
AWK

for t in "${TRAITS[@]}"; do
  src="$SRC_DIR/$t.fuma.tsv.gz"
  if [ ! -f "$src" ]; then echo "✗ 缺少源文件 $src"; continue; fi
  p3="$OUT_DIR/$t.grch38_3col.tsv.gz"
  p6="$OUT_DIR/$t.grch38_6col.tsv.gz"
  t3="$p3.tmp"; t6="$p6.tmp"

  echo "→ $t 生成 3col ..."
  gzip -dc "$src" | awk -F'\t' "$AWK3" | gzip -c > "$t3"

  echo "→ $t 生成 6col ..."
  gzip -dc "$src" | awk -F'\t' "$AWK6" | gzip -c > "$t6"

  # 只有两个流都完整关闭才落成正式文件名（原子替换）
  if gzip -t "$t3" 2>/dev/null && gzip -t "$t6" 2>/dev/null; then
    mv -f "$t3" "$p3"; mv -f "$t6" "$p6"
    printf "✓ %s  3col=%.1fMB  6col=%.1fMB\n" "$t" \
      "$(stat -c%s "$p3" | awk '{print $1/1048576}')" \
      "$(stat -c%s "$p6" | awk '{print $1/1048576}')"
  else
    echo "✗ $t 输出 gzip 校验失败"; rm -f "$t3" "$t6"
  fi
done

echo
echo "提交要点：勾选 GRCh38；Section 1 列名框全部留空；Section 2 填整数 N：SEN 480520 / CON 435524 / NOISE 474575。"
