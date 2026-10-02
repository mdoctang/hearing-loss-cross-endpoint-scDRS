#!/usr/bin/env bash
# =============================================================================
# run_sensitivity_a.sh — 只跑 A 段：换同源桥接表（scDRS builtin）在 whole_cochlea
#                        全图谱上的 CON / NOISE + 跨端点 compare
#
# 为何单独拆出来：
#   2026-10-01 实测发现，旧的合并脚本（A 段 + gs2000 B 段）在后台被系统回收后，
#   留下一个**只写了一半**的 score_CON.csv.gz（75 MB vs 完整 477 MB）。
#   而 run_scdrs_pipeline.py 的续跑逻辑只检查 os.path.exists，会静默复用这个
#   截断文件 → 下游 compare 报「需要 ≥2 个端点，实见 ['SEN']」。
#   ⇒ 每次重跑前必须先确认/删除不完整的 score_*.csv.gz。本脚本自带该检查。
#
# 另外一个后台任务正在跑 gs2000（B 段），本脚本不得再碰 gs2000 目录，
# 否则两个进程同时写 score_SEN.csv.gz 会互相覆盖。
# =============================================================================
set -uo pipefail

ROOT="C:/Users/docta/WorkBuddy/纯生信/project"
PY="C:/Users/docta/.workbuddy/binaries/python/envs/bioinfo/Scripts/python.exe"

ATLAS="$ROOT/analysis/scrna_integrated/whole_cochlea_atlas_labeled.h5ad"
OUT="$ROOT/analysis/sensitivity/ortholog_builtin/whole_cochlea_atlas"
GS="$ROOT/analysis/sensitivity/ortholog_builtin/hl_endpoints_builtin.gs"

export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-1}"
export TQDM_DISABLE=1
export NUMBA_OPT="${NUMBA_OPT:-0}"
export PYTHONUNBUFFERED=1

log() { printf '[%s] %s\n' "$(date '+%H:%M:%S')" "$*"; }

# ---------- 0. 清理不完整的打分文件 ----------
# 完整的 score_*.csv.gz 约 477 MB；小于 400 MB 一律视为被中断的残file，删除重跑。
for f in "$OUT"/score_*.csv.gz; do
  [ -e "$f" ] || continue
  sz=$(stat -c%s "$f")
  mb=$((sz / 1048576))
  if [ "$sz" -lt 419430400 ]; then
    log "删除不完整打分文件: $(basename "$f")（${mb} MB < 400 MB）"
    rm -f "$f"
  else
    log "保留完整打分文件: $(basename "$f")（${mb} MB）"
  fi
done

# ---------- 1. 打分 ----------
mkdir -p "$OUT"
CACHE="$OUT/preprocessed.h5ad"
CACHE_ARG=()
[ -f "$CACHE" ] && CACHE_ARG=(--cache-adata "$CACHE")

log "scDRS: CON,NOISE（homolog=builtin）"
"$PY" "$ROOT/scripts/run_scdrs_pipeline.py" \
  --atlas "$ATLAS" --gs "$GS" --out "$OUT" \
  --traits CON,NOISE --n-ctrl 1000 --weight-opt vs \
  --homolog builtin --reuse --group-cols none "${CACHE_ARG[@]}"
log "scDRS rc=$?"

# ---------- 2. 跨端点对比 ----------
log "compare"
"$PY" "$ROOT/scripts/compare_endpoints.py" \
  --atlas "$ATLAS" --dir "$OUT" --group-col cell_type_analysis
log "compare rc=$?"

log "A 段 DONE"
