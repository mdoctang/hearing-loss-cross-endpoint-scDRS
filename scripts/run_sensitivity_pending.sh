#!/usr/bin/env bash
# =============================================================================
# run_sensitivity_pending.sh — 补跑两项未完成的敏感性分析
#
#   A 段：换同源桥接表（scDRS builtin）在 whole_cochlea 全图谱上的 CON / NOISE
#         （SEN 已于 2026-10-01 完成并落盘，脚本会复用）
#   B 段：基因集规模 2000 在 whole_cochlea **全图谱**（101,082 细胞）上的三端点
#         （原计划改到 19,630 细胞子集上跑，已作废——非随机子集上重新打分会
#           因 scDRS 内部以当次纳入细胞为参照的分布匹配而严重失真）
#
# 设计要点：
#   - run_scdrs_pipeline.py 支持断点续跑（score_*.csv.gz 已存在则复用），
#     故长任务被回收后可安全重跑本脚本。
#   - --group-cols none：跳过 scDRS 自带的 group analysis（本敏感性只需细胞级
#     打分 + 自研 compare_endpoints.py 的跨端点检验），省机时。
# =============================================================================
set -uo pipefail

ROOT="C:/Users/docta/WorkBuddy/纯生信/project"
PY="C:/Users/docta/.workbuddy/binaries/python/envs/bioinfo/Scripts/python.exe"

ATLAS="$ROOT/analysis/scrna_integrated/whole_cochlea_atlas_labeled.h5ad"
SENS="$ROOT/analysis/sensitivity"

export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-1}"
export TQDM_DISABLE=1
export NUMBA_OPT="${NUMBA_OPT:-0}"

log() { printf '[%s] %s\n' "$(date '+%H:%M:%S')" "$*"; }

run_pipe () {  # $1=gs  $2=outdir  $3=traits  $4=homolog
  local gs="$1" outdir="$2" traits="$3" homolog="$4"
  mkdir -p "$outdir"
  log "scDRS: outdir=$outdir traits=$traits homolog=$homolog"
  "$PY" "$ROOT/scripts/run_scdrs_pipeline.py" \
    --atlas "$ATLAS" --gs "$gs" --out "$outdir" \
    --traits "$traits" --n-ctrl 1000 --weight-opt vs \
    --homolog "$homolog" --reuse --group-cols none
  local rc=$?
  log "scDRS rc=$rc"
  return $rc
}

run_cmp () {  # $1=outdir
  local outdir="$1"
  log "compare: $outdir"
  "$PY" "$ROOT/scripts/compare_endpoints.py" \
    --atlas "$ATLAS" --dir "$outdir" --group-col cell_type_analysis
  local rc=$?
  log "compare rc=$rc"
  return $rc
}

# ---------------- A 段：换桥接表 ----------------
log "===== A 段：同源桥接表敏感性（scDRS builtin） ====="
run_pipe "$SENS/ortholog_builtin/hl_endpoints_builtin.gs" \
         "$SENS/ortholog_builtin/whole_cochlea_atlas" "CON,NOISE" "builtin"
run_cmp  "$SENS/ortholog_builtin/whole_cochlea_atlas"

# ---------------- B 段：基因集规模 2000 ----------------
log "===== B 段：基因集规模 2000（全图谱） ====="
run_pipe "$SENS/gs2000/hl_endpoints_gs2000.gs" \
         "$SENS/gs2000/whole_cochlea_atlas" "SEN,CON,NOISE" "builtin"
run_cmp  "$SENS/gs2000/whole_cochlea_atlas"

log "ALL DONE"
