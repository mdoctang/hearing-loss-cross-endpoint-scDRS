#!/usr/bin/env bash
# =============================================================================
# run_sensitivity_gs2000.sh — 基因集规模 2000（whole_cochlea 全图谱）补跑 + compare
#
# 背景（2026-10-01）：
#   旧后台任务在写 score_NOISE.csv.gz 时被系统回收，留下 436 MB 的截断文件
#   （完整约 477 MB）。compare_endpoints.py 读它直接 EOFError：
#     "Compressed file ended before the end-of-stream marker was reached"
#   score_SEN / score_CON 是完整的，可复用；只需重打 NOISE。
#
#   已在 run_scdrs_pipeline.py 里加了完整性校验（行数 == adata.n_obs 等），
#   本脚本直接依赖该校验来判定该删哪个文件，不再用体积硬阈值。
#   另加一道：用 gzip -t 预检，坏档直接删，省得让 pandas 去撞 EOFError。
# =============================================================================
set -uo pipefail

ROOT="C:/Users/docta/WorkBuddy/纯生信/project"
PY="C:/Users/docta/.workbuddy/binaries/python/envs/bioinfo/Scripts/python.exe"

ATLAS="$ROOT/analysis/scrna_integrated/whole_cochlea_atlas_labeled.h5ad"
CACHE="$ROOT/analysis/scdrs_main/whole_cochlea_atlas/preprocessed.h5ad"
OUT="$ROOT/analysis/sensitivity/gs2000/whole_cochlea_atlas"
GS="$ROOT/analysis/sensitivity/gs2000/hl_endpoints_gs2000.gs"

export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-1}"
export TQDM_DISABLE=1
export NUMBA_OPT="${NUMBA_OPT:-0}"
export PYTHONUNBUFFERED=1

log() { printf '[%s] %s\n' "$(date '+%H:%M:%S')" "$*"; }

# ---------- 0. gzip 完整性预检 ----------
for f in "$OUT"/score_*.csv.gz; do
  [ -e "$f" ] || continue
  if gzip -t "$f" 2>/dev/null; then
    log "gzip 完整: $(basename "$f")（$(stat -c%s "$f") bytes）"
  else
    log "gzip 截断 → 删除重打: $(basename "$f")（$(stat -c%s "$f") bytes）"
    rm -f "$f"
  fi
done

# ---------- 1. 打分（SEN/CON 复用，NOISE 重打） ----------
mkdir -p "$OUT"
log "scDRS: SEN,CON,NOISE（gs=2000）"
"$PY" "$ROOT/scripts/run_scdrs_pipeline.py" \
  --atlas "$ATLAS" --gs "$GS" --out "$OUT" \
  --traits SEN,CON,NOISE --n-ctrl 1000 --weight-opt vs \
  --homolog builtin --reuse --group-cols none \
  --cache-adata "$CACHE"
log "scDRS rc=$?"

# ---------- 2. 跨端点对比 ----------
log "compare"
"$PY" "$ROOT/scripts/compare_endpoints.py" \
  --atlas "$ATLAS" --dir "$OUT" --group-col cell_type_analysis
log "compare rc=$?"

log "gs2000 DONE"
