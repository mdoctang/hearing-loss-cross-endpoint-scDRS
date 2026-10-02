#!/usr/bin/env bash
# =============================================================================
# run_gs2000_noise_retry.sh — gs2000 的 NOISE 端点：带重试的打分 + 只在三端点齐全时才 compare
#
# 为什么需要重试：
#   Windows + OpenBLAS 在 101k 细胞规模上会**随机**段错误（SIGSEGV, rc=139，无 traceback），
#   run_main_scdrs.sh 顶部注释已记录。2026-10-01 实测：gs2000 的 SEN / CON 一次性成功，
#   而 NOISE 连续两次 segfault（一次跑到 436 MB 被回收、一次只撑 8 秒）。
#   ⇒ 单纯重跑有较高概率成功，故做最多 5 次尝试。
#
# 为什么必须校验端点齐全：
#   上一次失败时 bash -c 里的 `$?` 取值错位，日志打印 "scDRS rc=0"，
#   导致 compare_endpoints.py 拿**只有 CON 和 SEN 两个端点**的目录跑出了一份
#   28 项（而非 84 项）的产物，外表完全正常。⇒ compare 前必须硬性核对三个
#   score_*.csv.gz 都存在且 gzip -t 通过，否则拒绝 compare 并退出非 0。
# =============================================================================
set -uo pipefail

ROOT="C:/Users/docta/WorkBuddy/纯生信/project"
PY="C:/Users/docta/.workbuddy/binaries/python/envs/bioinfo/Scripts/python.exe"
ATLAS="$ROOT/analysis/scrna_integrated/whole_cochlea_atlas_labeled.h5ad"
CACHE="$ROOT/analysis/scdrs_main/whole_cochlea_atlas/preprocessed.h5ad"
OUT="$ROOT/analysis/sensitivity/gs2000/whole_cochlea_atlas"
GS="$ROOT/analysis/sensitivity/gs2000/hl_endpoints_gs2000.gs"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export TQDM_DISABLE=1
export NUMBA_OPT=0
export PYTHONUNBUFFERED=1

log() { printf '[%s] %s\n' "$(date '+%H:%M:%S')" "$*"; }

ok_file () {  # $1=path
  [ -f "$1" ] || return 1
  gzip -t "$1" 2>/dev/null || return 1
  return 0
}

for attempt in 1 2 3 4 5; do
  if ok_file "$OUT/score_NOISE.csv.gz"; then
    log "NOISE 已完整，跳过打分"
    break
  fi
  rm -f "$OUT/score_NOISE.csv.gz"
  log "第 $attempt 次尝试：NOISE 打分"
  "$PY" "$ROOT/scripts/run_scdrs_pipeline.py" \
    --atlas "$ATLAS" --gs "$GS" --out "$OUT" \
    --traits NOISE --n-ctrl 1000 --weight-opt vs \
    --homolog builtin --reuse --group-cols none \
    --cache-adata "$CACHE"
  rc=$?
  log "第 $attempt 次 rc=$rc"
  if [ $rc -eq 0 ] && ok_file "$OUT/score_NOISE.csv.gz"; then
    log "第 $attempt 次成功"
    break
  fi
  log "第 $attempt 次失败（rc=$rc / 文件不完整），30 秒后重试"
  rm -f "$OUT/score_NOISE.csv.gz"
  sleep 30
done

# ---------- compare 前的硬性门检 ----------
missing=0
for t in SEN CON NOISE; do
  if ok_file "$OUT/score_${t}.csv.gz"; then
    log "端点 $t ✓（$(stat -c%s "$OUT/score_${t}.csv.gz") bytes）"
  else
    log "端点 $t ✗ 缺失或不完整 → 拒绝 compare"
    missing=1
  fi
done
if [ "$missing" -ne 0 ]; then
  log "ABORT：端点不齐，compare 会产生错误的检验项数（2 端点 28 项 vs 3 端点 84 项）"
  exit 9
fi

log "compare（三端点齐全）"
"$PY" "$ROOT/scripts/compare_endpoints.py" \
  --atlas "$ATLAS" --dir "$OUT" --group-col cell_type_analysis
rc=$?
log "compare rc=$rc"
exit $rc
