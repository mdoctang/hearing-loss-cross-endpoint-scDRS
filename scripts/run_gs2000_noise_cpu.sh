#!/usr/bin/env bash
# =============================================================================
# run_gs2000_noise_cpu.sh — gs2000 NOISE 端点：切换 numba 代码生成策略后重试
#
# 失败史（2026-10-01）：
#   5 次重试全部失败，退出码两种：
#     rc=132 = 128+4 = **SIGILL（非法指令）**——7 秒内必崩
#     rc=139 = 128+11 = SIGSEGV——撑到 8–12 分钟后崩
#   SEN / CON 同规模一次跑通，NOISE 必崩 ⇒ 不是规模问题，也不是纯随机。
#   SIGILL 早崩 + numba 参与打分 ⇒ 典型症状是 numba/LLVM 按本机 CPU 特性生成了
#   本机不支持（或 Windows 沙箱下不可用）的指令。
#
# 对策（分两组，逐组重试）：
#   组 1（第 1–3 次）：NUMBA_CPU_NAME=generic + NUMBA_CPU_FEATURES=  → 强制通用代码生成
#   组 2（第 4–6 次）：NUMBA_DISABLE_JIT=1                          → 退回纯 Python（慢但最稳）
#
# compare 前仍执行三端点门检（见 run_gs2000_noise_retry.sh 的说明）。
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
ok_file () { [ -f "$1" ] || return 1; gzip -t "$1" 2>/dev/null || return 1; return 0; }

for attempt in 1 2 3 4 5 6; do
  if ok_file "$OUT/score_NOISE.csv.gz"; then log "NOISE 已完整"; break; fi
  rm -f "$OUT/score_NOISE.csv.gz"

  if [ "$attempt" -le 3 ]; then
    export NUMBA_CPU_NAME=generic
    export NUMBA_CPU_FEATURES=
    unset NUMBA_DISABLE_JIT 2>/dev/null || true
    mode="numba generic CPU"
  else
    export NUMBA_DISABLE_JIT=1
    unset NUMBA_CPU_NAME NUMBA_CPU_FEATURES 2>/dev/null || true
    mode="numba JIT 关闭（纯 Python，慢）"
  fi

  log "第 $attempt 次：模式 = $mode"
  "$PY" "$ROOT/scripts/run_scdrs_pipeline.py" \
    --atlas "$ATLAS" --gs "$GS" --out "$OUT" \
    --traits NOISE --n-ctrl 1000 --weight-opt vs \
    --homolog builtin --reuse --group-cols none \
    --cache-adata "$CACHE"
  rc=$?
  log "第 $attempt 次 rc=$rc"
  if [ $rc -eq 0 ] && ok_file "$OUT/score_NOISE.csv.gz"; then log "第 $attempt 次成功"; break; fi
  rm -f "$OUT/score_NOISE.csv.gz"
  sleep 20
done

missing=0
for t in SEN CON NOISE; do
  if ok_file "$OUT/score_${t}.csv.gz"; then log "端点 $t ✓"; else log "端点 $t ✗"; missing=1; fi
done
if [ "$missing" -ne 0 ]; then
  log "ABORT：端点不齐，拒绝 compare（会产生 28 项而非 84 项的错误产物）"
  exit 9
fi

log "compare（三端点齐全）"
"$PY" "$ROOT/scripts/compare_endpoints.py" \
  --atlas "$ATLAS" --dir "$OUT" --group-col cell_type_analysis
rc=$?
log "compare rc=$rc"
exit $rc
