#!/usr/bin/env bash
# =============================================================================
# run_main_scdrs.sh — CAND-C v2 主分析一键驱动
#
# 前置条件：
#   FUMA 网页平台三个任务已完成，作者下载各任务的 magma.genes.out 放入
#   analysis/fuma_output/ ，文件名需包含对应端点 ID（见下方 ENDPOINTS）。
#
# 流程：
#   Stage 1  三个 magma.genes.out → 合并 zscore 矩阵（首列 GENE）
#   Stage 2  合成 zscore 矩阵 → scDRS .gs 基因集
#   Stage 3  三个图谱 × 三个端点 → scDRS 细胞级打分 + 细胞类型富集（FDR）
#   Stage 4  汇总日志
#
# 用法：
#   bash scripts/run_main_scdrs.sh            # 全跑
#   bash scripts/run_main_scdrs.sh --dry-run  # 只检查输入是否齐备
# =============================================================================
set -uo pipefail

ROOT="C:/Users/docta/WorkBuddy/纯生信/project"
PY="C:/Users/docta/.workbuddy/binaries/python/envs/bioinfo/Scripts/python.exe"
SCDRS_CLI="C:/Users/docta/.workbuddy/binaries/python/envs/bioinfo/Scripts/scdrs"

FUMA="$ROOT/analysis/fuma_output"
ANA="$ROOT/analysis"
INT="$ANA/scrna_integrated"
OUTROOT="$ANA/scdrs_main"
GS="$ANA/hl_endpoints.gs"
ZMAT="$ANA/gene_z.tsv"

# ---- 端点配置：短名=报告名|FUMA 文件名关键字|图谱分组列 ----
ENDPOINTS=(
  "SEN|H8_HL_SEN_NAS|cell_type_analysis,dataset,condition"
  "CON|H8_HL_CON_NAS|cell_type_analysis,dataset,condition"
  "NOISE|H8_NOISEINNER|cell_type_analysis,dataset,condition"
)

# ---- 图谱配置：路径|显示名 ----
ATLASES=(
  "$INT/sv_atlas_labeled.h5ad|sv_atlas"
  "$INT/whole_cochlea_atlas_labeled.h5ad|whole_cochlea_atlas"
  "$INT/sgn_atlas.h5ad|sgn_atlas"
)

N_CTRL=1000
N_MAX=1000
WEIGHT_OPT="vs"

# BLAS 线程限制：32 核机器上 OpenBLAS 默认开满线程会在 Windows 下随机段错误
# （EXIT=139，无 traceback）。run_scdrs_pipeline.py 内部已设默认 1，
# 这里再兜一层，防止通过其他入口调用时漏掉。
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-1}"
export TQDM_DISABLE=1
# numba/llvmlite 的 LLVM 优化 pass 在 sv_atlas（50,757 细胞）规模上触发
# access violation（OSError: exception: access violation writing 0x8，2026-10-01 实测）。
# NUMBA_OPT=0 跳过优化 pass，绕开该 bug；sgn_atlas（8,916 细胞）小规模路径本就不触发。
export NUMBA_OPT="${NUMBA_OPT:-0}"

# =============================================================================
log() { printf '[%s] %s\n' "$(date '+%H:%M:%S')" "$*"; }

# ---------- Stage 0: 输入体检 ----------
log "Stage 0 — 检查 FUMA 输出与图谱"
declare -a EP_NAMES=() EP_KEYS=() EP_FILES=() EP_GROUPS=()
ok=1

for spec in "${ENDPOINTS[@]}"; do
  IFS='|' read -r name key groups <<<"$spec"
  # 在 fuma_output 下按关键字找 magma.genes.out（兼容子目录与不同扩展名）
  found=""
  while IFS= read -r f; do
    case "$(basename "$f")" in
      *genes.out|*.genes.out.txt) found="$f"; break ;;
    esac
  done < <(find "$FUMA" -type f -iname "*${key}*" 2>/dev/null | sort)
  if [ -z "$found" ]; then
    log "  ✗ 缺失端点 $name（应为 *${key}*genes.out）"
    ok=0
  else
    log "  ✓ $name → $(basename "$found")"
    EP_NAMES+=("$name"); EP_KEYS+=("$key"); EP_FILES+=("$found"); EP_GROUPS+=("$groups")
  fi
done

declare -a AT_FILES=() AT_NAMES=()
for spec in "${ATLASES[@]}"; do
  IFS='|' read -r path disp <<<"$spec"
  if [ -f "$path" ]; then
    log "  ✓ 图谱 $disp ($(du -h "$path" | cut -f1))"
    AT_FILES+=("$path"); AT_NAMES+=("$disp")
  else
    log "  ✗ 图谱缺失: $path"; ok=0
  fi
done

if [ "$ok" -ne 1 ]; then
  log "输入不齐备，退出。请先把 FUMA 的 magma.genes.out 放进："
  log "  $FUMA"
  exit 2
fi
[ "${1:-}" = "--dry-run" ] && { log "dry-run：输入齐备，可开跑。"; exit 0; }

# ---------- Stage 1: 合成 zscore 矩阵 ----------
log "Stage 1 — 合并 FUMA 基因级 ZSTAT 为 zscore 矩阵"
args=()
for i in "${!EP_NAMES[@]}"; do
  args+=("${EP_NAMES[$i]}=${EP_FILES[$i]}")
done
"$PY" "$ROOT/scripts/prep_gene_z_from_fuma.py" "$ZMAT" "${args[@]}" || exit 3

# ---------- Stage 2: 生成 scDRS 基因集 ----------
log "Stage 2 — munge-gs 生成基因集 $GS"
rm -f "$GS"
"$SCDRS_CLI" munge-gs "$GS" \
  --zscore_file "$ZMAT" --weight zscore --n_max "$N_MAX" || exit 4
ls -la "$GS"

# ---------- Stage 3: 三图谱 × 三端点 scDRS ----------
TRAIT_LIST=$(IFS=,; echo "${EP_NAMES[*]}")
log "Stage 3 — scDRS 打分与细胞类型富集（traits=$TRAIT_LIST）"
mkdir -p "$OUTROOT"

fail=0
for i in "${!AT_FILES[@]}"; do
  atlas="${AT_FILES[$i]}"; aname="${AT_NAMES[$i]}"
  odir="$OUTROOT/$aname"
  mkdir -p "$odir"
  log "--- $aname ---"
  # 每个图谱一次 preprocess、逐 trait 打分（脚本内部已按此组织）
  "$PY" "$ROOT/scripts/run_scdrs_pipeline.py" \
    --atlas "$atlas" --gs "$GS" --out "$odir" \
    --traits "$TRAIT_LIST" --n-ctrl "$N_CTRL" \
    --weight-opt "$WEIGHT_OPT" \
    --group-cols "cell_type_analysis,dataset,condition" \
    > "$odir/run.log" 2>&1
  rc=$?
  if [ $rc -ne 0 ]; then
    log "  ✗ $aname 失败（exit=$rc），详见 $odir/run.log"
    fail=1
  else
    log "  ✓ $aname 完成 → $odir"
    # 跨端点对比（核心增量）：同图谱内多端点 q95×对照校准的双侧 MC 检验
    "$PY" "$ROOT/scripts/compare_endpoints.py" \
      --atlas "$atlas" --dir "$odir" --group-col cell_type_analysis \
      >> "$odir/run.log" 2>&1 \
      && log "  ✓ $aname 跨端点对比完成 → $odir/compare_summary.md" \
      || { log "  ✗ $aname 跨端点对比失败，详见 $odir/run.log"; fail=1; }
  fi
done

# ---------- Stage 4: 汇总 ----------
log "Stage 4 — 汇总"
{
  echo "# scDRS 主分析汇总"
  echo
  echo "- 时间: $(date '+%Y-%m-%d %H:%M:%S')"
  echo "- 端点: $TRAIT_LIST"
  echo "- 图谱: ${AT_NAMES[*]}"
  echo "- 参数: n_ctrl=$N_CTRL, n_max=$N_MAX, weight_opt=$WEIGHT_OPT"
  echo
  for i in "${!AT_NAMES[@]}"; do
    aname="${AT_NAMES[$i]}"
    echo "## $aname"
    ls -1 "$OUTROOT/$aname" 2>/dev/null | sed 's/^/  - /'
    echo
  done
} > "$OUTROOT/SUMMARY.md"
cat "$OUTROOT/SUMMARY.md"

if [ "$fail" -ne 0 ]; then
  log "有图谱运行失败，请检查各自 run.log"
  exit 5
fi
log "全部完成。下一步：跨端点对比（compare_endpoints.py）。"
