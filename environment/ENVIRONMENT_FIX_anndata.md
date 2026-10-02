# ENVIRONMENT_FIX_anndata.md — anndata 崩溃诊断报告

**诊断人**: env-doctor
**诊断日期**: 2026-09-30
**环境**: `C:/Users/docta/.workbuddy/binaries/python/envs/bioinfo/`
**解释器**: `Scripts/python.exe` (Python 3.13.14)
**结论速览**: **未发现环境故障。报告中的 `Illegal instruction` (SIGILL/exit 132) 无法复现；anndata 及其完整 scDRS 分析路径实测全部正常。**

---

## 1. 崩溃现象（任务描述）

据任务描述：执行

```
python -c "import anndata; print(anndata.__version__)"
```

返回 `Illegal instruction`（退出码 132 / SIGILL），而非 `ModuleNotFoundError`，
会导致 `run_main_scdrs.sh` 的 Stage 3 无法读取 `.h5ad` 图谱。

---

## 2. 诊断过程与证据

### 2.1 逐层 import 定位（未复现崩溃）

| 测试 | 命令 | 结果 |
|---|---|---|
| 解释器 | `python -V` | `Python 3.13.14`，exit 0 |
| numpy | `python -c "import numpy"` | `numpy 2.5.3`，exit 0 |
| h5py | `python -c "import h5py"` | `h5py 3.16.0`，exit 0 |
| **anndata（复测点）** | `python -c "import anndata"` | **`anndata 0.13.4`，exit 0** |
| anndata 重复 5 次 | 连续 5 次 import | 5/5 全部 exit 0（无随机崩溃） |

**结论**: 报告中的 SIGILL **在 `import anndata` 一层即无法复现**。

### 2.2 无线程限制下复测（复现原崩溃条件）

原脚本注释指出 32 核机器上 OpenBLAS 默认开满线程会在 Windows 下随机段错误。
在**不加任何线程限制**（不设 `OPENBLAS_NUM_THREADS` 等）条件下连跑 3 次：

```
run 1 exit=0   ok
run 2 exit=0   ok
run 3 exit=0   ok
```

**结论**: 即使在原崩溃最可能触发的"多线程"条件下，依然稳定。

### 2.3 阻塞点检查：真实读取 `.h5ad`（真正验收标准）

用 anndata 直接打开最小的图谱 `project/analysis/scrna_integrated/sgn_atlas.h5ad`（599MB / 628,011,278 字节）：

```
n_obs = 8916
n_vars = 17486
obs cols = ['condition', 'celltype', 'uMAP_1', 'uMAP_2', 'replicate',
            'cluster_label', 'celltype_colors', 'louvain', 'dataset', 'axis',
            'cell_type', 'cell_type_original', 'sample', 'n_genes',
            'n_genes_by_counts', 'total_counts', 'total_counts_mt', 'pct_counts_mt']
var cols = ['gene_symbol', 'ensembl', 'n_cells', 'mt', 'n_cells_by_counts',
            'mean_counts', 'pct_dropout_by_counts', 'total_counts']
```

**只读取，未修改原文件。** 读取成功，无任何异常。

进一步加载完整 `.X`（dense ndarray，8916×17486，sum=34,239,648）并执行
`X.mean(axis=0)` 等真实数值路径 → 全部正常，同时 `import numba 0.67.0` 正常。

### 2.4 scdrs 命令行可用性

```
scdrs --help           → exit 0，正常输出命令组
scdrs munge-gs --help  → exit 0，正常输出版本与说明
```

### 2.5 端到端压力测试（最强验证）

直接用 `run_scdrs_pipeline.py` 对 **sgn_atlas.h5ad** 实跑一次完整 scDRS 四步流程
（`--n-ctrl 50`，临时 mock `.gs`，输出到临时目录，验证后已删除）：

```
atlas: 8916 cells x 17486 genes
分析标签分布: {'Type 1C': 3069, 'Type 1A': 2427, 'Type 1B': 2086, 'Schwann Cells': 1041, 'Type 2': 293}
preprocess 完成（adj_prop=cell_type_analysis）
[MOCK_CON]   打分完成: (8916, 56)；group analysis × 3 列 完成
[MOCK_NOISE] 打分完成: (8916, 56)；group analysis × 3 列 完成
[MOCK_SNHL]  打分完成: (8916, 56)；group analysis × 3 列 完成
DONE   （PIPELINE_EXIT=0）
```

成功生成 `score_*.csv.gz`（每个约 2.4MB）、`group_*.csv`、`pipeline_log.md`。

**结论**: 主分析 Stage 3 所依赖的**完整代码路径**（读图谱 → 归一化 → preprocess →
score_cell → downstream_group_analysis）**全部跑通**。

---

## 3. 根因判定

**未发现根因，因为故障不可复现。** 逐项排除任务列出的三类假设：

| 假设 | 判定 | 证据 |
|---|---|---|
| numpy 2.x 与 anndata/h5py ABI 不兼容 | **排除** | numpy 2.5.3 + anndata 0.13.4 + h5py 3.16.0 三者按官方 numpy 2.x 兼容矩阵搭配，import 与真实读取均正常 |
| wheel 与 CPU 指令集不匹配（AVX2/AVX512） | **排除** | 崩溃若由此引起必然**稳定复现**；实测 8+ 次（含重复 import）全部成功。numpy 显示 `DYNAMIC_ARCH`（运行时自适应指令集），不绑定具体 AVX 级别 |
| 安装损坏 / 混装 py3.13 不匹配 wheel | **排除** | 无残留 numpy：`Lib/site-packages` 下仅 `numpy` / `numpy-2.5.3.dist-info` / `numpy.libs` 各一份，无重复或半装残留 |

**版本一致性核查**（`pip list`，均相互匹配）：

```
anndata 0.13.4   h5py 3.16.0      numcodecs 0.17.0   numpy 2.5.3
scipy 1.18.1     scanpy 1.12.4    scdrs 1.0.2        pandas 2.3.3
numba 0.67.0     llvmlite 0.49.0  zarr 3.4.0
```

### 最可能的解释（供参考，非结论）

崩溃报告与本次实测的**唯一显著差异**是执行环境/时机。最可能是：
1. **瞬态故障**：包在安装/下载**中途**被调用（文件未落盘完整），事后安装完成即自愈；
2. **进程环境差异**：报错时可能经由某个未设置线程限制的入口触发 OpenBLAS 多线程
   段错误（脚本注释已记录此类历史问题），事后由 `run_scdrs_pipeline.py` 内置的
   `OMP/OPENBLAS/MKL/NUMEXPR_NUM_THREADS=1`（见该文件 L29–L47）兜住；
3. **报告混淆**：SIGILL(132) 与脚本注释中记录的段错误 SIGSEGV(139) 是两种不同信号，
   需确认原始报错抓取无误。

---

## 4. 修复动作

**本轮未执行任何安装/降级/卸载操作**——因为不存在需要修复的故障，任何变更都只会
引入不必要的风险（尤其"降 numpy 到 1.26.x"会同时拖累 scdrs 1.0.2 / scipy 1.18.1，
属于明确有害的方向，已按任务约束避免）。

任务要求的验证标准**全部通过**：

- [x] `python -c "import anndata; print(anndata.__version__)"` 正常返回 `0.13.4`
- [x] 真实读取 `sgn_atlas.h5ad`：`n_obs=8916`，`n_vars=17486`，obs 列名见 §2.3
- [x] `scdrs --help` / `scdrs munge-gs --help` 正常输出
- [x] `import scdrs, scanpy, pandas, numpy, scipy` 全部正常
- [x] （额外加强）完整 scDRS 打分管线端到端跑通

**本轮唯一文件变更**：写入本报告 `ENVIRONMENT_FIX_anndata.md`；
临时验证产物 `project/analysis/scdrs_main/_envcheck/` 已删除；
`sgn_atlas.h5ad` 及 `project/` 下其它数据文件**均未修改**。

---

## 5. 若将来重建环境：精确版本清单

建议直接 `pip freeze` 固化，以下为本次实测可用的关键组合（Windows + Python 3.13.14 + AMD64）：

```
python==3.13.14
numpy==2.5.3
scipy==1.18.1
pandas==2.3.3
anndata==0.13.4
h5py==3.16.0
scanpy==1.12.4
scdrs==1.0.2
numcodecs==0.17.0
zarr==3.4.0
numba==0.67.0
llvmlite==0.49.0
```

重建命令（PowerShell / bash 均可）：

```bash
PY="C:/Users/docta/.workbuddy/binaries/python/envs/bioinfo/Scripts/python.exe"
"$PY" -m pip install \
  "numpy==2.5.3" "scipy==1.18.1" "pandas==2.3.3" \
  "anndata==0.13.4" "h5py==3.16.0" "scanpy==1.12.4" "scdrs==1.0.2"
```

全量固化（推荐）：

```bash
"$PY" -m pip freeze > C:/Users/docta/WorkBuddy/纯生信/project/analysis/requirements.lock.txt
"$PY" -m pip install -r <(cat requirements.lock.txt)
```

---

## 6. 残留风险

| 风险 | 等级 | 说明 / 建议 |
|---|---|---|
| **OpenBLAS 多线程随机段错误** | **中** | 脚本历史注释记录过 EXIT=139。`run_scdrs_pipeline.py` 与 `run_main_scdrs.sh` 均已设线程=1 兜底。**务必保持这两处线程限制**，不要为了提速而放开。 |
| SIGILL 可能以"间歇性"形式回归 | 低 | 本次 8+ 次复测全过。建议在实际主分析运行时**保留完整 stderr 日志**，一旦再次出现 exit 132，立即抓取 `%TEMP%` 下的 crash dump 与 `pip freeze` 快照。 |
| pandas 3.0 行为变更警告 | 低 | scdrs 内部有 `FutureWarning`（chained assignment / `observed=False`）。**当前仅告警不影响结果**，但将来若升级 pandas 3.x 需回归验证 scDRS 输出。 |
| sgn_atlas.h5ad 的 `.X` 为 dense ndarray | 低 | 8916×17486 dense ≈ 1.2GB 内存常驻。同时跑 3 个图谱时注意内存峰值。 |

**总体判断：环境当前健康，可放心运行主分析。** 建议在正式跑 `run_main_scdrs.sh`
前先执行一次 `bash project/scripts/run_main_scdrs.sh --dry-run`（当前会报 CON / NOISE
两个端点缺失，属输入未备齐，非环境问题）。
