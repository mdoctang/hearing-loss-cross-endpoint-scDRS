# 计算环境与版本锁定（CAND-C v2）

> 目的：投稿时 Methods 需要写明软件与版本；审稿人复现需要精确的环境快照。
> 所有分析均在本机隔离 venv 中运行，不污染系统环境。

## 1. 解释器与虚拟环境

| 项 | 值 |
|---|---|
| Python | 3.13.14 |
| venv 路径 | `C:\Users\docta\.workbuddy\binaries\python\envs\bioinfo` |
| 平台 | Windows (win32)，Git Bash |
| CPU | 32 逻辑核 |

## 2. 依赖版本快照（pip freeze 摘录，2026-09-30）

```
anndata==0.13.4
h5py==3.16.0
matplotlib==3.11.2
numba==0.67.0
numpy==2.5.3
pandas==2.3.3        # 见 §3.1，刻意钉在 2.x
pytz==2026.4
scanpy==1.12.4
scdrs==1.0.2
scikit-learn==1.9.1
scipy==1.18.1
seaborn==0.13.2
statsmodels==0.15.0
tqdm==4.70.1
```

## 3. 两处必须记录在案的兼容性处置

投稿时应在 Methods / Supplementary 中说明，否则审稿人用"最新版依赖"复现会失败。

### 3.1 pandas 必须锁在 2.x（不能升到 3.0）

`scdrs.pp.preprocess` 内部有一处链式赋值 + `inplace=True`：

```python
# scdrs/pp.py:252
temp_df["cell"].clip(lower=int(0.01 * temp_df["cell"].max()), inplace=True)
```

- pandas 3.0 强制 Copy-on-Write，该写法**静默失效**（只告警不报错），
  于是 `adj_prop` 下"小细胞类型计数下限裁剪"这一步不再生效，
  → `cell_weight` 被悄悄改变 → 下游所有富集统计量随之偏移。
- 因此本环境将 pandas 钉在 **2.3.3**（`scanpy 1.12.4` / `anndata 0.13.4` 均要求 `pandas>=2.3`，故 2.3.x 是唯一同时满足两者的区间）。
- 安装命令：`pip install "pandas>=2.3,<3"`

### 3.2 NumPy 2.x 需要 `np.float_` 垫片

`scdrs.method.gearys_c` 仍使用 NumPy 2.0 已移除的别名（两处）：

```python
# scdrs/method.py:1055
graph_data = graph.data.astype(np.float_, copy=False)
# scdrs/method.py:1062
vals = vals.astype(np.float_)
```

- 触发时机：`downstream_group_analysis` → `test_gearysc` → `gearys_c`，
  即**细胞级打分已完成之后**才崩，报 `AttributeError: np.float_ was removed in the NumPy 2.0 release`。
- 处置：在 `scripts/run_scdrs_pipeline.py` 顶部加垫片（**不修改 site-packages**，保证仓库自包含可复现）：

```python
if not hasattr(np, "float_"):
    np.float_ = np.float64
```

### 3.3 必须限制 BLAS 线程数，否则随机段错误（**最隐蔽的一个**）

**症状**：进程无任何 Python traceback 直接消失，shell 报
`Segmentation fault`，退出码 **139**。日志停在 `score_cell` 执行中途，
若碰巧在写 `score_*.csv.gz` 时崩，会留下一个 **gzip 流未闭合的截断文件**
（读取时报 `EOFError: Compressed file ended before the end-of-stream marker`），
极易被误判为"写盘失败"或"内存不足"。

**根因**：本机 32 逻辑核，OpenBLAS 默认按满核开线程。OpenBLAS 在 Windows 上
线程数过多时会出现栈耗尽类崩溃，且是概率性的（同一命令有时能过）。
**实测：与内存无关**——三个图谱的最坏内存占用仅 0.6 / 1.0 / 1.65 GB。

**处置**：在所有 BLAS 线程变量上锁 1（必须在 `import numpy/scipy` **之前**设置）：

```python
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")
```

已内置在 `scripts/run_scdrs_pipeline.py` 顶部，`run_main_scdrs.sh` 中另设一层
`export` 兜底。限制为 1 线程后连续多次运行均稳定通过。
如确需提速，可显式覆盖（`OPENBLAS_NUM_THREADS=8 ...`），但**不要**设为满核。

### 3.4 预处理缓存的写入需显式授权 nullable string

`anndata >= 0.11` 默认拒写 nullable string 列，抛
`RuntimeError: anndata.settings.allow_write_nullable_strings is None ...`。
已在脚本内设 `ad.settings.allow_write_nullable_strings = False`
（= 写 anndata<0.11 兼容格式，最大化下游工具互通）。
缓存采用 `.tmp` + `os.replace` **原子写入**，避免中途失败留下半截文件被 `--reuse` 误读。

### 3.5 两个"静默出错"级的工程约定（不涉及科学口径，但会毁掉一次运行）

**(a) 必须使用 venv 解释器，不要用 `binaries/python/versions/` 下的基础解释器。**
`C:\Users\docta\.workbuddy\binaries\python\versions\3.13.12\python.exe` 在本机会对任意脚本
概率性抛 `Illegal instruction`（SIGILL，无 traceback，退出码 132）。
**统一使用 `C:\Users\docta\.workbuddy\binaries\python\envs\bioinfo\Scripts\python.exe`。**

**(b) 大文件流式转换避开 Python，改用 awk + gzip 管道。**
对流式读写 ~2,000 万行的 gzip（典型场景：sumstats 列重排），Python 侧同样观测到概率性
SIGILL，且重试次数越多崩得越早（首次 4 分钟、后续 20 秒），与内存无关。
已验证稳定的替代通道：`gzip -dc … | awk … | gzip -c`（内存恒定，两千万行约 1.5 分钟/端点）。
参见 `scripts/make_fuma_input_grch38.sh` 与 `scripts/make_fuma_input_grch38.py` 的双实现。

**(c) 写 TSV 必须显式指定换行符。**
Python 在 Windows 上 `open(..., "wt")` / `gzip.open(..., "wt")` 的 `newline=None` 会把 `\n`
转成 `\r\n`，于是**表头最后一列的字节变成 `p_value\r`**。对表头做精确匹配的下游工具
（如 FUMA）会因此报格式错误。本次已实际踩到：`\r` 是 FUMA `ERROR:001` 的第二个成因。
**约定：所有对外交付的 TSV/CSV 一律以 `newline="\n"` 写出，并在交付前用 `gzip.open(p,"rb")`
读回首行核对字节。**


## 4. 其他运行时约定

- `TQDM_DISABLE=1`：关闭 scDRS 内部 tqdm。1000 个对照基因集会刷出上万行进度条，
  会把 `run.log` 撑爆（实测单次空跑日志 > 200 KB），也干扰 traceback 定位。
- 大图谱的内存安全策略：`sc.pp.scale` / PCA 一律**只在 HVG 子集上做**，
  绝不在全基因空间稠密化（`whole_cochlea_atlas` 101,082 细胞 × 26,525 基因，
  全矩阵稠密 float32 ≈ 10.7 GB，会直接爆内存）。

## 5. 运行入口

```bash
# 空跑验证（模拟基因集，不依赖 FUMA）
python scripts/run_scdrs_pipeline.py \
  --atlas analysis/scrna_integrated/sgn_atlas.h5ad \
  --gs analysis/test_mock.gs \
  --out analysis/scdrs_test --traits MOCK_SNHL

# 主分析一键驱动（需 FUMA 的 magma.genes.out 已就位）
bash scripts/run_main_scdrs.sh --dry-run   # 先体检输入
bash scripts/run_main_scdrs.sh             # 正式跑
```
