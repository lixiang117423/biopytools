# coverage - BAM/SAM 覆盖度分析 | BAM/SAM Depth Analysis

一句话理解：**算出一个或多个 BAM/SAM 在每个碱基位置上的测序深度，输出"样本 × 染色体 × 位置 × 深度"的长表，并附带每个样本/染色体的深度统计，可选滑窗深度。**

## 功能概述 | Overview

- 基于 `samtools depth -a` 计算深度，`-a` 会输出所有位置（**包括 0 覆盖度**，方便画连续的深度曲线）
- 支持多文件、多目录输入（`-i` 可重复），自动识别 `.bam` / `.sam`
- 支持限定染色体（`-c`）或区间（`-r start:end`，1-based），需要时自动补建索引
- 自动生成统计报告（每样本、每染色体的平均/中位/最大/最小深度）
- 可选滑窗分析（`--enable-windows`），按窗口统计平均深度
- 无断点续传（每次重算）

## 快速开始 | Quick Start

~~~bash
biopytools coverage -i sample.bam -o depth_results.txt
~~~

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| 深度(depth) | 某个碱基位置上有多少条读段覆盖，深=测得多 |
| 0 覆盖度 | 该位置没有任何读段覆盖，本工具会明确输出 0（而不是省略） |
| 碱基质量 / 比对质量 | 单个碱基的可靠度 / 整条读段比对位置的可靠度 |
| 滑窗 | 把染色体切成固定大小的"窗格"逐段统计，看大尺度趋势更平滑 |
| 步长(step) | 相邻窗口起点相距多远；等于窗口大小则无重叠，更小则重叠 |

## 输入 | Input

- 一个或多个 BAM/SAM 文件、目录（目录自动收集所有 `*.bam` 和 `*.sam`），`-i` 可重复或一次传多个
- 可选 `-c 染色体`（默认 `all`）、`-r 区间`（格式 `start:end`，1-based，默认 `all`）

~~~bash
# 指定染色体区间
biopytools coverage -i sample.bam -o results.txt -c chr12 -r 136491092:138554123

# 多文件并行
biopytools coverage -i sample1.bam sample2.bam -o results.txt -t 16
~~~

## 参数说明 | Parameters

### 输入输出 | Input & output

**通俗理解|In plain words:** `-i` 是输入（可多个文件/目录），`-o` 是输出（文件路径，或目录则自动生成文件名）。`-c` 限定染色体、`-r` 限定区间，都默认 `all`（全基因组）。

### 质量过滤 | Quality filters

**通俗理解|In plain words:** `-q`（默认 0）是最低碱基质量，`-Q`（默认 0）是最低比对质量。调大=只要更可靠读段，深度变低但更可信；默认不过滤。

### 滑窗分析 | Sliding window

**通俗理解|In plain words:** `--enable-windows` 开启后，`--window-size`（默认 1000bp）定窗口大小，`--window-step`（默认 0）定步长。**步长 0 表示无重叠**（相邻窗口紧挨着）；步长小于窗口大小则窗口有重叠、结果更平滑但文件更大。

### 其他 | Others

**通俗理解|In plain words:** `--samtools-path` 指定 samtools 位置（默认 PATH 里的 `samtools`）；`--output-format` 目前只支持 `txt`；`--compress` 参数已声明但当前版本未实际压缩输出（见 FAQ）。一般都不用动。

## 分析流程 | Pipeline

~~~text
输入 BAM/SAM(文件/目录, 可多个)
    │
    ▼
1. 检查 samtools 可用
2. 收集/去重/排序输入文件; 需区间时补建索引
3. 逐文件: samtools depth -a 提取深度 → temp_{sample}.depth
4. 过滤(染色体/区间/质量)后追加到输出长表
5. (可选)滑窗统计
6. 生成统计报告 .stats.txt
~~~

## 输出 | Output

~~~text
depth_results.txt                     # 主结果长表: Sample Chromosome Position Depth
depth_results.stats.txt               # 统计报告(每样本/每染色体深度统计)
depth_results_windows_1000bp.txt      # 滑窗结果(开启 --enable-windows 时)
depth_analysis.log                    # 运行日志
~~~

- 主结果列：`Sample  Chromosome  Position  Depth`（长表，多样本堆叠，每行一个位置）
- 滑窗结果列：`Sample  Chromosome  Window_Start  Window_End  Window_Center  Average_Depth  Data_Points`
- 统计报告：总体统计 + 各样品统计 + 各染色体统计（平均/中位/最小/最大/标准差）
- `-o` 传目录时，自动按染色体/区间生成文件名（如 `depth_chr1.txt`、`depth_results.txt`）

## 结果解读 | Interpreting Results

**通俗理解|In plain words:** 主结果表一行一个位置，`Depth` 列就是该位置的深度；统计报告 `.stats.txt` 直接给出每个样本/染色体的平均、中位深度。

- `Depth = 0` 的位置大量出现：低测序区、重复区或比对不上区（因为用了 `-a` 才被显式列出）
- 平均深度与中位深度差距大：深度分布不均
- 滑窗结果：把 `Average_Depth` 按窗口画曲线，可直观看出整条染色体的覆盖起伏、缺口
- 多样本对比：主结果表按 `Sample` 列拆分，或直接比较统计报告里的平均深度

## 参数选择建议 | Parameter Guidance

- **全基因组深度**：默认 `biopytools coverage -i sample.bam`
- **只看某区域**：`-c chr12 -r 136491092:138554123`
- **要连续深度曲线（含 0）**：默认即可（内置 `-a`）
- **看大尺度趋势**：`--enable-windows --window-size 100000 --window-step 10000`
- **多样本一起比**：一次传多个 `-i`

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input, -i` | 必填 |  | BAM/SAM文件｜Input BAM/SAM file paths or directories |
| `--output, -o` | `./depth_results.txt` | Path | 输出文件｜Output file path or directory (default: ./depth_results.txt) |
| `--chromosome, -c` | `all` |  | 目标染色体｜Target chromosome name (default: all) |
| `--region, -r` | `all` |  | 染色体区域:start:end｜Chromosome region, format: start:end (default: all) |
| `--threads, -t` | `12` | int | 线程数｜Number of threads (default: 88) |
| `--quality, -q` | `0` | int | 最小碱基质量｜Minimum base quality threshold (default: 0) |
| `--mapping-quality, -Q` | `0` | int | 最小比对质量｜Minimum mapping quality threshold (default: 0) |
| `--samtools-path` | `samtools` |  | samtools路径｜samtools software path (default: samtools) |
| `--output-format` | `txt` | txt | 输出格式｜Output format (default: txt) |
| `--compress` | — |  | 压缩输出文件｜Compress output file |
| `--enable-windows` | — |  | 启用滑动窗口分析｜Enable sliding window analysis |
| `--window-size` | `1000` | int | 窗口大小(bp)｜Window size (bp) (default: 1000) |
| `--window-step` | `0` | int | 窗口步长(bp)｜Window step size (bp) (default: 0) |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 输入BAM/SAM文件路径或文件夹 (支持多个文件/文件夹，文件夹时自动识别所有BAM/SAM文件)｜Input BAM/SAM file paths or directories (supports multiple files/directories, auto-detect all BAM/SAM files in directories) |
| `-o, --output` | `./depth_results.txt` |  | 输出文件路径或目录 (目录时自动生成文件名)｜Output file path or directory (auto-generate filename when directory) |
| `-c, --chromosome` | `all` |  | 目标染色体名称 (默认all表示所有染色体)｜Target chromosome name (default all means all chromosomes) |
| `-r, --region` | `all` |  | 染色体区间，格式: start:end (如 100:1235, 1-based坐标)｜Chromosome region, format: start:end (e.g., 100:1235, 1-based coordinates) |
| `-t, --threads` | `88` | int | 线程数｜Number of threads |
| `-q, --quality` | `0` | int | 最小碱基质量阈值｜Minimum base quality threshold |
| `-Q, --mapping-quality` | `0` | int | 最小比对质量阈值｜Minimum mapping quality threshold |
| `--samtools-path` | `samtools` |  | samtools软件路径｜samtools software path |
| `--output-format` | `txt` | txt | 输出格式｜Output format |
| `--compress` | — | store_true | 压缩输出文件｜Compress output file |
| `--enable-windows` | — | store_true | 启用滑窗分析｜Enable sliding window analysis |
| `--window-size` | `1000` | int | 窗口大小(bp)｜Window size (bp) |
| `--window-step` | `0` | int | 窗口步长(bp)，0表示无重叠，小于窗口大小则重叠｜Window step size (bp), 0 means no overlap, less than window size means overlap |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- `samtools`（通过 `--samtools-path` 指定，默认 PATH 中的 `samtools`）

## 常见问题 | FAQ

**Q1：会不会断点续传？**
不会。每次运行都重新提取深度（BAM 索引会复用/补建，但深度计算不缓存）。

**Q2：为什么输出里有很多 Depth=0 的行？**
本工具用 `samtools depth -a` 输出所有位置（含 0 覆盖度），是为了画出连续的深度曲线。若不需要 0 覆盖位置，可在结果里按 `Depth > 0` 过滤。

**Q3：`--compress` 参数有用吗？**
当前版本该参数已声明并存入配置，但实际处理逻辑未实现压缩（输出仍为普通文本）。建议暂不使用，如需压缩可事后手动 `gzip`。

**Q4：区间 `-r` 的坐标从几开始数？**
1-based（第一个碱基是 1），格式 `start:end`，且 start 必须小于 end。

**Q5：支持哪些输入？**
BAM 和 SAM 都支持；可以是单个文件、多个文件、目录（自动识别 `*.bam`/`*.sam`），也能混传。目录里找不到文件会直接报错。
