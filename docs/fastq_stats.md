# FASTQ 统计 | FASTQ Statistics

一句话理解：**批量统计一批 FASTQ 文件的基本指标（reads 数、长度、碱基总量、文件大小等），自动配对双末端文件，汇总成一张 CSV 或 Excel 表**，方便一眼看清所有样本的数据量。

## 功能概述 | Overview

- 基于 seqkit 的高性能统计，支持批量多样本
- 自动识别并配对双末端（R1 / R2）文件
- 统计 reads 数、最短/最长/平均长度、碱基总量、文件大小（压缩前后）
- 输出 CSV 或 Excel（xlsx）两种格式，结果按样本名排序

## 快速开始 | Quick Start

```bash
biopytools fq-stats -i /data/fastq/ -o results.csv -p "*_1.clean.fq.gz"
```

最小输入：一个 FASTQ 目录（或单个文件）+ 一个输出文件名（.csv 或 .xlsx）。

## 零基础概念速览 | Concepts in plain words

| 术语<br>Term | 通俗理解<br>Plain meaning |
|------|----------|
| reads 数 | 一个 FASTQ 里有多少条读段，反映「测了多少量」 |
| 双末端 paired-end | 一条 DNA 片段两头各测一次，产出 R1、R2 两个文件，通常是配对的 |
| R1 / R2 | 双末端测序的两个方向文件，习惯叫 read1 / read2 |
| 碱基总量 total bases | 所有 read 的碱基加起来有多少，等于「reads 数 × 平均长度」 |
| seqkit | 序列处理工具，本工具的统计引擎 |

## 输入 | Input

### 输入文件或目录

- 目录：扫描其中符合模式（pattern）的 FASTQ 文件
- 单个文件：直接统计这一个文件
- 支持 `.fastq` / `.fq` / `.fastq.gz` / `.fq.gz`

### 匹配模式 pattern（可选）

用 `*` 代表样本名，如 `*_1.clean.fq.gz`。程序会把 `*` 匹配到的部分当样本名，并自动推导对应的 R2 文件（支持 `_1.`→`_2.`、`_R1.`→`_R2.`、`.1.`→`.2.`、`.R1.`→`.R2.` 等命名）。

## 参数说明 | Parameters

### 必需参数 | Required

**通俗理解|In plain words:** `-i` 指定「统计谁」（目录或单文件），`-o` 指定「结果存哪、存成什么格式」。输出扩展名只能是 `.csv` 或 `.xlsx`，其他格式会直接报错。

### 匹配模式 | Pattern

**通俗理解|In plain words:** `-p` 是「按什么规则认出样本和配对」。模式里必须带一个 `*`（代表样本名），程序靠它抽样本名、并自动找 R2。不写 `-p` 时，会把每个 FASTQ 文件各自当成一个样本（不配对）。目录里文件命名规范时，用 `-p "*_1.clean.fq.gz"` 这类模式最省事。

### 线程数 | Threads

**通俗理解|In plain words:** `-t` 控制 seqkit 统计的并行度，默认 12。文件多、机器核多时可以调大加速；一般不用动。

## 分析流程 | Pipeline

```text
输入目录/文件 + pattern
  |
  v
1. 检查 seqkit 是否可用
  |
  v
2. 按 pattern 查找样本、配对 R1/R2(或单端)
  |
  v
3. seqkit stats -T 批量统计
  |
  v
4. 汇总成 CSV 或 Excel -> 输出文件
```

## 输出 | Output

```text
results.csv / results.xlsx     # 汇总统计表
```

CSV 列（Excel 顺序略有调整）：

- `sample_name`：样本名（从文件名 pattern 提取）
- `R1_file` / `R2_file`：对应的文件名（单端时 R2 为空）
- `R1_reads` / `R2_reads`：各自 reads 数
- `R1_min_length` / `R1_max_length` / `R1_mean_length`：R1 最短 / 最长 / 平均长度（R2 同理）
- `R1_total_bases` / `R2_total_bases`：碱基总量
- `R1_file_size_formatted` / `R1_uncompressed_size_formatted`：文件大小（压缩 / 解压后）
- `total_reads` / `total_bases`：双末端合计（R1 + R2）

## 结果解读 | Interpreting Results

**通俗理解|In plain words:** 看这张表主要关注三件事：每个样本有多少 reads、平均长度是否合理、双末端 R1/R2 是否一致。

- **reads 数与碱基总量**：反映测序量。reads 数异常少（相比同批其他样本）可能该样本测序失败
- **R1 / R2 一致性**：正常双末端样本 R1、R2 的 reads 数应相等或接近；差很多说明有一端数据缺失或配对出错
- **平均长度**：短读通常 100-150 bp，长读（HiFi）可达几千 bp；明显偏离预期要留意数据是否有问题
- **好坏判据**：各样本 reads 数在同一量级、R1/R2 数量一致、长度符合平台特征，即为正常

## 参数选择建议 | Parameter Guidance

- **统计目录下全部文件**：不写 `-p`，每个文件当独立样本
- **目录里有大量无关文件**：用 `-p "*_1.clean.fq.gz"` 精确圈定
- **单文件**：直接 `-i sample_R1.fastq.gz`，当单样本统计
- **要 Excel 给同事看**：`-o stats.xlsx`（需 pandas + openpyxl）

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input, -i` | 必填 |  | 输入FASTQ文件或目录｜Input FASTQ file or directory |
| `--output, -o` | 必填 | Path | 输出文件路径(.csv或.xlsx)｜Output file path (.csv or .xlsx) |
| `--pattern, -p` | — |  | FASTQ文件匹配模式，如"*_1.clean.fq.gz"｜FASTQ file matching pattern, e.g., "*_1.clean.fq.gz" |
| `--threads, -t` | `12` | int | 线程数｜Number of threads |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 输入FASTQ文件或包含FASTQ文件的目录｜Input FASTQ file or directory containing FASTQ files |
| `-o, --output` | 必填 |  | 输出文件路径 (.csv 或 .xlsx)｜Output file path (.csv or .xlsx) |
| `-p, --pattern` | — |  | FASTQ文件匹配模式，如 "*_1.clean.fq.gz"，*代表样品名称｜FASTQ file matching pattern, e.g., "*_1.clean.fq.gz", * represents sample name |
| `-t, --threads` | `12` | int | 线程数｜Number of threads (默认｜default: 12) |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- seqkit（`seqkit stats`、`seqkit version`）
- pandas + openpyxl（仅 `-o *.xlsx` 时）

无固定 conda 环境名；seqkit 直接通过 PATH 查找（`conda install -c bioconda seqkit` 可安装）。

## 常见问题 | FAQ

**Q1：输出必须是 .csv 或 .xlsx 吗？**
是的，其他扩展名会在初始化时直接报错「输出文件格式不支持」。要么改后缀，要么换用别的统计工具。

**Q2：pattern 不带 * 会怎样？**
程序会校验 pattern 必须包含 `*` 通配符，否则报错。`*` 就是样本名的占位符。

**Q3：为什么我设了 pattern 却配不上 R2？**
程序只识别几种常见命名转换（`_1.`→`_2.`、`_R1.`→`_R2.`、`.1.`→`.2.`、`.R1.`→`.R2.`）。若你的 R1/R2 命名不在这几种里（比如 R1 叫 `read1`、R2 叫 `read2`），配不上时会被当单端处理。

**Q4：xlsx 报 pandas 未安装？**
xlsx 需要 pandas 和 openpyxl，缺了会报错退出（与 ena_downloader 不同，这里不退回 CSV）。装 `pip install pandas openpyxl` 或改用 `.csv` 输出。

**Q5：没有断点续传，重跑会怎样？**
本工具单次跑完，重跑直接覆盖输出文件，没有已存在就跳过的逻辑，放心重跑即可。
