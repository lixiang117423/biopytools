# MEME Parser - 序列模体发现与解析 | MEME Parser Motif Discovery & Parsing

一句话理解：**在一堆序列里自动找出反复出现的「特征短片段」（模体/motif），并整理成表格**。输入一个蛋白或 DNA FASTA，运行 MEME 发现 motif，再把结果解析成 TSV/CSV/Excel，并提取每个 motif 命中的具体序列片段。

## 功能概述 | Overview

- 两种模式：运行模式（给 FASTA，先跑 MEME 再解析）与解析模式（`--parse-only`，直接解析已有的 MEME xml/txt 输出）
- 支持蛋白质（默认）与 DNA 输入，完整透传 MEME 核心参数（motif 数量、宽度范围、分布模式、目标函数等）
- 解析 MEME XML 或 TXT 输出，输出 motif 位置表（序列、起点、终点、宽度、p-value、E-value 等）
- 自动提取每个 motif 位点的序列片段，追加到表格并另存 FASTA
- 输出 TSV / CSV / Excel 三种格式（可分别关闭）
- 断点续传：MEME 输出 meme.xml 已存在则跳过运行，直接解析

## 快速开始 | Quick Start

```bash
biopytools meme-parser -i proteins.fa -protein
```

最小输入：一个 FASTA 文件。默认在 `./meme_results_meme_out/` 生成 MEME 原始输出，并在当前目录输出 `meme_results.tsv/csv/xlsx`。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| motif（模体） | 一组序列里反复出现、位置保守的「短特征片段」，常是蛋白/DNA 的功能位点 |
| motif 宽度 | 这个特征片段多长（几个氨基酸/碱基） |
| zoops | 假设每条序列最多出现一次该 motif；anr=可出现多次；oor=有序列可能完全没有 |
| p-value / E-value | 该位点「纯属巧合」的可能性，越小越可信 |
| 目标函数(objfun) | MEME 打分时用什么标准衡量 motif「好坏」，classic 最常用 |
| Markov 阶数 | 背景序列模型，用来校正随机出现的假 motif；0 阶最基础 |

## 输入 | Input

### 运行模式

- FASTA 文件（蛋白质或 DNA），配合 `-protein` 或 `-dna` 声明类型

### 解析模式（`--parse-only`）

- MEME 输出文件 `meme.xml` 或 `meme.txt`；若要提取 motif 序列，另用 `--input-fasta` 提供原始序列

```text
>seq1
MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ...
>seq2
MKLAFLKSQRQVAFVKSHFSRQLEDRLGL...
```

## 参数说明 | Parameters

### 输入与模式 | Input & mode

**通俗理解|In plain words:** `-i` 可以是 FASTA（运行模式）也可以是 MEME 输出（加 `--parse-only` 解析模式）。序列类型用 `-protein`（默认）或 `-dna` 声明，**二者互斥**。

### 输出配置 | Output

**通俗理解|In plain words:** 决定结果文件名前缀和目录。默认前缀 `meme_results`、目录当前目录（`.`），**一般不用动**；批量处理想区分结果时改前缀。

### MEME 核心参数 | MEME core options

**通俗理解|In plain words:** 这是 MEME 本身的关键参数。`-nmotifs` 决定找几个 motif（默认 10）；`-minw`/`-maxw` 决定 motif 长度范围（默认 6–50）；`-mod` 决定 motif 出现次数假设；`-objfun` 决定打分方式。**新手用默认即可**，只有对 motif 长度/数量有明确预期时才调。

### 序列提取与格式 | Extraction & format

**通俗理解|In plain words:** 默认会提取 motif 序列片段、输出 TSV/CSV/Excel。解析模式要提取序列时，务必用 `--input-fasta` 指明原始 FASTA。

## 分析流程 | Pipeline

**通俗理解|In plain words:** 运行模式先跑 MEME（有 meme.xml 就跳过），再解析 motif 位置，提取片段序列，最后保存表格。

```text
FASTA(蛋白/DNA)
    │
    ▼
[运行模式] 运行 MEME → {prefix}_meme_out/meme.xml + meme.txt
    │  (断点续传: meme.xml 已存在则跳过)
    ▼
解析 XML/TXT → motif 位置表(sequence/motif/start/end/width/pvalue/evalue)
    │
    ▼
(可选) 提取 motif 序列片段 → 追加 motif_sequence 列 + 写 FASTA
    │
    ▼
保存 {prefix}.tsv / .csv / .xlsx
```

## 输出 | Output

```text
.
├── meme_results.tsv                       # motif 位置主表
├── meme_results.csv                       # 同内容 CSV
├── meme_results.xlsx                      # 同内容 Excel
├── meme_results_motifs.fasta              # 提取的 motif 序列片段
└── meme_results_meme_out/                 # MEME 原始输出目录(运行模式)
    ├── meme.xml                           # MEME XML 结果
    └── meme.txt                           # MEME 文本结果
```

- 主表列：`sequence_id / length / motif_id / start / end / width / strand / pvalue / evalue / num_sites`（开启提取时追加 `motif_sequence` 列）
- 每个 motif 位点一行；同一序列可命中多个 motif

## 结果解读 | Interpreting Results

### 1. motif 位置表（TSV/CSV/Excel）

**通俗理解|In plain words:** 一行=一个「motif 出现在某序列某个位置」的命中。看 `evalue` 判断可信度，看 `start/end` 定位。

- `evalue` 越小越可信；`pvalue` 同理。通常挑 evalue 最小的若干位点作为「真实 motif 位点」
- `motif_id` 区分不同的 motif（MEME 按显著性从强到弱编号）；`num_sites` 是该 motif 命中的总位点数
- `start/end` 是 motif 在序列上的 1-based 坐标；`width` 是 motif 长度

### 2. motif 序列 FASTA（`meme_results_motifs.fasta`）

**通俗理解|In plain words:** 把每个命中的片段截出来存成 FASTA，方便后续做 motif logo 或比对。fasta 头形如 `序列名_Motif编号`。

### 3. 好坏判据 | Judgment

- motif 越靠前（编号越小）越显著，evalue 也越小；前 3–5 个 motif 通常是主要信号
- 若大量 motif 的 evalue 都很大（接近 1），说明输入里没有明显保守片段，考虑换输入或调整参数

## 参数选择建议 | Parameter Guidance

- **`-nmotifs`**：不知道有几个 motif 就保持默认 10 先扫一遍，再按显著性挑靠前的；有明确预期时改成该数量
- **`-mod zoops`**（默认）：大多数「每条序列至多一个 motif」的场景；motif 可能重复出现用 `anr`
- **`--parse-only`**：已经跑过 MEME、只想换解析方式或提取序列时用，省去重跑 MEME
- **`--input-fasta`**：解析模式要提取 motif 序列时必须给，否则只能拿到位置表

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input-file` | 必填 |  | 输入FASTA文件或MEME输出文件(xml/txt)｜Input FASTA file or MEME output file (xml/txt) |
| `--parse-only` | — |  | 仅解析模式，不运行MEME｜Parse-only mode, do not run MEME |
| `-o, --output-prefix` | `meme_results` |  | 输出文件前缀｜Output file prefix |
| `--output-dir` | `.` |  | 输出目录｜Output directory |
| `--no-tsv` | — |  | 不输出TSV文件｜Do not output TSV file |
| `--no-csv` | — |  | 不输出CSV文件｜Do not output CSV file |
| `--no-excel` | — |  | 不输出Excel文件｜Do not output Excel file |
| `--meme-path` | `~/miniforge3/envs/protein/bin/meme` |  | MEME可执行文件路径｜MEME executable path |
| `-protein` | `True` |  | 输入序列为蛋白质(默认)｜Input sequences are protein (default) |
| `-dna` | — |  | 输入序列为DNA｜Input sequences are DNA |
| `-mod` | `zoops` | zoops/anr/oor | Motif分布模式｜Motif distribution mode |
| `-nmotifs` | `10` | int | Motif数量｜Number of motifs |
| `-minw` | `6` | int | 最小motif宽度｜Minimum motif width |
| `-maxw` | `50` | int | 最大motif宽度｜Maximum motif width |
| `-objfun` | `classic` | classic/de/ce/cd | 目标函数｜Objective function |
| `-markov-order` | `0` | int | Markov链阶数｜Markov chain order |
| `--no-extract-motifs` | — |  | 不提取motif序列｜Do not extract motif sequences (default: extract) |
| `--input-fasta` | — |  | 原始FASTA文件路径（解析模式时用于提取motif序列）｜Original FASTA file path (for motif extraction in parse-only mode) |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input-file` | 必填 |  | 输入FASTA文件或MEME输出文件(xml/txt)｜Input FASTA file or MEME output file (xml/txt) |
| `--parse-only` | — | store_true | 仅解析模式，不运行MEME｜Parse-only mode, do not run MEME |
| `-o, --output-prefix` | `meme_results` |  | 输出文件前缀｜Output file prefix |
| `--output-dir` | `.` |  | 输出目录｜Output directory |
| `--no-tsv` | — | store_false | 不输出TSV文件｜Do not output TSV file |
| `--no-csv` | — | store_false | 不输出CSV文件｜Do not output CSV file |
| `--no-excel` | — | store_false | 不输出Excel文件｜Do not output Excel file |
| `--meme-path` | `~/miniforge3/envs/protein/bin/meme` |  | MEME可执行文件路径｜MEME executable path |
| `-protein` | `True` | store_true | 输入序列为蛋白质｜Input sequences are protein (default: True) |
| `-dna` | — | store_true | 输入序列为DNA｜Input sequences are DNA |
| `-mod` | `zoops` | zoops/anr/oor | Motif分布模式｜Motif distribution mode |
| `-nmotifs` | `10` | int | Motif数量｜Number of motifs |
| `-minw` | `6` | int | 最小motif宽度｜Minimum motif width |
| `-maxw` | `50` | int | 最大motif宽度｜Maximum motif width |
| `-objfun` | `classic` | classic/de/ce/cd | 目标函数｜Objective function |
| `-markov-order` | `0` | int | Markov链阶数｜Markov chain order |
| `--no-extract-motifs` | — | store_false | 不提取motif序列｜Do not extract motif sequences |
| `--input-fasta` | — |  | 原始FASTA文件路径（解析模式时用于提取motif序列）｜Original FASTA file path (for motif extraction in parse-only mode) |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- MEME Suite 的 meme（默认 `~/miniforge3/envs/protein/bin/meme`，通过 conda 环境自动包装调用；仅运行模式需要）
- Python 3：pandas、biopython（SeqIO）、openpyxl（Excel）

## 常见问题 | FAQ

**Q1：`-protein` 和 `-dna` 能同时用吗？**
不能，二者互斥。`-protein` 是默认；给 DNA 时用 `-dna`（此时自动关闭 protein）。

**Q2：解析模式提示「找不到输入 FASTA」？**
解析模式要提取 motif 序列片段时需要原始 FASTA。用 `--input-fasta` 显式指定即可；只想要位置表可忽略该警告。

**Q3：MEME 为什么没重跑？**
断点续传按 `meme_results_meme_out/meme.xml` 是否存在判断。想强制重跑，删掉该 xml（或整个 `_meme_out` 目录）再跑。

**Q4：输出目录里为什么多了 `_meme_out` 目录？**
那是 MEME 的原始输出（meme.xml/meme.txt），保留它便于以后用 `--parse-only` 复解析。不想要可直接删。