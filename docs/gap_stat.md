# Gap 统计 | Genome Gap Statistics

一句话理解：**数出基因组里有多少段「NNNN…」空洞、每段在哪条序列哪个位置、有多长，还能把一个文件夹里的多个基因组一起统计**。

## 功能概述 | Overview

- 纯 Python 实现，无任何外部软件依赖
- 单文件模式：统计一个 FASTA 里的所有 gap
- 文件夹批量模式：一次统计目录下所有 FASTA，合并输出到单个文件（带 sample 列区分）
- 输出每个 gap 的序列名、起止坐标、长度

## 快速开始 | Quick Start

```bash
biopytools gap-stat -i genome.fa -o gaps.txt
```

输入一个 FASTA（或一个装 FASTA 的文件夹），把 gap 列表写到文件；不写 -o 则打印到终端。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| gap / N 区段 | 组装时没测通、用一串 N 占位的「空洞」 |
| 1-based 坐标 | 从 1 开始数位置（基因组坐标的常见约定），第 1 个碱基是位置 1 |
| 连续 N | 一串紧挨着的 N；--min-n 决定「至少多少个 N 连在一起才算一个 gap」 |

## 输入 | Input

一个 FASTA 文件，或一个包含多个 FASTA 文件的文件夹（批量模式自动识别 .fa/.fasta/.fna 等扩展名）。

```text
>Chr1
ACGTNNNNNACGTNNNTT
>Chr2
NNNNGGGGNN
```

## 参数说明 | Parameters

### 输入输出 | Input & output

**通俗理解|In plain words:** -i 是要统计的文件或文件夹；-o 是结果写到哪个文件（不写就打印到终端）。批量模式下 -o 必须给一个文件路径，所有基因组的 gap 会合并写进去。

### 过滤参数 | Filtering

**通俗理解|In plain words:** --min-n 决定「至少几个 N 连在一起才算 gap」。默认 1 意味着单个 N 也算；如果你只关心大的空洞（比如 >=100 bp），把它调大即可，调大后小串 N 会被忽略。

## 输出 | Output

单文件模式输出一个 TSV 表（4 列）：

```text
seq_name	start	end	gap_length
Chr1	5	9	5
Chr1	14	16	3
Chr2	1	4	4
```

批量模式多一列 sample（样品/文件名），并在文件末尾追加每个样品的汇总统计：

```text
sample	seq_name	start	end	gap_length
genomeA	Chr1	5	9	5
genomeB	Chr1	1	4	4
...
统计信息|Statistics ...
```

此外，汇总统计（总 gap 数、总 gap 长度、含 gap 的序列数）会打印到 stderr（终端）。

## 结果解读 | Interpreting Results

**通俗理解|In plain words:** 每一行就是一个空洞的位置和大小。

- start/end：该 gap 在序列上的起止坐标（1-based，含两端）；
- gap_length：这段空洞多长（bp）；
- **gap 越多、越长，组装连续性越差**；配合 N50、gap 占比一起评估组装质量。

## 参数选择建议 | Parameter Guidance

- 一般统计用默认 --min-n 1；
- 只关心大空洞时把 --min-n 调到 100 或更大，能过滤掉测序/组装产生的零星小 N。

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 输入FASTA文件或文件夹｜Input FASTA file or directory |
| `-o, --output` | — |  | 输出文件路径｜Output file path |
| `--min-n` | `1` | int | 最少N数量｜Minimum consecutive N count |

<!-- END PARAMS:auto -->
## 依赖 | Dependencies

- 无外部软件依赖，纯 Python 标准库实现

## 常见问题 | FAQ

**Q1：支持断点续传吗？**
不支持，也不需要——统计是纯读文件、无中间状态，重跑即得同样结果。

**Q2：-o 不写会怎样？**
gap 列表直接打印到终端（stdout）；汇总统计仍打印到 stderr。想留档就写 -o。

**Q3：文件夹批量模式怎么识别 FASTA？**
按扩展名识别 .fa/.fasta/.FA/.FASTA/.fna/.fnasta/.FNA/.FNASTA，其余文件忽略，识别到的文件按文件名排序后合并统计。

**Q4：能直接测 .gz 压缩的 FASTA 吗？**
不能，需先解压成普通 FASTA。

