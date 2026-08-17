# 基因组挂载率统计 | Genome Mount Rate

一句话理解：**算出基因组里「最长的 N 条序列」加起来占总长度的百分之几，用来衡量组装把序列「挂载」到染色体的完成度**。

## 功能概述 | Overview

- 纯 Python 实现，无外部依赖
- 计算前 N 条（或最长 N 条）序列占总基因组长度的百分比
- 支持按文件原顺序或按长度降序两种统计口径
- 结果直接打印到终端

## 快速开始 | Quick Start

```bash
biopytools genome-mount-rate -i genome.fa -n 10 --sort
```

输入一个 FASTA，-n 指定统计几条，--sort 表示按长度取最长 N 条。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| 挂载率(mount rate) | 「最长 N 条序列」占全基因组的比例；比例越高说明序列越集中、越接近染色体级 |
| 序列(sequence) | FASTA 里一条以 > 开头的记录（可能是染色体、scaffold 或 contig） |
| 前 N 条 vs 最长 N 条 | 前者按文件里的出现顺序取，后者先按长度排序再取；通常「最长 N 条」才有意义 |

## 输入 | Input

一个 FASTA 文件（.fa / .fasta）。

```text
>Chr1
ACGT...（很长）
>Chr2
ACGT...
>contig1
ACGT...（很短）
```

## 参数说明 | Parameters

### 输入与数量 | Input & number

**通俗理解|In plain words:** -i 是要统计的 FASTA；-n 是「取几条序列来算占比」。-n 常设为该物种的染色体数（如人 24、水稻 12）。

### 排序开关 | Sort

**通俗理解|In plain words:** --sort 决定「取哪 N 条」。**强烈建议加上**——加了之后先按长度从大到小排序，取的是「最长 N 条」（通常对应染色体）；不加则按文件里的原始顺序取前 N 条，如果文件里小 contig 排在前面，算出来的挂载率会严重偏低、没有意义。

## 输出 | Output

结果直接打印到终端（不写文件）：

```text
----------------------------------------
总序列数|Total sequences: 1352
总基因组大小|Total genome size: 812,345,678 bp
统计目标|Target: 最长|longest 12 条序列|sequences
目标序列总长|Target sequences total length: 790,000,000 bp
----------------------------------------
占比|Mount rate: 97.25%
----------------------------------------
```

## 结果解读 | Interpreting Results

**通俗理解|In plain words:** 最后一行「占比|Mount rate」就是要的数。

- 占比接近 90-100%：大部分基因组都集中到了少数几条长序列上，挂载良好；
- 占比偏低（如 <80%）：大量序列散落成小 contig，组装连续性差，需要继续挂载（如 Hi-C）或检查组装。

## 参数选择建议 | Parameter Guidance

- 始终加 --sort，取「最长 N 条」；
- -n 用该物种的染色体数（或预期染色体数）；
- 也可用 -n 100 这类较大值，看「最长 100 条」能覆盖多少，评估 scaffold 挂载情况。

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input, -i` | 必填 |  | 输入FASTA文件｜Input FASTA file path |
| `--number, -n` | 必填 | int | 序列数量｜Number of sequences |
| `--sort` | — |  | 按长度从大到小排序(计算最长N条)｜Sort by length descending (calculate longest N) |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 输入的FASTA文件路径｜Input FASTA file path |
| `-n, --number` | 必填 | int | 要计算的前N条序列数量｜Number of top N sequences to calculate |
| `--sort` | — | store_true | [推荐]开启此选项后，会先按长度从大到小排序，再计算前N条（即计算最长N条的占比）｜[Recommended] Sort sequences by length (descending) before calculating (i.e., calculate longest N) |

<!-- END PARAMS:auto -->
## 依赖 | Dependencies

- 无外部软件依赖，纯 Python 标准库实现

## 常见问题 | FAQ

**Q1：为什么加了 --sort 和没加结果差很多？**
没加 --sort 时按文件里序列的**原始顺序**取前 N 条；若小 contig 排前面，取到的就是一堆短序列，占比自然很低。加 --sort 后先按长度排序，取的是真正最长的 N 条，才有参考意义。

**Q2：结果写到哪里？**
直接打印到终端（stdout），不生成文件；需要留档可重定向，如 biopytools genome-mount-rate -i g.fa -n 10 --sort > result.txt。

**Q3：支持断点续传吗？**
不需要也不支持——纯读文件即时计算，无中间状态。

**Q4：-n 超过了序列总数会怎样？**
不会报错，程序会自动取 min(-n, 总序列数)，即最多取全部序列（此时占比恒为 100%）。

