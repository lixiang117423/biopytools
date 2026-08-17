# Protein Stats - 蛋白质理化性质分析 | Protein Stats Physicochemical Properties

一句话理解：**给每条蛋白算一堆「理化指标」——多长、多重、等电点、疏水性等**。输入一个蛋白质 FASTA，输出一张表格，每条蛋白一行、每项指标一列，方便做蛋白组整体统计。

## 功能概述 | Overview

- 基于 Biopython ProtParam 计算蛋白质理化性质
- 默认计算：序列长度、分子量、等电点(pI)及酸碱性判断
- 可选计算：氨基酸组成、不稳定指数、脂肪指数(GRAVY 疏水性)、芳香性
- 输出 TSV / CSV / Excel 三种格式，可按需关闭基础项（`--no-length` 等）
- 纯 Python 计算，无外部命令行调用，无需 conda 环境包装

## 快速开始 | Quick Start

```bash
biopytools protein-stats -i proteins.fa
```

最小输入：一个蛋白质 FASTA 文件。结果默认写到当前目录的 `protein_stats.tsv`。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| 分子量(MW) | 一个蛋白分子有多「重」，单位道尔顿(Da)，越大分子越大 |
| 等电点(pI) | 蛋白「不带电」的 pH 值：pI<7 偏酸性，>7 偏碱性，用于选择纯化/分离条件 |
| 不稳定指数 | 衡量蛋白在试管里「容不容易降解」，>40 通常判为不稳定 |
| GRAVY | 整体疏水性：正值越疏水（亲油），负值越亲水 |
| 芳香性 | 含苯环类氨基酸(色氨酸、酪氨酸、苯丙氨酸)的比例，与紫外吸收有关 |

## 输入 | Input

蛋白质 FASTA 文件：

```text
>protein_01
MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRV...
>protein_02
MSSVPTGIPQETGFSTGGH...
```

- 序列必须是氨基酸（蛋白质）；含 B/J/O/U/X/Z 等模糊氨基酸时会先剔除它们再算分子量/等电点，长度仍按原始序列计

## 参数说明 | Parameters

### 必需参数 | Required

**通俗理解|In plain words:** 唯一必填是输入蛋白 FASTA 文件。

### 输出配置 | Output

**通俗理解|In plain words:** 决定结果写到哪、用什么格式。默认 `protein_stats.tsv`，格式 tsv；Excel 需要 openpyxl。

### 基础计算项（默认开）| Basic metrics

**通俗理解|In plain words:** 长度、分子量、等电点默认都算，**一般不用动**；只想算某几项可用 `--no-length`/`--no-mw`/`--no-pi` 关掉对应的列。

### 高级计算项（默认关）| Advanced metrics

**通俗理解|In plain words:** 氨基酸组成、不稳定指数、GRAVY、芳香性默认不计算（省时）。需要哪个就用对应开关开启，**不需要就保持默认**。

## 输出 | Output

```text
.
└── protein_stats.tsv        # 结果表格(每蛋白一行)
    # (另生成 protein_stats_analysis.log 运行日志于当前目录)
```

- 输出是**单个文件**（路径由 `-o` 指定），不是目录
- 默认列：`id / length / mw_da / pi / acid_base / 酸碱性`；开启高级项后追加 `aa_*`、`instability_index`、`gravy`、`aromaticity` 等列

## 结果解读 | Interpreting Results

**通俗理解|In plain words:** 一行一个蛋白，每个指标单独看含义即可；做整组统计时可对每列求平均/分布。

- `length`：氨基酸残基数；`mw_da`：分子量（Da）
- `pi`：等电点；`acid_base`/`酸碱性`：pI<7 酸性、>7 碱性、=7 中性
- `instability_index`：>40 判为不稳定蛋白（经验阈值），<40 相对稳定
- `gravy`：正值疏水（常见于膜蛋白），负值亲水（常见于可溶性蛋白）
- `aromaticity`：芳香族氨基酸占比，0–1 之间的小数

## 参数选择建议 | Parameter Guidance

- **只看基础指标**：用默认即可，最快
- **`--gravy`**：做「膜蛋白 vs 可溶性蛋白」的粗略区分时开启，看疏水性分布
- **`--instability-index`**：做蛋白表达/稳定性评估时开启
- **`--aa-composition`**：需要每种氨基酸占比时开启（会加 20 列左右）
- **`--output-format excel`**：要在 Excel 里做筛选/排序时用，需 openpyxl

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --protein-fasta` | 必填 |  | 蛋白序列FASTA文件｜Protein sequence FASTA file |
| `-o, --output-file` | `protein_stats.tsv` |  | 输出文件路径｜Output file path |
| `--output-format` | `tsv` | tsv/csv/excel | 输出文件格式｜Output file format |
| `--no-length` | — |  | 不计算序列长度｜Do not calculate sequence length |
| `--no-mw` | — |  | 不计算分子量｜Do not calculate molecular weight |
| `--no-pi` | — |  | 不计算等电点｜Do not calculate isoelectric point |
| `--aa-composition` | — |  | 计算氨基酸组成｜Calculate amino acid composition |
| `--instability-index` | — |  | 计算不稳定指数｜Calculate instability index |
| `--gravy` | — |  | 计算脂肪指数(疏水性)｜Calculate gravy (hydropathy) |
| `--aromaticity` | — |  | 计算芳香性｜Calculate aromaticity |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --protein-fasta` | 必填 |  | 蛋白序列FASTA文件｜Protein sequence FASTA file |
| `-o, --output-file` | `protein_stats.tsv` |  | 输出文件路径｜Output file path |
| `--output-format` | `tsv` | tsv/csv/excel | 输出文件格式｜Output file format |
| `--no-length` | — | store_false | 不计算序列长度｜Do not calculate sequence length |
| `--no-mw` | — | store_false | 不计算分子量｜Do not calculate molecular weight |
| `--no-pi` | — | store_false | 不计算等电点｜Do not calculate isoelectric point |
| `--aa-composition` | — | store_true | 计算氨基酸组成｜Calculate amino acid composition |
| `--instability-index` | — | store_true | 计算不稳定指数｜Calculate instability index |
| `--gravy` | — | store_true | 计算脂肪指数(疏水性)｜Calculate gravy (hydropathy) |
| `--aromaticity` | — | store_true | 计算芳香性｜Calculate aromaticity |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- Python 3：biopython（SeqIO + ProtParam）、pandas
- openpyxl（仅 `--output-format excel` 时需要）

## 常见问题 | FAQ

**Q1：为什么 pI/分子量和预期差一点？**
序列里的模糊氨基酸（B/J/O/U/X/Z）会被剔除后再算分子量和 pI，可能导致与包含模糊残基的完整序列略有差异；长度仍按原始序列计。

**Q2：能一次算多个文件吗？**
一次只处理一个 FASTA。多个文件可写循环逐文件调用，或用 `-o` 分别指定输出文件名。

**Q3：本工具有断点续传吗？**
没有。它是单步纯计算、速度很快，重跑会直接覆盖输出文件，无需续传。

**Q4：输出文件里没有日志目录？**
本工具不产生 `99_logs/` 目录，运行日志 `protein_stats_analysis.log` 写在当前工作目录。