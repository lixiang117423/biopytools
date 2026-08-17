# Resistify NLR 分析 | Resistify NLR Analysis

一句话理解：**运行 Resistify 软件识别植物抗病基因(NLR)，把它的原始结果解析成干净整齐的 TSV/CSV/Excel 汇总表，还能按分类/长度筛选，并按需提取 NLR 和 NB-ARC 序列。**

## 功能概述 | Overview { #overview }

- 运行 `resistify nlr` 识别 NLR(植物细胞内免疫受体)基因
- 解析 `results.tsv` + `domains.tsv`，整合结构域信息(TIR / NB-ARC / CC / LRR / RPW8)
- 支持单文件或多样本目录批量处理(每样本一个子目录)
- 按分类(如 TN / CNL / NL)、序列长度、LRR 长度筛选
- 可选提取 NLR 全序列和 NB-ARC 结构域序列
- 输出 TSV / CSV / Excel 三种格式

## 快速开始 | Quick Start { #quick-start }

```bash
biopytools resistify -i proteins.pep.fa -o ./resistify_output
```

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语 | 通俗理解 |
|------|----------|
| NLR | 植物细胞内的免疫受体，抗病关键，像「细胞内的哨兵」 |
| NB-ARC | NLR 的核心结构域，几乎所有 NLR 都有 |
| TIR / CC | NLR N 端的两种结构域类型，像「哨兵的手」 |
| LRR | 富亮氨酸重复，负责识别病原，像「哨兵的眼睛」 |
| RPW8 | 另一类 N 端结构域 |
| 分类 TN / CNL / NL | 按结构域组合给 NLR 起的名(T=TIR，C=CC，N=NB-ARC，L=LRR) |

## 输入 | Input { #input }

一个蛋白质 FASTA 文件，或一个装 FASTA 文件的目录(自动识别 `.fa` `.fasta` `.faa` `.pep` `.fna`)：

- 单文件：结果直接写到 `--output-dir` 下
- 目录：每个 FASTA 文件视为一个样本，Resistify 结果写到 `<输出目录>/<样本名>/`，汇总表合并到输出目录根

```text
>gene1
MASSS...
>gene2
MAAQ...
```

## 参数说明 | Parameters { #parameters }

### 输入输出 | Input & output { #parameters-io }

**通俗理解|In plain words:** `-i` 给蛋白 FASTA(文件或目录)，`-o` 给输出目录。`--output-prefix` 决定解析结果文件的名字前缀(默认 resistify_results)。

### Resistify 运行 | Resistify run { #parameters-resistify }

**通俗理解|In plain words:** `--resistify-path` 是 Resistify 程序位置，默认走 conda 环境 protein，一般不用动。`--skip-resistify` 跳过运行、只解析已有结果，适合在 Resistify 已经跑完时复用。`--threads` 控制并行度。

### 输出格式 | Output format { #parameters-output }

**通俗理解|In plain words:** 默认 TSV + CSV + Excel 三种都出，`--no-tsv` / `--no-csv` / `--no-excel` 分别关掉。不需要 Excel 就关掉，能省一点时间(还要装 openpyxl)。

### 序列提取 | Sequence extraction { #parameters-extract }

**通俗理解|In plain words:** `--extract-nlr` 提取筛选后每条 NLR 的全序列，`--extract-nbarc` 提取 NB-ARC 结构域序列。默认都关，需要序列时才开。

### 筛选 | Filtering { #parameters-filter }

**通俗理解|In plain words:** `--filter-classification` 按分类名包含匹配(如 TN、CNL、NL)；`--min-length` / `--max-length` 按全长过滤；`--min-lrr-length` 按 LRR 长度过滤。都是「缩小候选集」用，默认不过滤。

## 分析流程 | Pipeline { #pipeline }

```text
输入蛋白 FASTA
    │
    ▼
resistify nlr  ──►  results.tsv / domains.tsv / nlr.fasta / nbarc.fasta
    │
    ▼
解析并整合结构域信息(has_TIR / has_NB_ARC / has_CC / has_LRR / has_RPW8)
    │
    ▼
按分类 / 长度 / LRR 长度筛选
    │
    ▼
提取 NLR / NB-ARC 序列(可选)
    │
    ▼
保存 TSV / CSV / Excel 汇总
```

## 输出 | Output { #output }

```text
resistify_output/
├── results.tsv                     # Resistify 原始结果(每条 NLR 一行)
├── domains.tsv                     # Resistify 结构域结果
├── nlr.fasta                       # Resistify 输出的 NLR 序列
├── nbarc.fasta                     # Resistify 输出的 NB-ARC 序列
├── resistify_results.tsv           # 解析整合后的汇总(核心结果)
├── resistify_results.csv
├── resistify_results.xlsx
├── resistify_results_nlr.fasta     # 仅 --extract-nlr
├── resistify_results_nbarc.fasta   # 仅 --extract-nbarc
└── 99_logs/resistify.log           # 运行日志
```

汇总表关键列：`Sequence`(序列 ID)、`Length`(全长)、`LRR_Length`、`Domains`、`Classification`、`NBARC_motifs`、`has_TIR`、`has_NB_ARC`、`has_CC`、`has_LRR`(批量模式额外有 `Sample` 列)。

## 结果解读 | Interpreting Results { #interpreting-results }

- `resistify_results.tsv` 是最核心的汇总表：一行一条 NLR，`Classification` 是结构域组合名(TN / CNL / NL 等)
- `has_TIR` / `has_CC` / `has_NB_ARC` / `has_LRR` 表示这条 NLR 是否含对应结构域，可用来快速分类统计
- 配合 `--filter-classification CNL` 等可只保留某类 NLR 做下游分析
- 提取的 `_nlr.fasta` / `_nbarc.fasta` 可直接用于建树、多序列比对等
- 若汇总表为空：检查输入是不是蛋白 FASTA(不是核酸)、Resistify 是否正常跑出 results.tsv

## 参数选择建议 | Parameter Guidance { #parameter-guidance }

- 常规「跑一遍拿全表」：默认即可，直接看 resistify_results.tsv / xlsx
- 只要某类 NLR(如 CNL)：`--filter-classification CNL`
- 只要完整可信的 NLR：`--min-length 500` 左右过滤掉碎片
- 需要下游序列分析：加 `--extract-nlr`(和/或 `--extract-nbarc`)
- 已有 Resistify 结果只想重新汇总/筛选：`--skip-resistify -i <含 results.tsv 的目录>`

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 输入蛋白质FASTA文件或目录｜Input protein FASTA file or directory |
| `-o, --output-dir` | `./resistify_output` |  | 输出目录｜Output directory |
| `-t, --threads` | `12` | int | 并行线程数｜Number of parallel threads |
| `--resistify-path` | — |  | Resistify可执行文件路径｜Resistify executable path |
| `--skip-resistify` | — |  | 跳过Resistify运行，仅解析已有结果｜Skip Resistify run, only parse existing results |
| `--output-prefix` | `resistify_results` |  | 解析结果文件前缀｜Parsed results file prefix |
| `--no-tsv` | — |  | 不输出TSV文件｜Do not output TSV file |
| `--no-csv` | — |  | 不输出CSV文件｜Do not output CSV file |
| `--no-excel` | — |  | 不输出Excel文件｜Do not output Excel file |
| `--extract-nlr` | — |  | 提取NLR序列｜Extract NLR sequences |
| `--extract-nbarc` | — |  | 提取NB-ARC序列｜Extract NB-ARC sequences |
| `--filter-classification` | — |  | 按分类筛选(如TN, CNL, NL等)｜Filter by classification (e.g., TN, CNL, NL) |
| `--min-length` | — | int | 最小序列长度｜Minimum sequence length |
| `--max-length` | — | int | 最大序列长度｜Maximum sequence length |
| `--min-lrr-length` | — | int | 最小LRR长度｜Minimum LRR length |
| `--include-motifs` | — |  | 包含motifs详情｜Include motifs details |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 输入蛋白质FASTA文件或目录｜Input protein FASTA file or directory |
| `-o, --output-dir` | `./resistify_output` |  | 输出目录｜Output directory |
| `-t, --threads` | `12` | int | 并行线程数｜Number of parallel threads |
| `--resistify-path` | — |  | Resistify可执行文件路径｜Resistify executable path |
| `--skip-resistify` | — | store_true | 跳过Resistify运行，仅解析已有结果｜Skip Resistify run, only parse existing results |
| `--output-prefix` | `resistify_results` |  | 解析结果文件前缀｜Parsed results file prefix |
| `--no-tsv` | — | store_false | 不输出TSV文件｜Do not output TSV file |
| `--no-csv` | — | store_false | 不输出CSV文件｜Do not output CSV file |
| `--no-excel` | — | store_false | 不输出Excel文件｜Do not output Excel file |
| `--extract-nlr` | — | store_true | 提取NLR序列｜Extract NLR sequences |
| `--extract-nbarc` | — | store_true | 提取NB-ARC序列｜Extract NB-ARC sequences |
| `--filter-classification` | — |  | 按分类筛选(如TN, CNL, NL等)｜Filter by classification (e.g., TN, CNL, NL) |
| `--min-length` | — | int | 最小序列长度｜Minimum sequence length |
| `--max-length` | — | int | 最大序列长度｜Maximum sequence length |
| `--min-lrr-length` | — | int | 最小LRR长度｜Minimum LRR length |
| `--include-motifs` | — | store_true | 包含motifs详情｜Include motifs details |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies { #dependencies }

- resistify：`~/miniforge3/envs/protein/bin/resistify`(conda 环境 protein)
- Python 3 + pandas + openpyxl(Excel 输出)+ biopython(序列提取)

## 常见问题 | FAQ { #faq }

**Q1：支持断点续传吗？**
没有自动断点续传。但 `--skip-resistify` 可复用已有 Resistify 结果，只重跑解析/筛选/提取部分。

**Q2：`--skip-resistify` 需要什么前提？**
输入目录里要有 Resistify 生成的 `results.tsv` 和 `domains.tsv`；单文件模式下 `-i` 指向的就是含这些文件的目录。

**Q3：`--include-motifs` 为什么没效果？**
这是预留参数，motifs.tsv 的解析目前尚未实现，传了也不会多出结果列，属正常现象。

**Q4：Excel 输出失败？**
需要 openpyxl。没装时会跳过 Excel(只出 TSV/CSV)；也可用 `--no-excel` 直接关掉。

**Q5：批量模式结果在哪？**
每个样本的 Resistify 原始结果在 `<输出目录>/<样本名>/`，整合后的汇总表在输出目录根(resistify_results.tsv 等，带 Sample 列区分来源)。