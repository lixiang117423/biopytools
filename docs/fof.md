# FOF 文件生成 | FOF File Generation

一句话理解：**扫描一个目录（或单个文件），自动生成一张「样本名 -> 文件绝对路径」的两列映射表（FOF 文件）**，供下游批量流程按样本名去取文件。

## 功能概述 | Overview

- 从文件或目录生成「样品名 -> 绝对路径」的 tab 分隔映射表
- 支持按后缀过滤文件、递归扫描子目录、识别软链接
- 自动从文件名剥离已知后缀得到样本名
- 排除输出文件自身、检测重复样本名并告警
- 纯 Python 实现，无外部依赖

## 快速开始 | Quick Start

```bash
biopytools fof -i ./data/ -o samples.fof
```

最小输入：一个输入目录（或文件）+ 一个输出 FOF 文件路径。

## 零基础概念速览 | Concepts in plain words

| 术语<br>Term | 通俗理解<br>Plain meaning |
|------|----------|
| FOF 文件 | 「File Of Files」的缩写，即一份文件清单，本工具产出的是「样本名 + 路径」两列的表格 |
| 样本名 sample name | 从文件名剥掉后缀得到的前半部分，如 `sample1_R1.fastq.gz` 剥出 `sample1_R1` |
| 后缀 suffix | 文件名末尾的标识，如 `.fastq.gz`，用来圈定要扫描的文件类型 |
| 递归扫描 | 是否连同所有子目录一起扫，还是只扫最外层 |

## 输入 | Input

### 输入文件或目录

- 单个文件：直接作为一条记录
- 目录：扫描其中文件（是否递归由 `-r` 决定）
- 支持软链接（按链接目标判断是否为文件）

## 参数说明 | Parameters

### 必需参数 | Required

**通俗理解|In plain words:** `-i` 指定「扫哪里」（文件或目录），`-o` 指定「映射表写到哪」。输出文件会自动排除在扫描结果之外，避免把自己写进清单。

### 后缀过滤 | Suffix filter

**通俗理解|In plain words:** `-s` 只保留指定后缀的文件，可多次指定（如 `-s .fastq.gz -s .fq.gz`）。不写 `-s` 就扫目录里所有文件。后缀必须以 `.` 开头，否则报错。

### 递归扫描 | Recursive

**通俗理解|In plain words:** `-r` 是否连子目录一起扫。目录层级多、文件散落在子目录时加 `-r`；只想扫最外层一层就不用加。

## 分析流程 | Pipeline

```text
输入文件/目录 + 后缀/递归选项
  |
  v
1. 收集符合条件的所有文件(排除输出文件自身,排序)
  |
  v
2. 从每个文件名剥离后缀得到样本名
  |
  v
3. 写出两列 FOF 文件(样本名 	 绝对路径)
  |
  v
4. 检测重复样本名并告警
```

## 输出 | Output

```text
samples.fof      # tab 分隔的映射表,每行: 样本名 	 文件绝对路径
fof.log          # 运行日志(位于输出文件的父目录)
```

FOF 文件内容示例：

```text
sample1	/abs/path/sample1.fastq.gz
sample2	/abs/path/sample2.fastq.gz
```

- 两列：第一列样本名，第二列文件的绝对路径（软链接会解析成真实路径）
- 样本名提取规则：优先按你给的后缀精确剥离（最长匹配）；否则先剥压缩后缀（.gz 等）再剥一层已知数据扩展（.fastq/.bam/.vcf 等）；兜底剥最后一层扩展名
- 结果为空时仍会写一个空的 FOF 文件（0 行）

## 结果解读 | Interpreting Results

**通俗理解|In plain words:** 打开 FOF 文件看两点：行数是否等于你预期的文件数、每行两列是否完整、样本名是否符合预期。

- **行数**：应等于被扫描到的文件数。明显偏少说明后缀过滤太严或目录不对
- **样本名**：检查剥离是否干净（比如 `sample1_R1.fastq.gz` 若只想要 `sample1`，需用 `-s _R1.fastq.gz` 这类后缀精确剥离）
- **重复样本名告警**：日志里若提示「检测到重复样品名」，说明多个文件剥出了同一个名字，下游按名取文件会歧义，需检查命名
- **好坏判据**：行数符合预期、无重复样本名告警、路径真实存在，即为正常

## 参数选择建议 | Parameter Guidance

- **目录里只有一种文件**：直接 `biopytools fof -i dir -o out.fof`
- **目录里混了多种文件**：用 `-s` 精确圈定（如只要 FASTQ 就 `-s .fastq.gz`）
- **文件散落在子目录**：加 `-r` 递归
- **想要更干净的样本名**：给 `-s` 一个更长的后缀，如 `-s _R1.fastq.gz`，剥出的名字就更短

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 输入文件或目录｜Input file or directory |
| `-o, --output` | 必填 | Path | 输出FOF文件路径｜Output FOF file path |
| `-s, --suffix` | — |  | 文件后缀过滤（可多次指定，如 -s .fastq.gz -s .fq.gz）｜File suffix filter (can be specified multiple times) |
| `-r, --recursive` | `False` |  | 递归扫描子目录｜Recursively scan subdirectories |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 输入文件或目录｜Input file or directory |
| `-o, --output` | 必填 |  | 输出FOF文件路径｜Output FOF file path |
| `-s, --suffix` | — | append | 文件后缀过滤（可多次指定，如 -s .fastq.gz -s .fq.gz）｜File suffix filter (can be specified multiple times) |
| `-r, --recursive` | `False` | store_true | 递归扫描子目录｜Recursively scan subdirectories |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- 仅 Python 标准库（pathlib / os），无外部软件、无 conda 环境名

## 常见问题 | FAQ

**Q1：后缀为什么要以 . 开头？**
程序会校验后缀必须以 `.` 开头，否则报错。这是为了明确「这是文件后缀」而不是任意子串。

**Q2：输出文件会不会把自己也写进清单？**
不会。程序会解析输出文件的真实路径并从扫描结果里排除，避免自引用。

**Q3：软链接能识别吗？**
能。输入路径和扫描到的文件都支持软链接，按链接目标判断是否为文件。

**Q4：样本名剥不干净怎么办？**
样本名提取有优先级：先按 `-s` 后缀精确剥离（最长匹配），所以给个足够长的后缀最可控；不写 `-s` 时用内置的「压缩后缀 + 数据扩展」规则剥一层，遇到自定义命名可能剥不干净，此时建议显式指定 `-s`。
