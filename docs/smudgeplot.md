# Smudgeplot 基因组倍体分析工具 | Smudgeplot Genome Ploidy Analysis Tool

## 简介 | Introduction

Smudgeplot 是一个基于 k-mer 分析评估基因组倍体结构和杂合度特征的工具。该工具使用 FastK k-mer 数据库进行基因组倍体推断，通过分析杂合 k-mer 对的覆盖度模式来识别基因组的倍体结构（如二倍体、三倍体、四倍体等）。

本模块集成了 **FastK** 和 **Smudgeplot** 的完整工作流程，支持从 FASTQ 文件一键运行全部分析步骤。

Smudgeplot is a tool for evaluating genome ploidy structure and heterozygosity based on k-mer analysis. It uses FastK k-mer databases to infer genome ploidy by analyzing coverage patterns of heterozygous k-mer pairs.

This module integrates the complete **FastK** and **Smudgeplot** workflow, supporting one-click execution of all analysis steps from FASTQ files.

## 功能 | Features

- **FastK k-mer 计数**: 从 FASTQ 文件生成 k-mer 数据库
- **k-mer 对提取** (hetmers): 从 FastK k-mer 数据库中提取所有杂合 k-mer 对
- **倍体推断** (plot/all): 基于 k-mer 对覆盖度模式推断基因组倍体
- **可视化**: 生成 smudgeplot 图表展示倍体结构
- **灵活配置**: 支持多种参数自定义和输出格式

## 软件路径 | Software Paths

- **FastK**: `/share/org/YZWL/yzwl_lixg/miniforge3/envs/smudgeplot/bin/FastK`
- **Smudgeplot**: `/share/org/YZWL/yzwl_lixg/miniforge3/envs/smudgeplot/bin/smudgeplot`

## 安装 | Installation

Smudgeplot 已集成在 BioPyTools 中，无需额外安装。

Smudgeplot is integrated into BioPyTools, no additional installation required.

## 使用方法 | Usage

### 方式一：从文件夹自动识别 FASTQ 文件 | Method 1: Auto-detect FASTQ Files from Directory

```bash
# 自动递归搜索文件夹中所有 FASTQ 文件
# Auto-recursively find all FASTQ files in directory
# 支持：.fastq, .fq, .fastq.gz, .fq.gz
biopytools smudgeplot \\
    -i fastq_folder/ \\
    -o results \\
    -n 50
```

**支持的文件格式**：
- 未压缩：`.fastq`, `.fq`
- 压缩：`.fastq.gz`, `.fq.gz`
- 递归搜索子目录

### 方式二：从 FASTQ 文件运行完整流程 | Method 2: Run Complete Workflow from FASTQ Files

```bash
# 双端测序数据
# Paired-end sequencing data
biopytools smudgeplot \\
    -i sample_R1.fastq.gz \\
    -I sample_R2.fastq.gz \\
    -o results \\
    -n 50
```

### 方式三：从已有 FastK ktab 文件运行 | Method 3: Run from Existing FastK ktab File

```bash
# 使用默认参数运行完整分析
# Run complete analysis with default parameters
biopytools smudgeplot -i fastk_table.ktab -o results
```

### 指定单倍体覆盖度 | Specify Haploid Coverage

```bash
# 当已知或已估算单倍体覆盖度时指定
# Specify when haploid coverage is known or estimated
biopytools smudgeplot -i fastk_table.ktab -o results -n 50
```

### 自定义 FastK 参数 | Custom FastK Parameters

```bash
# 从 FASTQ 运行时自定义 FastK 参数
# Custom FastK parameters when running from FASTQ
biopytools smudgeplot \\
    -i sample_R1.fastq.gz \\
    -I sample_R2.fastq.gz \\
    -o results \\
    -k 21 \\
    --fastk-threads 8 \\
    --fastk-memory 32
```

### 自定义 Smudgeplot 参数 | Custom Smudgeplot Parameters

```bash
# 自定义 k-mer 阈值和线程数
# Custom k-mer threshold and thread count
biopytools smudgeplot \\
    -i fastk_table.ktab \\
    -o results \\
    -L 12 \\
    --hetmers-threads 8 \\
    -n 50
```

### 绘图参数 | Plot Parameters

```bash
# 自定义绘图样式
# Custom plot style
biopytools smudgeplot \\
    -i fastk_table.ktab \\
    -o results \\
    -n 50 \\
    --ylim 70 \\
    --col-ramp magma \\
    --format pdf
```

### 分步运行 | Step-by-step Execution

```bash
# 只运行 FastK 步骤 (从 FASTQ 生成 ktab)
# Run only FastK step (generate ktab from FASTQ)
biopytools smudgeplot \\
    -i sample_R1.fastq.gz \\
    -I sample_R2.fastq.gz \\
    -o results \\
    --step fastk

# 只运行 hetmers 步骤
# Run only hetmers step
biopytools smudgeplot -i fastk_table.ktab -o results --step hetmers

# 只运行 plot 步骤 (需要已有 .smu 文件)
# Run only plot step (requires existing .smu file)
biopytools smudgeplot -i fastk_table.ktab -o results -n 50 --step plot
```

## 参数说明 | Parameters

### 必需参数 | Required Parameters

| 参数 | 说明 | 示例 |
|------|------|------|
| `-i, --input` | 输入路径：FastK ktab文件、FASTQ文件或包含FASTQ的目录 | `fastq_folder/` 或 `sample_R1.fq.gz` 或 `fastk_table.ktab` |
| `-I, --input2` | 第二个 FASTQ 文件 (仅单文件输入时使用) | `sample_R2.fastq.gz` |
| `-o, --output` | 输出目录 | `results/` |

**输入说明**：
- **文件夹输入**：自动递归搜索所有 FASTQ 文件（支持 `.fastq`, `.fq`, `.fastq.gz`, `.fq.gz`）
- **单文件输入**：指定单个 FASTQ 文件，可使用 `-I` 添加第二个文件
- **ktab 文件**：直接使用已有的 FastK k-mer 数据库文件

### FastK 选项 | FastK Options

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `-k, --kmer-size` | 31 | k-mer 大小 |
| `--fastk-threads` | 4 | FastK 线程数 |
| `--fastk-memory` | 16 | FastK 内存限制 (GB) |
| `--fastk-path` | 预设路径 | FastK 可执行文件路径 |

### hetmers 选项 | hetmers Options

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `-L, --l-threshold` | None | k-mer 计数阈值 (低于此值视为错误) |
| `--hetmers-threads` | 4 | hetmers 线程数 |

### plot 选项 | plot Options

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `-n, --haploid-coverage` | None | 预期单倍体覆盖度 |
| `--title` | None | 图表标题 |
| `--ylim` | None | Y 轴上限 (覆盖度和) |
| `--col-ramp` | viridis | 调色板 (viridis/magma/mako/grey.colors) |
| `--invert-cols` | False | 反转调色板 |
| `--format` | png | 输出格式 (pdf/png/svg) |
| `--json-report` | False | 生成 JSON 报告 |

### all 命令选项 | all Command Options

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--cov-min` | 6 | 最小覆盖度探索范围 |
| `--cov-max` | 100 | 最大覆盖度探索范围 |
| `--cov` | None | 假定的覆盖度 (不进行 1n 覆盖度推断) |
| `-d, --distance` | None | Manhattan 距离 (用于局部聚合) |

### 其他选项 | Other Options

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--tmp-dir` | /tmp | 临时文件目录 |
| `--smudgeplot-path` | 预设路径 | Smudgeplot 可执行文件路径 |
| `-v, --verbose` | 0 | 详细输出模式 (-v: INFO, -vv: DEBUG) |
| `--quiet` | False | 静默模式 (只输出 ERROR) |
| `--log-file` | smudgeplot.log | 日志文件名 |
| `--log-dir` | ./logs | 日志目录 |
| `--dry-run` | False | 模拟运行 (不实际执行) |

## 输出文件 | Output Files

分析完成后，将在输出目录中生成以下文件：

After analysis completion, the following files will be generated in the output directory:

| 文件 | 说明 |
|------|------|
| `{prefix}_kmerpairs.smf` | k-mer 对二进制文件 |
| `{prefix}_kmerpairs.smu` | k-mer 对汇总文件 |
| `{prefix}_kmerpairs_text.smu` | k-mer 对文本文件 (covA, covB, freq) |
| `{prefix}_log_linear.png` | 线性刻度 smudgeplot |
| `{prefix}_log_log.png` | 对数刻度 smudgeplot |
| `{prefix}_smudge_summary.txt` | 倍体推断摘要 |

## 完整工作流程示例 | Complete Workflow Example

### 方式一：使用 BioPyTools 一键运行完整流程

```bash
# 从 FASTQ 文件一键运行完整流程 (FastK + Smudgeplot)
# Run complete workflow from FASTQ files in one command
biopytools smudgeplot \\
    -i data/sample_R1.fastq.gz \\
    -I data/sample_R2.fastq.gz \\
    -o smudgeplot_results \\
    -n 50 \\
    -L 12 \\
    --ylim 70
```

### 方式二：分步运行

#### 1. 运行 FastK 生成 k-mer 数据库

```bash
# 使用 BioPyTools 运行 FastK
# Run FastK using BioPyTools
biopytools smudgeplot \\
    -i data/sample_R1.fastq.gz \\
    -I data/sample_R2.fastq.gz \\
    -o smudgeplot_results \\
    --step fastk

# 或直接使用 FastK 命令
# Or use FastK command directly
/share/org/YZWL/yzwl_lixg/miniforge3/envs/smudgeplot/bin/FastK \\
    -v -t4 -k31 -M16 -T4 \\
    data/sample_R1.fastq.gz data/sample_R2.fastq.gz \\
    -N smudgeplot_results/FastK_Table
```

#### 2. 运行 Smudgeplot hetmers

```bash
biopytools smudgeplot \\
    -i smudgeplot_results/FastK_Table.ktab \\
    -o smudgeplot_results \\
    --step hetmers \\
    -L 12
```

#### 3. 运行 Smudgeplot plot

```bash
biopytools smudgeplot \\
    -i smudgeplot_results/FastK_Table.ktab \\
    -o smudgeplot_results \\
    -n 50 \\
    --ylim 70 \\
    --step plot
```

### 解读结果

- **smudgeplot 图表**: 每个倍体结构会在图表上形成独特的 "smudge" (斑点)
- **二倍体**: 会在 (0.25, 0.5) 和 (0.5, 0.5) 位置形成 smudge
- **三倍体**: 会在 (0.33, 0.67) 位置形成 smudge
- **四倍体**: 会在多个位置形成 smudge

## 注意事项 | Notes

1. **覆盖度要求**: 建议测序覆盖度至少 30x 以获得可靠的倍体推断结果
2. **k-mer 大小**: 推荐使用 k=21 或 k=31
3. **单倍体覆盖度**: 如果不知道确切的覆盖度，可以先用 GenomeScope 估算
4. **内存需求**: hetmers 步骤可能需要大量内存，建议至少 32GB RAM

## 参考文献 | References

Ranallo-Benavidez, T.R., Jaron, K.S. & Schatz, M.C. GenomeScope 2.0 and Smudgeplot for reference-free profiling of polyploid genomes. *Nature Communications* **11**, 1432 (2020). https://doi.org/10.1038/s41467-020-14998-3

## 官方资源 | Official Resources

- GitHub: https://github.com/KamilSJaron/smudgeplot
- FAQ Wiki: https://github.com/KamilSJaron/smudgeplot/wiki/FAQ
- GenomeScope Web Server: http://genomescope.org/genomescope2.0
