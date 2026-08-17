# Bismark 全基因组甲基化分析 | Bismark Whole-Genome Methylation Analysis

一句话理解：**把重亚硫酸盐测序（WGBS）读段比对到基因组，逐位点读出 DNA 上哪些 C 被加了甲基化修饰**，解决「想知道基因组哪些位置发生了甲基化」的问题。

## 功能概述 | Overview

- 基于 Bismark（Bowtie2 后端）的 WGBS 全流程：建索引 → 比对 → 甲基化提取 → 结果拆分汇总
- 自动为每个样品配对 R1/R2 读段（按后缀模式），支持多样品批量处理
- 提取 CG、CHG、CHH 三种甲基化背景，输出可读的 CX_report 和按背景拆分的文本
- 比对用 `--non_directional` 非定向模式，默认排除重叠读段（`--no_overlap`）
- 断点续传：索引已建、某样品最终结果已存在时自动跳过
- 每个样品独立处理，单样失败不影响其他样品

## 快速开始 | Quick Start

```bash
biopytools bismark -g genome.fa -i ./raw_data -o ./bismark_results
```

最小输入：一个基因组 FASTA、一个存放配对 FASTQ 的目录、一个输出目录。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| 重亚硫酸盐测序（WGBS） | 一种专门测「DNA 甲基化」的技术：未甲基化的 C 会变成 T，甲基化的 C 保持不变 |
| 甲基化 | DNA 上的一种化学标记，像「开关标签」，常调控基因表达 |
| CpG 背景（CG） | C 后面紧跟 G 的位点，动植物里最常见的甲基化位点 |
| CHG / CHH | C 后面跟 H（非 G）的位点，植物里常见 |
| CX_report | 汇总每个 C 位点甲基化状态的报告文件 |
| 非定向模式 | 不管读段方向，两头都能比对；适合全基因组建库 |
| 重叠读段 | 一对读段（R1/R2）中间重叠的部分，默认只算一次（避免重复计数） |

## 输入 | Input

### 基因组

FASTA 格式（.fa / .fasta）。Bismark 索引会建在**基因组文件所在目录**下的 `Bisulfite_Genome/` 子目录里。

### 原始 FASTQ 目录

存放配对读段的目录。R1 文件必须匹配 `*<pattern>`（默认 `*_1_clean.fq.gz`），R2 由程序自动推断（把模式里的 `1` 换成 `2`，或 `R1` 换成 `R2`）。

```text
raw_data/
├── sampleA_1_clean.fq.gz
├── sampleA_2_clean.fq.gz
├── sampleB_1_clean.fq.gz
└── sampleB_2_clean.fq.gz
```

## 参数说明 | Parameters

### 必需参数 | Required

**通俗理解|In plain words:** 三个参数分别是「基因组、原始数据目录、结果放哪」。输出目录必须能写入，程序会先做可写性检查。

### 读段匹配 | Read matching

**通俗理解|In plain words:** `--pattern` 决定怎么认出「这是 R1 文件」。默认 `_1_clean.fq.gz` 对应「清洗后的 R1」。如果你的文件名是别的后缀（如 `_R1.fastq.gz`），改成对应模式；程序会据此自动找配对的 R2。

### 资源与线程 | Resources

**通俗理解|In plain words:** 线程数只影响速度，不影响结果。注意比对阶段实际用的是「线程数除以 4」跑 `--multicore`，提取阶段用全部线程。`--sort-buffer` 是甲基化提取排序时的内存缓冲，大数据不够时可调大，一般用默认 400G。

### 重叠读段 | Overlap

**通俗理解|In plain words:** 默认排除 R1/R2 重叠部分的重复计数（更严谨）。加 `--include-overlap` 才保留重叠读段，通常不需要动。

## 分析流程 | Pipeline

```text
输入基因组FASTA + 配对FASTQ目录
  -> Step 1: 构建Bismark索引（bismark_genome_preparation，若已存在则跳过）
  -> Step 2: 逐样品比对（bismark --bowtie2 --non_directional --multicore）
  -> Step 3: 甲基化提取（bismark_methylation_extractor --paired-end --CX_context）
  -> Step 4: 整理结果、按CG/CHG/CHH拆分CX_report
  -> Step 5: 生成总结报告 bismark_pipeline_summary.txt
```

## 输出 | Output

```text
bismark_results/
├── bismark_pipeline_summary.txt   # 运行总结报告
├── mapping/
│   └── <sample>_bismark_bt2_pe.bam  # 各样品比对BAM
├── result/
│   └── <sample>/
│       ├── <...>.CX_report.txt      # 全背景甲基化报告（核心）
│       ├── <...>.CpG.txt            # 按CG背景拆分
│       ├── <...>.CHG.txt            # 按CHG背景拆分
│       ├── <...>.CHH.txt            # 按CHH背景拆分
│       ├── <...>.bedGraph           # 甲基化bedGraph
│       ├── <...>_splitting_report.txt  # 拆分统计报告
│       └── <...>M-bias.txt          # M-bias 报告
├── tmp/
│   └── <sample>/                    # 中间文件
└── logs/
    └── bismark_pipeline_<时间戳>.log
```

关键文件说明：

- **CX_report.txt**：每行一个 C 位点，含染色体、位置、正负链、甲基化/未甲基化计数、背景（CG/CHG/CHH），是下游分析的核心输入。
- **CpG/CHG/CHH.txt**：从 CX_report 按背景拆分的文本，直接拿去做差异甲基化分析。
- **M-bias.txt**：反映「读段两端甲基化是否偏高」，用于判断建库质量。

## 结果解读 | Interpreting Results

**通俗理解|In plain words:** 看甲基化结果好坏，先看比对报告里的比对率，再看 M-bias 是否正常。

- **CX_report.txt**：第三、四列通常是「甲基化 C 数 / 未甲基化 C 数」，两者相加就是这个位点的覆盖度；覆盖度很低的位点结论不可靠。
- **CpG 甲基化水平**：动植物基因组里 CpG 甲基化水平通常较高（动物常在 70%~90%），若异常接近 0 或 100%，检查建库是否失败。
- **M-bias.txt**：正常样本读段两端甲基化略高，但不应出现一端极高、另一端极低；M-bias 曲线严重倾斜提示建库或酶切问题。
- **splitting_report.txt**：三种背景各拆出多少行，行数为 0 的背景说明该背景没有覆盖（植物 CHG/CHH 应有数据，动物可能极少）。

## 参数选择建议 | Parameter Guidance

- **--pattern**：按实际文件名改，务必让 R1 模式里含 `1` 或 `R1`，否则程序无法推断 R2。
- **-t 线程数**：建议给机器核数；比对与提取都吃 CPU，线程越多越快。
- **--sort-buffer**：处理大样本或高覆盖数据时若报内存/排序错，可调小（如 100G）以适配机器内存。
- **--include-overlap**：一般不需要；只有当你明确要统计重叠区域的甲基化时才开。

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--version` | — |  | 显示版本信息｜Show version information |
| `--genome, -g` | 必填 |  | 基因组FASTA文件｜Genome FASTA file |
| `--input, -i` | 必填 |  | 原始FASTQ数据目录｜Raw FASTQ data directory |
| `--output-dir, -o` | 必填 | Path | 主输出目录｜Main output directory |
| `--pattern, -p` | `_1_clean.fq.gz` |  | R1文件后缀模式｜Suffix pattern for R1 files |
| `--threads, -t` | `12` | int | 线程数｜Number of threads |
| `--sort-buffer` | `400G` | str | 排序缓冲区大小｜Sort buffer size |
| `--include-overlap` | — |  | 包含重叠读段｜Include overlapping reads |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-g, --genome` | 必填 |  | 基因组FASTA文件路径｜Genome FASTA file |
| `-i, --input` | 必填 |  | 原始FASTQ数据目录｜Raw FASTQ data directory |
| `-o, --output-dir` | 必填 |  | 主输出目录｜Main output directory |
| `-p, --pattern` | `_1_clean.fq.gz` |  | R1文件的后缀模式｜Suffix pattern for R1 files |
| `-t, --threads` | `88` | int | 使用的线程数｜Number of threads |
| `--sort-buffer` | `400G` | str | 提取步骤的排序缓存大小｜Sort buffer size |
| `--no-no-overlap` | — | store_false | 不忽略重叠的reads｜Do NOT ignore overlapping reads |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- bismark（需在 PATH 中）
- bismark_genome_preparation（需在 PATH 中）
- bowtie2（需在 PATH 中，索引构建时自动定位其目录）
- bismark_methylation_extractor（需在 PATH 中）
- 无内置 conda 环境名，所有 Bismark 相关命令依赖调用环境 PATH

## 常见问题 | FAQ

**Q1：断点续传怎么判断？换数据重跑会不会跳过？**
索引看 `<基因组目录>/Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa` 是否存在；每个样品看 `result/<sample>/*CX_report.txt` 是否存在。若换了数据或参数想重跑某样品，先删掉对应样品的 result 目录。

**Q2：提示「未找到配对的 FASTQ 文件」？**
说明在原始数据目录里没找到匹配 `*<pattern>` 的 R1 文件，或找到了 R1 但缺对应 R2。检查文件命名是否与 --pattern 一致，以及 R1/R2 是否成对存在。

**Q3：--pattern 必须以 .gz 结尾吗？**
是的。校验要求模式以 `.gz` 结尾（程序按压缩 FASTQ 处理）。若你的数据是未压缩的 .fastq，需先压缩或改文件名。

**Q4：比对率偏低怎么办？**
先看日志里 bismark 的比对报告；常见原因是非定向数据用了定向参数、或基因组与样本物种不匹配。程序已固定用 `--non_directional`，确认建库类型即可。

