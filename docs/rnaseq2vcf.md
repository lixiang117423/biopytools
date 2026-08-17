# 转录组变异检测(到 VCF) | RNA-seq Variant Calling (to VCF)

一句话理解：**从 RNA-seq 测序读段里找出基因上的变异位点（SNP/INDEL），走 HISAT2 比对 + GATK 多样本联合检测，最终产出一个过滤后的多样本 VCF 文件**；注释（ANNOVAR）需另跑 `biopytools annovar`。

## 功能概述 | Overview

- 端到端 RNA-seq 变异检测全流程，产出到 VCF 为止（不含注释）
- 多样本联合 calling：每样本 gVCF → 联合基因分型 → 一个多样本 VCF
- 流程：HISAT2 比对 → 每样本 gVCF（HaplotypeCaller）→ CombineGVCFs → GenotypeGVCFs → VariantFiltration → bcftools PASS
- 支持原始 FASTQ（自动跑 fastp 质控）或已清洗 FASTQ（跳过质控）两种输入
- 自动构建共享基因组索引（faidx / dict / HISAT2 索引 / 可选剪接位点）
- 断点续传：索引、质控、比对、call、联合检测各步均可跳过已完成项

## 快速开始 | Quick Start

```bash
biopytools rnaseq2vcf -g genome.fa --gff3 anno.gff3 -i reads/ -o out/
```

最小输入：参考基因组 FASTA + 原始 FASTQ 目录（GFF3 可选）。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| 变异位点 variant | 基因组上某个位置，样本与参考序列「不一样」的地方（SNP 单碱基换、INDEL 插入缺失） |
| VCF | 记录变异位点的标准文本格式 |
| gVCF | 「每个位置都记录」的中间格式（含无变异位置），用于多样本联合检测 |
| 联合 calling | 把所有样本的信息放一起判断，比单样本更准 |
| 比对 alignment | 把读段放到参考基因组上找位置 |
| 质控 QC | 用 fastp 去掉低质量读段和接头 |
| PASS | 通过过滤、可信度高的变异标记 |
| 注释 annotation | 给变异位点加上「落在哪个基因、会不会改蛋白」等信息（本工具不包含，需另跑） |

## 输入 | Input

### 参考文件

- `-g/--genome`：参考基因组 FASTA（必需）
- `--gff3`：参考 GFF3（可选，提供则用于提取 HISAT2 剪接位点；不提供则 HISAT2 自行发现剪接位点）

### 测序数据（二选一）

1. **原始 FASTQ**（`-i/--input`）：目录内成对 FASTQ，会自动跑 fastp 质控
2. **已清洗 FASTQ**（`--clean-fastq-dir`）：已质控，跳过 QC 步骤

默认按 `_1.clean.fq.gz`/`_1.fq.gz`/`_1.fastq.gz` 等常见后缀自动识别成对样本（R1/R2），也可用 `--read1-pattern`/`--read2-pattern` 手动指定后缀。

## 参数说明 | Parameters

### 输入输出 | Input & output

**通俗理解|In plain words:** `-g` 是参考基因组；测序数据要么给 `-i`（原始，会先质控）要么给 `--clean-fastq-dir`（已清洗，跳过质控），两者必须给一个。`--gff3` 可选，有基因注释时提供能帮助比对更准。`-o` 是结果目录。

- `-g/--genome`：参考基因组 FASTA
- `--gff3`：参考 GFF3（可选）
- `-i/--input`：原始 FASTQ 目录（跑 fastp）
- `--clean-fastq-dir`：已清洗 FASTQ 目录（跳过 QC）
- `-o/--output-dir`：输出目录（默认当前目录 `.`）

### 变异检测与过滤 | Variant calling & filtering

**通俗理解|In plain words:** `--min-conf` 是 HaplotypeCaller 判定变异的最低置信度，调高更保守（更少但更可信），调低更敏感。后四项是 GATK 标准硬过滤阈值：`--fs-threshold`（链偏好，太高说明一条链的假信号）、`--qd-threshold`（质量/深度比，太低说明质量不可靠）、`--cluster-window`/`--cluster-size`（在多少 bp 窗口内聚了太多变异视为簇，常见于比对错误区域）。**这些都是 GATK 官方推荐默认值，一般不用动。**

- `--min-conf`：HaplotypeCaller 最小置信度（默认 20）
- `--fs-threshold`：FS 过滤阈值（默认 30.0）
- `--qd-threshold`：QD 过滤阈值（默认 2.0）
- `--cluster-window`：SNP cluster 过滤窗口 bp（默认 35）
- `--cluster-size`：SNP cluster 过滤数量（默认 3）

### 样本命名 | Sample naming

**通俗理解|In plain words:** 控制「怎么从文件名里识别 R1/R2 和样本名」。默认自动识别 `_1.clean.fq.gz`/`_1.fq.gz` 等常见后缀，**只有你的命名不标准时才需要手动指定**。

- `--read1-pattern` / `--read2-pattern`：R1/R2 后缀（默认自动识别）

### 运行与断点续传 | Runtime & checkpoint

**通俗理解|In plain words:** `-t` 是线程数。`-s/--step 0` 表示「只建索引就停」（想提前把共享索引建好时用）；省略则跑全流程。`--no-checkpoint` 关断点续传；`-f/--force` 忽略断点强制重跑；`--skip-qc` 跳过 fastp；`--dry-run` 只打印命令。**正常跑全流程这些都不用加。**

- `-t/--threads`：线程数（默认 12）
- `-s/--step`：0 = 仅建索引（省略 = 全流程）
- `--no-checkpoint`：关闭断点续传
- `-f/--force`：忽略断点重跑
- `--skip-qc`：跳过 fastp
- `--dry-run`：只打印命令
- `--log-file` / `--log-level`：日志文件 / 级别

## 分析流程 | Pipeline

```text
预检(工具存在性 + 磁盘空间 >= 200GB)
    |
    v
共享索引: samtools faidx + gatk CreateSequenceDictionary + hisat2-build (+ 可选剪接位点)
    |
    v
逐样本:
   +- fastp 质控(原始 FASTQ 时)
   +- HISAT2 --dta 比对 -> 排序 BAM
   +- GATK: 加读组 -> MarkDuplicates -> SplitNCigarReads -> HaplotypeCaller(gVCF)
    |
    v
联合 calling:
   +- CombineGVCFs(合并所有样本 gVCF)
   +- GenotypeGVCFs(联合基因分型)
   +- VariantFiltration(FS/QD/cluster 过滤)
   +- bcftools view -f PASS(提取通过项)
    |
    v
一个多样本 VCF + 分析报告
```

## 输出 | Output

```text
output_dir/
+-- genome_index/                      # 共享基因组索引(faidx/dict/hisat2)
+-- 01_qc/{sample}_1.clean.fq.gz      # fastp 质控后读段
+-- 02_align/
|   +-- {sample}.sorted.bam           # 比对结果
|   +-- {sample}.hisat2.log           # 比对日志
+-- 03_calling/                        # 每样本 calling 中间产物
|   +-- {sample}.rg.bam / .dedup.bam / .split.bam
|   +-- {sample}.g.vcf.gz (+ .tbi)    # 每样本 gVCF
|   +-- {sample}.metrics
+-- 04_joint/                          # 联合 calling 结果
|   +-- combined.g.vcf.gz             # 合并后 gVCF
|   +-- joint.vcf.gz                  # 联合基因分型 VCF
|   +-- joint.filtered.vcf.gz         # 过滤后(带 FILTER 标记)
|   +-- all_samples.pass.vcf.gz       # 最终 PASS VCF(主要产物)
|   +-- all_samples.pass.vcf.gz.tbi   # 索引
+-- 00_pipeline_info/checkpoints/      # 断点标记
+-- 99_logs/pipeline.log               # 日志
+-- tmp/                               # GATK 临时文件
+-- ANALYSIS_REPORT.txt                # 分析报告
```

## 结果解读 | Interpreting Results

### 1. 最终 VCF（`04_joint/all_samples.pass.vcf.gz`）

**通俗理解|In plain words:** 这是最终成果——一个包含全部样本、已经过滤的多样本 VCF。可直接送 ANNOVAR 注释或做下游群体/进化分析。

### 2. 分析报告（`ANALYSIS_REPORT.txt`）

**通俗理解|In plain words:** 一键看懂「过滤前后各有多少变异、每个样本有多少变异」。

- 过滤前总数（`joint.vcf.gz`）与过滤后 PASS 数（`all_samples.pass.vcf.gz`）对比，以及被过滤掉的数量与百分比
- 每样本 PASS 变异数（非参考基因型计数）
- 被滤变异可在 `joint.filtered.vcf.gz` 里复查（FILTER 列标有 FS/QD/SNPCluster 原因）

### 3. 过滤后仍保留的记录怎么读

VCF 的 FILTER 列为 `PASS` 的即通过项。被标记 FS（链偏好）、QD（低质量）、SNPCluster（成簇）的记录仍在 `joint.filtered.vcf.gz` 中，供人工复核。

## 参数选择建议 | Parameter Guidance

- **`--min-conf`**：标准 20；若想更保守减少假阳性可调到 30
- **`-i` vs `--clean-fastq-dir`**：数据未质控用 `-i`（自动 fastp）；已用 fastp 等质控过用 `--clean-fastq-dir` 省时间；`--skip-qc` 也可强制跳过质控
- **`--gff3`**：有参考基因注释建议提供，HISAT2 会用已知剪接位点，比对更准
- **`-s/--step 0`**：想先把共享索引建好（如多批样本复用同一参考）时用
- **`-f/--force`**：改参数后想彻底重跑时加，会绕过断点

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-g, --genome` | 必填 |  | 参考基因组 FASTA｜Reference genome FASTA |
| `--gff3` | — |  | 参考 GFF3(可选,HISAT2 剪接位点;不提供则 HISAT2 de novo 发现 junction)｜Reference GFF3 (optional, HISAT2 splice sites) |
| `-i, --input` | — |  | 原始 FASTQ 目录(跑 fastp)｜Raw FASTQ dir (runs fastp) |
| `--clean-fastq-dir` | — |  | 已清洗 FASTQ 目录(跳过 QC)｜Clean FASTQ dir (skip QC) |
| `-o, --output-dir` | `.` |  | 输出目录｜Output dir |
| `-t, --threads` | `12` | int | 线程数｜Threads |
| `--min-conf` | `20` | int | HaplotypeCaller 最小置信度｜min confidence |
| `--fs-threshold` | `30.0` | float | FS 过滤阈值(标记大于此值)｜FS filter threshold (mark if greater) |
| `--qd-threshold` | `2.0` | float | QD 过滤阈值(标记小于此值)｜QD filter threshold (mark if lower) |
| `--cluster-window` | `35` | int | SNP cluster 过滤窗口(bp)｜SNP cluster filter window (bp) |
| `--cluster-size` | `3` | int | SNP cluster 过滤数量｜SNP cluster filter count |
| `--read1-pattern` | — |  | R1 后缀(默认自动识别 _1.clean.fq.gz/_1.fq.gz 等)｜R1 suffix (auto-detected by default) |
| `--read2-pattern` | — |  | R2 后缀(默认自动识别)｜R2 suffix (auto-detected by default) |
| `-s, --step` | — | IntRange | 0=仅建索引｜index only;省略=全流程｜omit for full pipeline |
| `--no-checkpoint` | — |  | 关闭断点续传｜Disable checkpoint |
| `-f, --force` | — |  | 忽略断点重跑｜Force rerun |
| `--dry-run` | — |  | 只打印命令｜Dry run |
| `--skip-qc` | — |  | 跳过 fastp｜Skip QC |
| `--log-file` | — |  | 日志文件｜Log file |
| `--log-level` | `INFO` |  | 日志级别｜Log level |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-g, --genome` | 必填 |  | 参考基因组 FASTA｜Reference genome FASTA |
| `--gff3` | — |  | 参考 GFF3(可选,HISAT2 剪接位点)｜Reference GFF3 (optional, HISAT2 splice sites) |
| `-i, --input` | — |  | 原始 FASTQ 目录｜Raw FASTQ dir (runs fastp) |
| `--clean-fastq-dir` | — |  | 已清洗 FASTQ 目录(跳过 QC)｜Clean FASTQ dir (skip QC) |
| `-o, --output-dir` | `.` |  | 输出目录｜Output dir |
| `-t, --threads` | `12` | int | 线程数｜Threads (default 12) |
| `--min-conf` | `20` | int | HaplotypeCaller 最小置信度｜min confidence |
| `--fs-threshold` | `30.0` | float |  |
| `--qd-threshold` | `2.0` | float |  |
| `--cluster-window` | `35` | int |  |
| `--cluster-size` | `3` | int |  |
| `--read1-pattern` | — |  | R1 后缀(默认自动识别 _1.clean.fq.gz/_1.fq.gz 等)｜R1 suffix (auto) |
| `--read2-pattern` | — |  | R2 后缀(默认自动识别)｜R2 suffix (auto) |
| `-s, --step` | — | 0/1/2/3/4 | 0=仅建索引｜index only;省略=全流程｜omit for full pipeline |
| `--no-checkpoint` | — | store_true |  |
| `-f, --force` | — | store_true |  |
| `--dry-run` | — | store_true |  |
| `--skip-qc` | — | store_true |  |
| `--log-file` | — |  |  |
| `--log-level` | `INFO` |  |  |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- HISAT2（`hisat2`、`hisat2-build`、`hisat2_extract_splice_sites.py`）
- fastp（通过 `biopytools fastp` 调用）
- GATK（HaplotypeCaller / CombineGVCFs / GenotypeGVCFs / VariantFiltration 等）
- SAMtools（`samtools`）、bcftools
- 工具路径默认按 conda 功能域环境解析（hisat2/fastp 在 rna 域、samtools/bcftools/gatk 在 align 域），也可用环境变量覆盖

## 常见问题 | FAQ

**Q1：支持断点续传吗？**
支持。共享索引用 `00_pipeline_info/checkpoints/*.done` 标记；每样本质控/比对/call 按输出文件存在性跳过；联合 calling 以 `all_samples.pass.vcf.gz` 存在为完成标志。`--no-checkpoint` 关闭，`-f/--force` 强制重跑。

**Q2：为什么提示磁盘空间不足？**
预检要求输出目录所在盘有至少 200 GB 可用空间（GATK 中间文件很大）。空间不足会在日志中警告；请清理或换盘。

**Q3：断点续传时怎么判断 gVCF 是否完整？**
只有 `{sample}.g.vcf.gz` 与其索引 `.tbi` **同时存在**才视为已完成。若中断导致只有 `.g.vcf.gz` 缺 `.tbi`，程序会视为残缺并重跑该样本，避免用到损坏结果。

**Q4：输出 VCF 里有注释吗？**
没有。本工具只做到 VCF（含过滤），ANNOVAR 注释需另跑 `biopytools annovar`。

**Q5：`-s/--step` 能指定别的步骤吗？**
CLI 入口只允许 `-s 0`（仅建索引）或省略（全流程）。没有「只跑到第 N 步」的用法；省略即跑完整流程。

**Q6：样本命名不标准怎么办？**
用 `--read1-pattern`/`--read2-pattern` 手动指定 R1/R2 后缀。默认已支持 `_1.clean.fq.gz`/`_1.fq.gz`/`_1.fastq.gz` 等常见命名。
