# Parabricks GPU 加速 WGS 分析 | GPU-accelerated WGS Analysis

一句话理解：**用 NVIDIA GPU(通过 Parabricks 容器)把一批 FASTQ 快速跑完「比对→变异检测→多样本联合」，产出 BAM、个体 gVCF 和联合 VCF**，解决「WGS 分析太慢、要用 GPU 加速」的问题。

## 功能概述 | Overview { #overview }

- 输入一个 FASTQ 目录(按 R1/R2 模式配对)，通过 apptainer/singularity 容器运行 Parabricks 的 `pbrun` 命令
- 工作流可选：`fq2bam`(比对)、`haplotypecaller`(变异检测)、`genotypegvcf`(联合)、`all`(完整流程，默认)
- 默认输出 gVCF 格式并做 Joint Calling，产出个体 gVCF 和联合 VCF
- 断点续传：已完成的样本跳过，联合结果已存在也跳过
- 生成分析总结报告 `parabricks_analysis_summary.txt`

## 快速开始 | Quick Start { #quick-start }

```bash
biopytools parabricks -i /data/fastq -o results/ -r /ref/genome.fa
```

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语<br>Term | 通俗理解<br>In plain words |
|------|------|
| GPU 加速 | 用显卡并行计算替代 CPU，把生信流程大幅提速 |
| Parabricks | NVIDIA 出的 GPU 生信套件，核心命令叫 `pbrun` |
| 容器(.sif) | 把软件和依赖打包成单个文件，用 apptainer/singularity 运行 |
| fq2bam | 把 FASTQ 比对成 BAM 的步骤 |
| haplotypecaller | 从 BAM 检测变异、产出 VCF/gVCF 的步骤 |
| genotypegvcf | 把多个 gVCF 合并成一个联合 VCF 的步骤(联合基因分型) |
| gVCF | 每个位点都有记录(含未变异位点)的 VCF，适合多样本联合 |
| joint calling | 多样本联合变异检测，提高低频变异检出率 |

## 输入 | Input { #input }

### FASTQ 目录

`-i` 指向一个目录，默认按 `*_1.clean.fq.gz` / `*_2.clean.fq.gz` 配对(可用 `--read1-pattern` / `--read2-pattern` 改)。

```text
fastq/
├── sample1_1.clean.fq.gz
├── sample1_2.clean.fq.gz
├── sample2_1.clean.fq.gz
└── sample2_2.clean.fq.gz
```

### 参考基因组

`-r` 指定参考基因组 FASTA。

## 参数说明 | Parameters { #parameters }

### 必需参数 | Required

**通俗理解|In plain words:** `-i` 是 FASTQ 目录，`-o` 是结果目录，`-r` 是参考基因组。三个都要给。

### 工作流参数 | Workflow

**通俗理解|In plain words:** `-w` 决定跑到哪一步。`all`(默认)从头跑到联合；`fq2bam` 只比对；`haplotypecaller` 只做变异检测；`genotypegvcf` 只做联合(此时不需要 FASTQ，会自动检测输出目录里已有的 gVCF)。**想补跑某一步时很有用，默认 `all` 最省心。**

### 运行参数 | Runtime

**通俗理解|In plain words:** `--parabricks-path` 是容器文件位置(默认 `~/software/containers/parabricks.sif`)；`--tmp-dir` 临时目录(默认输出目录下 `tmp/`)；`-t` 线程数(默认 12)为预留参数，当前未传入 `pbrun` 命令。**一般不用动。**

### GVCF 与 Joint Calling 参数 | GVCF & Joint calling

**通俗理解|In plain words:** `--gvcf/--no-gvcf`(默认输出 gVCF)决定个体结果是 gVCF 还是普通 VCF；`--joint-calling/--no-joint-calling`(默认开)决定是否做联合；`--combined-output`(默认 `combined.g.vcf`)是联合结果文件名。**联合需要至少 2 个样本，且只输出 gVCF 时才有意义，一般保持默认。**

### 质控参数 | Quality control

**通俗理解|In plain words:** `--min-confidence`(默认 30)、`--min-base-quality`(默认 20)、`--ploidy`(默认 2)、`--pcr-indel-model`(默认 CONSERVATIVE)这几个参数当前被接受并记录到配置，但暂未直接传入 `pbrun` 命令。**一般不用动。**

### 文件模式参数 | File patterns

**通俗理解|In plain words:** 告诉程序 R1/R2 怎么命名，只有和默认 `*_1.clean.fq.gz` / `*_2.clean.fq.gz` 不同时才需要改，且必须含 `*`。

## 分析流程 | Pipeline { #pipeline }

**通俗理解|In plain words:** 先检查容器和依赖，再逐个样本「比对→变异检测」，最后把所有 gVCF 合并成联合 VCF，收尾写总结。

```text
输入 FASTQ 目录 + 参考基因组
    |
    ▼
步骤0: 依赖检查(parabricks.sif + apptainer/singularity + bcftools)
    |
    ▼
步骤1: 初始化容器(apptainer/singularity exec --nv)
    |
    ▼
步骤2: 每样本 pbrun fq2bam(→ bam/<样本>.sorted.bam)
    |
    ▼
步骤3: 每样本 pbrun haplotypecaller(→ vcf/<样本>.g.vcf.gz)
    |
    ▼
步骤4: Joint Calling(pbrun genotypegvcf → vcf/combined.g.vcf.gz + bcftools index)
    |
    ▼
步骤5: 生成总结报告 parabricks_analysis_summary.txt
```

## 输出 | Output { #output }

```text
results/
├── bam/                             # 每样本比对结果
│   └── sample1.sorted.bam
├── vcf/                             # 变异结果
│   ├── sample1.g.vcf.gz             # 个体 gVCF(默认)
│   └── combined.g.vcf.gz            # 联合 VCF(默认 joint calling 时)
├── tmp/                             # 临时文件目录
├── parabricks_analysis_summary.txt  # 分析总结报告
└── parabricks_analysis.log          # 运行日志
```

## 结果解读 | Interpreting Results { #interpreting }

### 1. vcf/<样本>.g.vcf.gz(个体 gVCF)

**通俗理解|In plain words:** 每个样本自己的变异记录(gVCF 格式，含未变异位点)。用 `--no-gvcf` 时改为普通 `<样本>.vcf.gz`。

### 2. vcf/combined.g.vcf.gz(联合 VCF，核心)

**通俗理解|In plain words:** 把所有样本的 gVCF 合并后重新判读的群体变异结果，是多样本项目最常用的下游输入。

- 文件名由 `--combined-output` 决定(默认 `combined.g.vcf`，实际写成 `combined.g.vcf.gz`)
- 联合需要至少 2 个 gVCF；合并后会自动用 bcftools index 建索引

### 3. parabricks_analysis_summary.txt(总结报告)

**通俗理解|In plain words:** 一份收尾清单，含分析配置、质控参数、文件统计和目录结构。

## 参数选择建议 | Parameter Guidance { #guidance }

- **正常多样本流程**：默认 `-w all` + 默认 gVCF/joint calling，一次拿全
- **只想补做联合**：`-w genotypegvcf`(自动检测已有 gVCF，无需 FASTQ)
- **只想重新比对**：`-w fq2bam`
- **单样本、不需要联合**：`--no-joint-calling`
- **想要普通 VCF 而非 gVCF**：`--no-gvcf`

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input-dir, -i` | 必填 |  | 输入目录(FASTQ文件)｜Input directory containing FASTQ files |
| `--output-dir, -o` | 必填 | Path | 输出目录｜Output directory |
| `--reference, -r` | 必填 |  | 参考基因组文件｜Reference genome file |
| `--workflow, -w` | `all` | fq2bam/haplotypecaller/genotypegvcf/all | 工作流程｜Workflow: fq2bam, haplotypecaller, genotypegvcf, or all |
| `--threads, -t` | `12` | int | 线程数｜Number of threads |
| `--parabricks-path` | `~/software/containers/parabricks.sif` | str | Parabricks程序路径｜Parabricks program path |
| `--tmp-dir` | — | Path | 临时目录｜Temporary directory |
| `--gvcf/--no-gvcf` | `True` |  | 输出GVCF格式｜Output GVCF format |
| `--joint-calling/--no-joint-calling` | `True` |  | 启用Joint Calling｜Enable Joint Calling |
| `--combined-output` | `combined.g.vcf` |  | Joint Calling输出文件名｜Joint Calling output filename |
| `--min-confidence` | `30` | int | 最小置信度阈值｜Minimum confidence threshold |
| `--min-base-quality` | `20` | int | 最小碱基质量阈值｜Minimum base quality threshold |
| `--ploidy` | `2` | int | 倍性｜Ploidy |
| `--pcr-indel-model` | `CONSERVATIVE` | str | PCR indel模型｜PCR indel model |
| `--read1-pattern` | `*_1.clean.fq.gz` | str | R1文件模式｜R1 file pattern |
| `--read2-pattern` | `*_2.clean.fq.gz` | str | R2文件模式｜R2 file pattern |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input-dir` | 必填 |  | 输入目录 (FASTQ文件)｜Input directory (FASTQ files) |
| `-o, --output-dir` | 必填 |  | 输出目录｜Output directory |
| `-r, --reference` | 必填 |  | 参考基因组｜Reference genome |
| `-w, --workflow` | `all` | fq2bam/haplotypecaller/genotypegvcf/all | 工作流程｜Workflow: fq2bam(比对), haplotypecaller(变异检测), genotypegvcf(合并), all(完整流程) |
| `-t, --threads` | `88` | int | 线程数｜Threads |
| `--parabricks-path` | `~/software/containers/parabricks.sif` |  | parabricks路径｜parabricks path |
| `--tmp-dir` | — |  | 临时目录｜Temporary directory |
| `--min-confidence` | `30` | int | 最小置信度｜Min confidence |
| `--min-base-quality` | `20` | int | 最小碱基质量｜Min base quality |
| `--ploidy` | `2` | int | 倍性｜Ploidy |
| `--pcr-indel-model` | `CONSERVATIVE` |  | PCR indel模型｜PCR indel model |
| `--read1-pattern` | `*_1.clean.fq.gz` |  | R1文件模式｜R1 pattern |
| `--read2-pattern` | `*_2.clean.fq.gz` |  | R2文件模式｜R2 pattern |
| `--no-gvcf` | — | store_false | 输出VCF而非GVCF｜Output VCF instead of GVCF |
| `--no-joint-calling` | `True` | store_false | 禁用Joint Calling｜Disable Joint Calling |
| `--combined-output` | `combined.g.vcf` |  | Joint Calling输出文件名｜Joint Calling output filename |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies { #dependencies }

- NVIDIA Parabricks 容器(默认 `~/software/containers/parabricks.sif`)
- apptainer 或 singularity(容器运行时，需能跑 GPU 的 `--nv` 模式)
- bcftools(启用 joint calling 时，为联合 VCF 建索引)
- NVIDIA GPU(容器以 `--nv` 模式运行)

## 常见问题 | FAQ { #faq }

**Q1：报「未找到 Apptainer 或 Singularity」？**
Parabricks 通过容器运行，系统里必须能 `which apptainer` 或 `which singularity` 找到其中一个，且需要 GPU(容器用 `--nv` 模式)。

**Q2：支持断点续传吗？**
支持。单个样本按 `bam/<样本>.sorted.bam` 和 `vcf/<样本>.g.vcf.gz` 是否都存在判断跳过；联合结果 `vcf/combined.g.vcf.gz` 已存在也跳过。换参数重跑前需先删除旧产物。

**Q3：`genotypegvcf` 模式为什么不用给 FASTQ？**
这个模式只做联合，会自动扫描输出目录 `vcf/` 下已有的 `*.g.vcf.gz` 文件来合并，所以只需要参考基因组和已有 gVCF。

**Q4：`-t` 线程数为什么没效果？**
当前版本线程数为预留参数，未直接传入 `pbrun` 命令(Parabricks 主要靠 GPU 加速)。

**Q5：质控参数(如 `--min-confidence`)真的生效吗？**
当前版本这几个参数被接受并记录到配置/总结报告，但暂未直接传入 `pbrun` 命令，以源码当前实现为准。
