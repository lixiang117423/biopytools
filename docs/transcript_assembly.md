# 转录本组装分析模块 | Transcript Assembly Module

**支持 FASTQ 或 BAM 直入的转录本组装工具,输出 GFF3 基因结构 | Transcript assembly from FASTQ or BAM, outputs GFF3 gene structure**

## 功能概述 | Overview

转录本组装模块基于 HISAT2 比对和 StringTie 组装,支持两种输入路径:

- **FASTQ 路径**(二代):HISAT2 建索引 → 比对(--dta) → BAM 排序 → StringTie 组装
- **BAM 直入路径**(二代短读 或 三代长读):跳过比对,直接 StringTie 组装

两条路径在 StringTie 逐样本组装处汇合,逐 BAM 组装后 `--merge` 合并,最终输出 **GFF3 基因结构**(`merged.gff3`)。支持长读(`stringtie -L`)、可选 reference-guided(`-G` 参考注释)、多样本合并与断点续传。

## 主要特性 | Key Features

- **双输入模式**:FASTQ 目录(二代)或 BAM 文件直入(短读/长读),`-i` 与 `-b` 互斥
- **长读支持**:自动检测读长,三代 BAM(PacBio/ONT)用 `stringtie -L`
- **GFF3 主输出**:`gffread` 将合并 GTF 转为 GFF3;可选 `--transcripts` 额外输出 cDNA
- **可选 guided**:提供 `--guide-gff` 参考注释时,组装与合并均启用 `-G`(reference-guided)
- **多样本合并**:多个 BAM 逐个组装后 `stringtie --merge`;单 BAM 跳过合并
- **断点续传**:已完成的步骤自动跳过
- **灵活步骤控制**:6 个步骤独立运行或完整流程

## 快速开始 | Quick Start

### BAM 直入(短读或长读)→ GFF3 | BAM input → GFF3

```bash
# 单个 BAM(自动检测短/长读)
biopytools transcript-assembly -b sample.sorted.bam -o ./out

# 多个 BAM(逐个组装后合并);多个 -b 重复给出
biopytools transcript-assembly -b s1.bam -b s2.bam -o ./out
```

### FASTQ(二代)→ GFF3 | FASTQ input → GFF3

```bash
biopytools transcript-assembly -g genome.fasta -i ./clean_data -o ./out
```

### reference-guided(有参考注释)| With reference annotation

```bash
biopytools transcript-assembly -b sample.bam --guide-gff ref.gff3 -o ./out
```

## 参数说明 | Parameters

### 输入(二选一)| Input (mutually exclusive)

| 参数 | 描述 | 示例 |
|------|------|------|
| `-b, --bam` | 输入 BAM 文件(一个或多个,可多次 `-b`,与 `-i` 互斥) | `-b s1.bam -b s2.bam` |
| `-i, --input` | 输入 FASTQ 文件目录(与 `-b` 互斥) | `-i ./clean_data` |
| `-g, --genome` | 基因组 FASTA(FASTQ 模式或 `--transcripts` 时必需;纯 BAM→GFF3 不需要) | `-g genome.fasta` |
| `--guide-gff` | 参考注释 GTF/GFF3(启用 `-G` guided 组装与合并) | `--guide-gff ref.gff3` |
| `--read-type` | 读长类型:`auto`(默认,自动检测)/`short`/`long` | `--read-type long` |
| `--transcripts` | 额外输出 `transcripts.fa` cDNA(需 `-g`) | `--transcripts` |

### 必需 | Required

| 参数 | 描述 |
|------|------|
| `-o, --output` | 输出目录 |

### 可选 | Optional

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-p, --pattern` | `*_1.clean.fq.gz` | FASTQ 文件命名模式(`*` 为样本名占位符) |
| `-t, --threads` | `12` | 线程数 |
| `--sample-timeout` | `43200` | 单样本处理超时(秒),默认 12 小时 |

### 步骤控制 | Step Control

| 参数 | 描述 | BAM 模式 |
|------|------|----------|
| `-s 1` | 构建 HISAT2 基因组索引 | 不可用 |
| `-s 2` | HISAT2 比对(--dta) | 不可用 |
| `-s 3` | SAM 转排序 BAM | 不可用 |
| `-s 4` | StringTie 逐样本组装(长读加 `-L`,guided 加 `-G`) | 可用 |
| `-s 5` | StringTie 合并 GTF(单 BAM 跳过) | 可用 |
| `-s 6` | gffread 输出 GFF3(+可选 cDNA) | 可用 |

> BAM 模式仅支持 step 4-6;step 1-3 仅 FASTQ 模式可用。

### 高级选项 | Advanced Options

| 参数 | 描述 |
|------|------|
| `--force` | 强制重新处理已完成的步骤 |
| `-v` | 增加输出详细程度(可多次使用) |
| `--quiet` | 静默模式,仅输出错误信息 |

## 读长自动检测 | Read-type Auto-detection

BAM 模式下(`--read-type auto`,默认),每个 BAM 采样前 1000 条 read 的长度中位数:

- 中位数 < 500 bp → `short`(Illumina,~150bp)→ `stringtie` 正常模式
- 中位数 ≥ 500 bp → `long`(PacBio HiFi ~10-15kb / ONT ~1-10kb)→ `stringtie -L`

可用 `--read-type short|long` 强制指定,跳过检测。

## 输出结果 | Output Results

### 输出目录结构 | Output Directory Structure

```
out/
├── 00_pipeline_info/          # 流程元数据
│   └── software_versions.yml  # 软件版本 + 参数
├── 01_hisat2_index/           # 仅 FASTQ 路径
├── 02_hisat2_align/           # 仅 FASTQ 路径
├── 03_bam_sort/               # 仅 FASTQ 路径(BAM 直入原位引用,不复制)
├── 04_stringtie/              # 逐样本 GTF
│   ├── sample1.gtf
│   └── gtf_list.txt           # GTF 列表(多样本时)
├── 05_merge/                  # 合并 GTF(>1 样本时)
│   └── merged.gtf
├── 06_gff3/                   # 主输出
│   └── merged.gff3            # GFF3 基因结构
├── 06_transcripts/            # 仅 --transcripts
│   └── transcripts.fa         # cDNA 序列
└── 99_logs/
    └── pipeline.log
```

### 关键输出文件 | Key Output Files

| 文件 | 描述 |
|------|------|
| `06_gff3/merged.gff3` | **主输出**:GFF3 基因结构 |
| `05_merge/merged.gtf` | 合并后的 GTF(>1 样本时) |
| `06_transcripts/transcripts.fa` | 可选:cDNA 序列(需 `--transcripts -g`) |
| `00_pipeline_info/software_versions.yml` | 软件版本与运行参数 |

## 流程说明 | Pipeline Details

- **Step 1-3**(仅 FASTQ):HISAT2 建索引 → `--dta` 比对 → `samtools sort/index`。
- **Step 4**:逐样本 StringTie 组装。长读加 `-L`;提供 `--guide-gff` 时加 `-G`。
- **Step 5**:`stringtie --merge` 合并所有样本 GTF。`-G` 仅接参考注释(**不接基因组 FASTA**,已修复原 bug);单 BAM 跳过此步。
- **Step 6**:`gffread merged.gtf -o merged.gff3 -F` 输出 GFF3(主);`--transcripts` 时额外 `gffread -g genome -w transcripts.fa`。

## 行为变化(相对旧版本)| Behavior Changes

1. **主输出从 `transcripts.fa` 改为 `merged.gff3`**。需要 cDNA 加 `--transcripts`。
2. **`-s 6` 现产出 GFF3**(不再是 transcripts.fa)。
3. **`-g` 从必需变可选**(纯 BAM→GFF3 不需要 genome)。
4. **`-i` 从必需变可选**(可与 `-b` 二选一)。
5. **merge 不再接收 genome**(`-G` 仅接注释)——修复原 `--merge -G genome.fa` 误用。

## 依赖软件 | Dependencies

均在 conda 环境 `RNA_Seq` 中(路径可通过环境变量 `STRINGTIE_PATH`/`GFFREAD_PATH`/`HISAT2_PATH`/`HISAT2_BUILD_PATH` 覆盖;samtools 走 `SAMTOOLS_PATH`):

- **StringTie** - 转录本组装与合并(长读 `-L`)
- **gffread** - GTF→GFF3 转换与 cDNA 提取
- **HISAT2** - 二代比对(仅 FASTQ 路径)
- **SAMtools** - BAM 校验/索引/排序

## 常见问题 | Troubleshooting

**Q: BAM 模式报"BAM 非坐标排序"**

BAM 必须是 coordinate-sorted。先 `samtools sort -o sorted.bam input.bam` 再传入。模块会自动建索引(若 `.bai` 缺失)。

**Q: 如何同时跑短读+长读**

分别用 `-b` 给出多个 BAM(如 `-b short.bam -b long.bam`),模块逐个组装(各自用检测到的模式)再 `--merge`。注意:这是"逐个组装+合并",不是单跑 `--mix`。

**Q: 如何从中断处继续**

直接重跑相同命令,已完成步骤自动跳过。`--force` 强制重跑。

**Q: 找不到样本 FASTQ 文件对**

检查命名模式,用 `-p` 指定,如 `-p "*_R1.fq.gz"`。
