# 疫霉菌基因组注释 | Oomycete Genome Annotation

**疫霉菌基因组的T2T Augustus注释流程, 分阶段整合RNA-seq/蛋白/三代/LTR/效应子证据 | T2T Augustus annotation pipeline for oomycetes, phased integration of RNA-seq/protein/long-read/LTR/effector evidence**

## 功能概述 | Overview

oomycete_anno 模块针对疫霉菌(*Phytophthora* 等)基因组, 跑 T2T Augustus 注释流程, 分阶段整合多源证据:

- **Phase1 主注释**: GeneMark-ES + Augustus, 整合 RNA-seq 比对 hints
- **Phase2 证据增强**: 蛋白 hints / LTR 转座子 / 三代转录本(Iso-seq)证据
- **Phase3 效应子位点救援**: 用已知效应子全长蛋白 miniprot 比对, 替换/补回 Augustus 在效应子簇位点的错注(嵌合/截断)与漏注

针对疫霉效应子(RxLR/CRN)多落在 TE 区、易被注释工具误伤的特点做了专门处理。

## 快速开始 | Quick Start

```bash
# 基础用法(RNA-seq 证据)
biopytools oomycete-anno -g genome.fa -s phytophthora --rnaseq-dirs rna1/ rna2/ -o out/ -t 24

# 加蛋白 + 三代 + 效应子(启用 Phase2/3)
biopytools oomycete-anno -g genome.fa -s psojae --rnaseq-dirs rna/ \
  --prot-seq proteins.faa --isoseq iso.fq --effectors effectors.faa -o out/
```

## 参数说明 | Parameters

### 必需参数 | Required

| 参数 | 描述 |
|------|------|
| `-g, --genome` | 基因组 FASTA |
| `-s, --species` | 物种名(Augustus/GeneMark 用) |

### 常用可选参数 | Common Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-o, --output-dir` | `./oomycete_anno_output` | 输出目录 |
| `-t, --threads` | `12` | 线程数 |
| `--rnaseq-dirs` | — | RNA-seq 目录(可多个) |
| `--prot-seq` | — | 蛋白证据(Phase2 hints) |
| `--isoseq` | — | 三代转录本文件/目录(Phase2) |
| `--effectors` | — | 已知效应子蛋白(Phase3 救援) |
| `--read1-pattern` / `--read2-pattern` | `_1/_2.clean.fq.gz` | R1/R2 文件后缀 |
| `--rna-strandness` | — | RNA-seq 链特异性 |
| `--no-soft-masking` | `False` | 不做 soft-mask |
| `--rescue-min-identity` | `0.85` | Phase3 救援最低 identity |
| `--rescue-conflict-overlap` | `0.50` | Phase3 冲突重叠阈值 |
| `--gmes-petap-path` | 自动 | GeneMark-ES 路径 |
| `--skip-repeat` / `--skip-rna` / `--skip-iso` / `--skip-protein` / `--skip-ltr` / `--skip-rescue` | `False` | 各阶段跳过开关 |

(运行 `biopytools oomycete-anno -h` 查看完整参数列表)

## 输出 | Output

- Augustus/GeneMark 注释结果(GFF3/蛋白)
- 各阶段中间产物(RepeatMasker、RNA-seq BAM、hints 等)
- Phase3 效应子救援结果
- 流程元数据与日志

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-g, --genome` | 必填 |  | 基因组 FASTA｜Genome FASTA |
| `-s, --species` | 必填 |  | Augustus 物种名(简单字母, 如 phytophthora)｜Augustus species name |
| `--rnaseq-dirs` | — |  | 二代 RNA-seq 目录(可多个)｜Short RNA-seq dir(s) |
| `--prot-seq` | — |  | 同源蛋白文件(Phase2)｜Homologous proteins (P2) |
| `--isoseq` | — |  | 三代转录本文件(Phase2)｜Long-read transcripts (P2) |
| `--effectors` | — |  | 已知效应子蛋白(Phase3 救援)｜Known effectors (P3 rescue) |
| `--read1-pattern` | `_1.clean.fq.gz` |  | R1 文件后缀模式｜R1 suffix pattern |
| `--read2-pattern` | `_2.clean.fq.gz` |  | R2 文件后缀模式｜R2 suffix pattern |
| `--rna-strandness` | `` |  | 链特异性: ''(非链特异性)/FR/RF｜Strandness: ''/FR/RF |
| `-o, --output-dir` | `./oomycete_anno_output` |  | 输出目录｜Output directory |
| `-t, --threads` | `12` |  | 线程数｜Threads |
| `--no-soft-masking` | `False` |  | 禁用软屏蔽(改用硬屏蔽)｜Disable soft masking |
| `--gmes-petap-path` | — |  | GeneMark gmes_petap.pl 路径｜GeneMark gmes_petap.pl path |
| `--genemark-perl-env` | — |  | GeneMark perl 提供环境｜GeneMark perl provider env |
| `--skip-repeat` | `False` |  | 跳过重复屏蔽｜Skip repeat masking |
| `--skip-rna` | `False` |  | 跳过 RNA-seq 比对｜Skip RNA-seq alignment |
| `--skip-iso` | `False` |  | 跳过三代转录本(Phase2)｜Skip long-read (P2) |
| `--skip-protein` | `False` |  | 跳过蛋白证据(Phase2)｜Skip protein hints (P2) |
| `--skip-ltr` | `False` |  | 跳过 LTR 注解(Phase2)｜Skip LTR annotation (P2) |
| `--skip-rescue` | `False` |  | 跳过效应子救援(Phase3)｜Skip effector rescue (P3) |
| `--rescue-min-identity` | `0.85` | float | 效应子救援 miniprot 最低 identity｜Rescue min identity |
| `--rescue-conflict-overlap` | `0.5` | float | Augustus 与效应子模型重叠>此比例则替换｜Conflict overlap fraction |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-g, --genome` | 必填 |  | 基因组 FASTA｜Genome FASTA |
| `-s, --species` | 必填 |  | Augustus 物种名(简单字母, 如 phytophthora)｜Augustus species name (simple, e.g. phytophthora) |
| `--rnaseq-dirs` | — |  | 二代 RNA-seq 目录(可多个)｜Short RNA-seq dir(s) |
| `--prot-seq` | — |  | 同源蛋白文件(Phase2)｜Homologous proteins (P2) |
| `--isoseq` | — |  | 三代转录本文件(Phase2)｜Long-read transcripts (P2) |
| `--effectors` | — |  | 已知效应子蛋白(Phase3 救援, 直接当基因模型替换错注/漏注位点)｜Known effectors (P3 rescue, used as gene models to fix loci) |
| `--read1-pattern` | `_1.clean.fq.gz` |  | R1 文件后缀模式｜R1 suffix pattern |
| `--read2-pattern` | `_2.clean.fq.gz` |  | R2 文件后缀模式｜R2 suffix pattern |
| `--rna-strandness` | `` |  | 链特异性: ''(非链特异性,默认) / FR / RF｜Strandness: ''(unstranded,default)/FR/RF |
| `-o, --output-dir` | `./oomycete_anno_output` |  | 输出目录｜Output directory |
| `-t, --threads` | `12` | int | 线程数｜Threads |
| `--no-soft-masking` | — | store_false | 禁用软屏蔽(改用硬屏蔽)｜Disable soft masking (use hard masking) |
| `--gmes-petap-path` | — |  | GeneMark gmes_petap.pl 路径(默认 ~/software/GeneMark/...)｜GeneMark gmes_petap.pl path |
| `--genemark-perl-env` | — |  | GeneMark perl 提供环境(默认 braker_v.3.0.8)｜GeneMark perl provider env |
| `--skip-repeat` | — | store_true | 跳过重复屏蔽｜Skip repeat masking |
| `--skip-rna` | — | store_true | 跳过 RNA-seq 比对｜Skip RNA-seq alignment |
| `--skip-iso` | — | store_true | 跳过三代转录本(Phase2)｜Skip long-read (P2) |
| `--skip-protein` | — | store_true | 跳过蛋白证据(Phase2)｜Skip protein hints (P2) |
| `--skip-ltr` | — | store_true | 跳过 LTR 注解(Phase2)｜Skip LTR annotation (P2) |
| `--skip-rescue` | — | store_true | 跳过效应子救援(Phase3)｜Skip effector rescue (P3) |
| `--rescue-min-identity` | `0.85` | float | 效应子救援 miniprot 最低 identity(全长靠 Target起=1+stop 判)｜Rescue min identity |
| `--rescue-conflict-overlap` | `0.5` | float | Augustus 基因与效应子模型重叠>此比例则替换｜Conflict overlap fraction |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- **Augustus / BRAKER**: 基因预测
- **GeneMark-ES (gmes_petap)**: 自训练基因预测
- **minimap2 / samtools**: RNA-seq 与三代比对
- **miniprot**: Phase3 效应子比对
- **RepeatMasker / EDTA**: 重复序列(按阶段)

## 相关链接 | References

- [项目主页](https://github.com/lixiang117423/biopytools)
