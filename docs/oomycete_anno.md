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

## 依赖 | Dependencies

- **Augustus / BRAKER**: 基因预测
- **GeneMark-ES (gmes_petap)**: 自训练基因预测
- **minimap2 / samtools**: RNA-seq 与三代比对
- **miniprot**: Phase3 效应子比对
- **RepeatMasker / EDTA**: 重复序列(按阶段)

## 相关链接 | References

- [项目主页](https://github.com/lixiang117423/biopytools)
