# gatk-joint | GATK Joint Genotyping 流程

**对一批 GVCF 文件执行 GATK 联合基因分型（Joint Genotyping）并完成 SNP/INDEL 硬过滤 | GATK Joint Genotyping pipeline with built-in SNP/INDEL hard filtering**

## 功能概述 | Overview

本模块封装 GATK4 的 Best Practices Joint Genotyping 流程，将一组样本的 GVCF 合并为群体级别的 VCF。标准流程包含：检测输入文件类型、构建 sample map、GenomicsDBImport、GenotypeGVCFs 联合分型、SelectVariants 分离 SNP/INDEL、VariantFiltration 按 GATK 推荐阈值硬过滤，最后用 BCFtools 合并。

适合人群重测序、群体遗传学分析等需要对多样本同时分型的场景。默认 SNP/INDEL 过滤阈值遵循 GATK 硬过滤推荐值，用户可通过参数微调。

## 快速开始 | Quick Start

```bash
# 输入目录包含多个样本的 .g.vcf.gz 文件
biopytools gatk-joint -i ./gvcf_folder/ -g reference.fa -o ./joint_output/
```

## 参数说明 | Parameters

### 必需参数 | Required

| 参数 | 描述 |
|------|------|
| `-i, --input` | 输入目录，包含 GVCF/VCF 文件 |
| `-g, --genome` | 参考基因组文件（`.fasta` / `.fa`） |

### 常用可选参数 | Common Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-o, --output-dir` | `./joint_genotyping_output` | 输出目录 |
| `-t, --threads` | `12` | 线程数 |
| `-m, --memory` | `100g` | Java 堆内存设置 |
| `-L, --intervals` | - | 分析区间（染色体名或区间文件） |
| `--gatk-path` | `gatk` | GATK 可执行文件路径 |
| `--bcftools-path` | `bcftools` | BCFtools 可执行文件路径 |

### SNP 过滤阈值 | SNP Thresholds

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--snp-qd` | `2.0` | QD 阈值 |
| `--snp-fs` | `60.0` | FS 阈值 |
| `--snp-mq` | `40.0` | MQ 阈值 |
| `--snp-mqrs` | `-12.5` | MQRankSum 阈值 |
| `--snp-rprs` | `-8.0` | ReadPosRankSum 阈值 |
| `--snp-sor` | `3.0` | SOR 阈值 |

### INDEL 过滤阈值 | INDEL Thresholds

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--indel-qd` | `2.0` | QD 阈值 |
| `--indel-fs` | `200.0` | FS 阈值 |
| `--indel-rprs` | `-20.0` | ReadPosRankSum 阈值 |
| `--indel-sor` | `10.0` | SOR 阈值 |

（运行 `biopytools gatk-joint -h` 查看完整参数列表）

## 输出 | Output

```
joint_genotyping_output/
├── raw_variants.vcf.gz             # GenotypeGVCFs 原始联合分型结果
├── raw_variants.vcf.gz.snp.vcf.gz  # 分离出的 SNP
├── raw_variants.vcf.gz.indel.vcf.gz  # 分离出的 INDEL
├── filtered.snp.vcf.gz             # 过滤后 SNP
├── filtered.indel.vcf.gz           # 过滤后 INDEL
├── joint_genotyped.vcf.gz          # 合并后的最终 VCF
└── *.log                           # 运行日志
```

## 依赖 | Dependencies

- [GATK4](https://gatk.broadinstitute.org/)（必需，提供 GenomicsDBImport、GenotypeGVCFs、SelectVariants、VariantFiltration）
- [BCFtools](http://www.htslib.org/)（用于 SNP/INDEL 合并）
- Java 8+（GATK 运行时依赖）
- 参考基因组需已建立 `.dict` 与 `.fai` 索引，GVCF 需已建立 `.tbi` 索引

## 引用 | Citation

- Poplin R et al. Scaling accurate genetic variant discovery to tens of thousands of samples. *bioRxiv*, 2018. doi:10.1101/201178
- Van der Auwera GA & O'Connor BD. *Genomics in the Cloud: Using Docker, GATK, and WDL in Terra*. O'Reilly Media, 2020.

## 相关链接 | References

- [项目主页](https://github.com/lixiang117423/biopytools)
- [GATK Best Practices for Germline SNPs & Indels](https://gatk.broadinstitute.org/hc/en-us/articles/360035890471)
