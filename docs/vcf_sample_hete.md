# vcf-sample-hete | VCF 基因型统计分析

**对 VCF 文件中每个样本统计杂合率、缺失率、纯合率与深度/质量分布 | Per-sample genotype statistics (heterozygosity, missingness, homozygosity, depth/quality) from a VCF**

## 功能概述 | Overview

本模块面向群体重测序质控（QC）：逐样本扫描 VCF，统计每位样本的位点总数、杂合位点数、参考/替代纯合位点数、缺失基因型数与比例，以及深度和质量分数的分布。结果用于发现异常样本（如异常高杂合率提示污染、高缺失率提示低覆盖度或样本质量问题）。

输出分为"详细统计"（每样本逐指标一行）与"汇总统计"（全体均值/中位数/分位数），两者均可单独关闭。支持按最小深度与最小质量过滤，可选排除缺失基因型后再统计杂合率。

## 快速开始 | Quick Start

```bash
# 基本统计
biopytools vcf-sample-hete -i variants.vcf -o vcf_stats_output

# 排除缺失基因型，并要求最小深度 5
biopytools vcf-sample-hete -i variants.vcf -o vcf_stats_output -e -d 5
```

## 参数说明 | Parameters

### 必需参数 | Required

| 参数 | 描述 |
|------|------|
| `-i, --input` | 输入 VCF 文件路径 |

### 常用可选参数 | Common Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-o, --output` | `vcf_stats_output` | 输出目录 |
| `-d, --min-depth` | `0` | 最小深度过滤阈值 |
| `-q, --min-qual` | `0.0` | 最小质量分数过滤阈值 |
| `-e, --exclude-missing` | 关闭 | 排除缺失基因型 |
| `-D, --no-detailed` | 关闭 | 不输出详细统计结果 |
| `-S, --no-summary` | 关闭 | 不输出汇总统计结果 |
| `--verbose` (重复) | 关闭 | 增加输出详细程度 |
| `--quiet` | 关闭 | 静默模式 |
| `--log-level` | `INFO` | 日志级别 |
| `--log-file` | - | 日志文件路径 |
| `--dry-run` | 关闭 | 试运行模式 |

（运行 `biopytools vcf-sample-hete -h` 查看完整参数列表）

## 输出 | Output

```
vcf_stats_output/
├── detailed_stats.txt     # 每样本逐项统计（位点总数、Ho、He、缺失率、深度均值等）
├── summary_stats.txt      # 全体样本汇总（均值、中位数、分位数）
└── *.log                  # 日志
```

## 依赖 | Dependencies

- Python：标准库

## 引用 | Citation

本模块基于 VCF 规范直接解析基因型字段，请引用 VCF 与下游分析工具：

- Danecek P et al. The variant call format and VCFtools. *Bioinformatics*, 2011, 27(15): 2156-2158. doi:10.1093/bioinformatics/btr330
- Purcell S et al. PLINK: a tool set for whole-genome association and population-based linkage analyses. *American Journal of Human Genetics*, 2007, 81(3): 559-575. doi:10.1086/519795

## 相关链接 | References

- [项目主页](https://github.com/lixiang117423/biopytools)
- [VCF 规范](https://samtools.github.io/hts-specs/VCFv4.3.pdf)
