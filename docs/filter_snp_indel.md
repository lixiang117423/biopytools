# filter-snp-indel | VCF SNP/INDEL 过滤工具

**基于 GATK 最佳实践对 VCF 中的 SNP 与 INDEL 分别应用质控阈值进行过滤 | GATK best-practice VCF filtering with separate thresholds for SNPs and INDELs**

## 功能概述 | Overview

本模块封装 BCFtools，对 GATK 或其他变异检测流程产出的 VCF 文件进行硬过滤（hard filtering）。它遵循 GATK 推荐做法，自动将 SNP 与 INDEL 分离，分别应用不同的质控指标阈值（QD、FS、MQ、SOR、MQRankSum、ReadPosRankSum、DP、QUAL、MAF 等），过滤后再合并。同时也支持只处理 SNP 或只处理 INDEL。

默认参数来自 GATK 官方推荐的硬过滤阈值，可作为常规过滤的起点；用户也可针对低深度或高深度数据集自行调整。模块还提供 VCF 自动修复（列数不匹配等）和模拟运行（dry-run）能力。

## 快速开始 | Quick Start

```bash
# 使用默认 GATK 阈值过滤
biopytools filter-snp-indel -i variants.vcf.gz -o filtered_output/

# 自定义 SNP 与 INDEL 阈值
biopytools filter-snp-indel -i variants.vcf.gz \
    --snp-qd 3.0 --snp-fs 50 --indel-fs 200 \
    -o filtered_output/
```

## 参数说明 | Parameters

### 必需参数 | Required

| 参数 | 描述 |
|------|------|
| `-i, --input` | 输入 VCF 文件路径（支持 `.vcf` / `.vcf.gz`） |

### 常用可选参数 | Common Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-o, --output-dir` | `./filtered_vcf` | 输出目录 |
| `-t, --threads` | `12` | 线程数 |
| `--variant-type` | `both` | 输入 VCF 的变异类型：`both` / `snp_only` / `indel_only` |
| `--bcftools-path` | `bcftools` | BCFtools 可执行文件路径 |
| `--repair-vcf` | 关闭 | 自动修复损坏的 VCF 文件（列数不匹配等） |
| `-f, --force` | 关闭 | 强制覆盖已存在文件 |
| `--dry-run` | 关闭 | 模拟运行，不实际执行命令 |
| `-v` (重复) | INFO | 详细输出（`-v` INFO，`-vv` DEBUG） |
| `--quiet` | 关闭 | 仅输出 ERROR 日志 |
| `--log-level` | - | 日志级别（DEBUG/INFO/WARNING/ERROR/CRITICAL） |
| `--log-file` | - | 日志文件路径 |

### SNP 过滤阈值 | SNP Thresholds

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--snp-qual` | `30.0` | 最小 QUAL |
| `--snp-dp` | `10` | 最小测序深度 |
| `--snp-mq` | `40.0` | 最小比对质量 |
| `--snp-qd` | `2.0` | 最小 QD（QUAL/DP） |
| `--snp-fs` | `60.0` | 最大 FisherStrand |
| `--snp-sor` | `3.0` | 最大 StrandOddsRatio |
| `--snp-mqrs` | `-12.5` | 最小 MQRankSum |
| `--snp-rprs` | `-8.0` | 最小 ReadPosRankSum |
| `--snp-maf` | `0.05` | 最小次等位基因频率 |
| `--snp-biallelic` | 开启 | 只保留双等位位点（`--no-snp-biallelic` 关闭） |

### INDEL 过滤阈值 | INDEL Thresholds

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--indel-qual` | `30.0` | 最小 QUAL |
| `--indel-dp` | `10` | 最小测序深度 |
| `--indel-mq` | `40.0` | 最小比对质量 |
| `--indel-qd` | `2.0` | 最小 QD |
| `--indel-fs` | `200.0` | 最大 FisherStrand |
| `--indel-sor` | `10.0` | 最大 StrandOddsRatio |
| `--indel-rprs` | `-20.0` | 最小 ReadPosRankSum |

（运行 `biopytools filter-snp-indel -h` 查看完整参数列表）

## 输出 | Output

```
filtered_vcf/
├── {base}.filtered.snp.vcf.gz      # 过滤后的 SNP
├── {base}.filtered.indel.vcf.gz    # 过滤后的 INDEL
├── {base}.filtered.vcf.gz          # 合并后的最终 VCF（variant_type=both 时）
├── {base}.filtered.snp.biallelic.vcf.gz  # 双等位 SNP（若启用）
└── filtering_report.txt            # 过滤统计报告
```

## 依赖 | Dependencies

- [BCFtools](http://www.htslib.org/)（必需，用于变异分离、过滤与合并）
- Python：标准库（无额外第三方依赖）

## 引用 | Citation

模块仅封装 BCFtools，请引用：

- Li H. A statistical framework for SNP calling, mutation discovery, association mapping and population genetical parameter estimation from sequencing data. *Bioinformatics*, 2011, 27(21): 2987-2993. doi:10.1093/bioinformatics/btr509
- GATK 硬过滤建议：Van der Auwera GA & O'Connor BD. *Genomics in the Cloud*. O'Reilly, 2020.

## 相关链接 | References

- [项目主页](https://github.com/lixiang117423/biopytools)
- [BCFtools 文档](http://www.htslib.org/doc/bcftools.html)
- [GATK 硬过滤指南](https://gatk.broadinstitute.org/hc/en-us/articles/360035890471)
