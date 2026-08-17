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

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input, -i` | 必填 |  | 输入VCF文件路径(支持.vcf/.vcf.gz)｜Input VCF file path (supports .vcf/.vcf.gz) |
| `--output-dir, -o` | `./filtered_vcf` |  | 输出目录｜Output directory |
| `--threads, -t` | `12` | int | 线程数｜Number of threads |
| `--snp-qual` | `30.0` | float | SNP最小质量值｜SNP minimum QUAL |
| `--snp-dp` | `10` | int | SNP最小测序深度｜SNP minimum DP |
| `--snp-mq` | `40.0` | float | SNP最小比对质量｜SNP minimum MQ |
| `--snp-qd` | `2.0` | float | SNP最小质量/深度比｜SNP minimum QD |
| `--snp-fs` | `60.0` | float | SNP最大FisherStrand值｜SNP maximum FS |
| `--snp-sor` | `3.0` | float | SNP最大StrandOddsRatio｜SNP maximum SOR |
| `--snp-mqrs` | `-12.5` | float | SNP最小MappingQualityRankSum｜SNP minimum MQRankSum |
| `--snp-rprs` | `-8.0` | float | SNP最小ReadPosRankSum｜SNP minimum ReadPosRankSum |
| `--snp-maf` | `0.05` | float | SNP最小次等位基因频率｜SNP minimum MAF |
| `--snp-biallelic` | `True` |  | 只保留双等位位点SNP｜Keep only biallelic SNPs |
| `--no-snp-biallelic` | `False` |  | 禁用双等位点过滤｜Disable biallelic filtering |
| `--indel-qual` | `30.0` | float | INDEL最小质量值｜INDEL minimum QUAL |
| `--indel-dp` | `10` | int | INDEL最小测序深度｜INDEL minimum DP |
| `--indel-mq` | `40.0` | float | INDEL最小比对质量｜INDEL minimum MQ |
| `--indel-qd` | `2.0` | float | INDEL最小质量/深度比｜INDEL minimum QD |
| `--indel-fs` | `200.0` | float | INDEL最大FisherStrand值｜INDEL maximum FS |
| `--indel-sor` | `10.0` | float | INDEL最大StrandOddsRatio｜INDEL maximum SOR |
| `--indel-rprs` | `-20.0` | float | INDEL最小ReadPosRankSum｜INDEL minimum ReadPosRankSum |
| `--variant-type` | `both` | both/snp_only/indel_only | 输入VCF文件的变异类型｜Variant type in input VCF |
| `--bcftools-path` | `bcftools` |  | BCFtools软件路径｜BCFtools software path |
| `--verbose, -v` | — |  | 详细输出模式(-v: INFO, -vv: DEBUG)｜Verbose mode (-v: INFO, -vv: DEBUG) |
| `--quiet` | — |  | 静默模式(仅输出ERROR)｜Quiet mode (ERROR only) |
| `--log-level` | — |  | 日志级别(DEBUG/INFO/WARNING/ERROR/CRITICAL)｜Log level |
| `--force, -f` | — |  | 强制覆盖已存在文件｜Force overwrite existing files |
| `--dry-run` | — |  | 模拟运行(不实际执行)｜Dry run without execution |
| `--repair-vcf` | — |  | 自动修复损坏的VCF文件（列数不匹配等问题）｜Auto-repair corrupted VCF files (column mismatch, etc.) |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 输入VCF文件路径(支持压缩和未压缩)｜Input VCF file path (supports compressed and uncompressed) |
| `-o, --output-dir` | `./filtered_vcf` |  | 输出目录｜Output directory (default: ./filtered_vcf) |
| `-t, --threads` | `12` | int | 线程数｜Number of threads (default: 12) |
| `--variant-type` | `both` | both/snp_only/indel_only | 输入VCF文件的变异类型｜Variant type in input VCF (default: both) |
| `--snp-qual` | `30.0` | float | SNP最小质量值｜SNP minimum QUAL (default: 30.0) |
| `--snp-dp` | `10` | int | SNP最小测序深度｜SNP minimum DP (default: 10) |
| `--snp-mq` | `40.0` | float | SNP最小比对质量｜SNP minimum MQ (default: 40.0) |
| `--snp-qd` | `2.0` | float | SNP最小质量/深度比｜SNP minimum QD (default: 2.0) |
| `--snp-fs` | `60.0` | float | SNP最大FisherStrand值｜SNP maximum FS (default: 60.0) |
| `--snp-sor` | `3.0` | float | SNP最大StrandOddsRatio｜SNP maximum SOR (default: 3.0) |
| `--snp-mqrs` | `-12.5` | float | SNP最小MappingQualityRankSum｜SNP minimum MQRankSum (default: -12.5) |
| `--snp-rprs` | `-8.0` | float | SNP最小ReadPosRankSum｜SNP minimum ReadPosRankSum (default: -8.0) |
| `--snp-maf` | `0.05` | float | SNP最小次等位基因频率｜SNP minimum MAF (default: 0.05) |
| `--snp-biallelic` | `True` | store_true | 只保留双等位位点SNP｜Keep only biallelic SNPs (default: True) |
| `--no-snp-biallelic` | — | store_false | 禁用双等位位点过滤｜Disable biallelic filtering |
| `--indel-qual` | `30.0` | float | INDEL最小质量值｜INDEL minimum QUAL (default: 30.0) |
| `--indel-dp` | `10` | int | INDEL最小测序深度｜INDEL minimum DP (default: 10) |
| `--indel-mq` | `40.0` | float | INDEL最小比对质量｜INDEL minimum MQ (default: 40.0) |
| `--indel-qd` | `2.0` | float | INDEL最小质量/深度比｜INDEL minimum QD (default: 2.0) |
| `--indel-fs` | `200.0` | float | INDEL最大FisherStrand值｜INDEL maximum FS (default: 200.0) |
| `--indel-sor` | `10.0` | float | INDEL最大StrandOddsRatio｜INDEL maximum SOR (default: 10.0) |
| `--indel-rprs` | `-20.0` | float | INDEL最小ReadPosRankSum｜INDEL minimum ReadPosRankSum (default: -20.0) |
| `--bcftools-path` | `bcftools` |  | BCFtools软件路径｜BCFtools software path (default: bcftools) |
| `-v, --verbose` | `0` | count | 详细输出模式(-v: INFO, -vv: DEBUG)｜Verbose mode (-v: INFO, -vv: DEBUG) |
| `--quiet` | — | store_true | 静默模式(只输出ERROR)｜Quiet mode (ERROR only) |
| `--log-level` | — |  | 日志级别(DEBUG/INFO/WARNING/ERROR/CRITICAL)｜Log level (default: INFO) |
| `-f, --force` | — | store_true | 强制覆盖已存在文件｜Force overwrite existing files |
| `--dry-run` | — | store_true | 模拟运行(不实际执行)｜Dry run without execution |
| `--repair-vcf` | — | store_true | 自动修复损坏的VCF文件（列数不匹配等问题）｜Auto-repair corrupted VCF files (column mismatch, etc.) |
| `-V, --version` | — | version |  |

<!-- END PARAMS:auto -->

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
