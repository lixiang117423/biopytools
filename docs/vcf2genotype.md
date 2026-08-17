# vcf2genotype | VCF 基因型提取工具

**从 VCF 中提取每个样本在每个位点的基因型，输出为 txt/csv/excel 矩阵 | Extract per-sample genotypes from a VCF into txt/csv/excel matrices**

## 功能概述 | Overview

本模块将 VCF 的 GT 字段（如 `0/0`, `0/1`, `1/1`, `./.`）转换为易读的样本 × 位点矩阵，便于在 Excel 或 R/Python 中做群体分析。为了表示更紧凑，基因型被映射为 `0`（参考纯合）、`1`（杂合）、`2`（替代纯合）、`-`（缺失）等简单编码。

采用两遍扫描：第一遍扫描所有位点的 GT 类型，第二遍按确定的列顺序流式写出，避免一次性把整张表载入内存，对大体量 VCF 友好。可选按染色体拆分输出，可选只保留双等位位点，或限定最小/最大变异长度。

## 快速开始 | Quick Start

```bash
# 基本用法
biopytools vcf2genotype -i variants.vcf -o genotypes

# 只保留双等位位点，导出 Excel
biopytools vcf2genotype -i variants.vcf -o genotypes --biallelic-only -t excel

# 按染色体拆分输出
biopytools vcf2genotype -i variants.vcf -o genotypes -e yes
```

## 参数说明 | Parameters

### 必需参数 | Required

| 参数 | 描述 |
|------|------|
| `-i, --input` | 输入 VCF 文件路径 |

### 常用可选参数 | Common Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-o, --output` | `vcf2genotype` | 输出文件前缀 |
| `-s, --samples` | `all` | 样本选择：`all` 或逗号分隔的样本名 |
| `-e, --each` | `n` | 是否按染色体拆分输出：`yes/y` 或 `no/n` |
| `-t, --output-type` | `txt` | 输出格式：`txt` / `csv` / `excel` |
| `--output-dir` | `./` | 输出目录 |
| `--biallelic-only` | 关闭 | 只保留双等位位点 |
| `--min-length` | - | 最小变异长度（含） |
| `--max-length` | - | 最大变异长度（含） |
| `-v` (重复) | 关闭 | 增加输出详细程度 |
| `--quiet` | 关闭 | 静默模式 |
| `--log-level` | `INFO` | 日志级别（DEBUG/INFO/WARNING/ERROR/CRITICAL） |
| `--log-file` | - | 日志文件路径 |
| `--dry-run` | 关闭 | 试运行模式 |

（运行 `biopytools vcf2genotype -h` 查看完整参数列表）

## 输出 | Output

```
./
├── {prefix}.txt                 # 主基因型矩阵（或 .csv / .xlsx）
├── {prefix}_chr1.txt            # 按染色体拆分输出（-e yes 时）
├── {prefix}_chr2.txt
└── ...
```

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input, -i` | 必填 | Path | 输入VCF文件路径｜Input VCF file path |
| `--output, -o` | `vcf2genotype` | str | 输出文件前缀｜Output file prefix |
| `--samples, -s` | `all` | str | 样本选择：all(所有样本)或逗号分隔的样本名称｜Sample selection: all or comma-separated sample names |
| `--biallelic-only` | — |  | 只保留双等位位点｜Keep only biallelic sites |
| `--min-length` | — | int | 最小变异长度（包含），默认不限制｜Minimum variant length (inclusive), default no limit |
| `--max-length` | — | int | 最大变异长度（包含），默认不限制｜Maximum variant length (inclusive), default no limit |
| `--each, -e` | `n` | yes/y/no/n | 按染色体拆分输出文件：yes/y或no/n｜Split output files by chromosome: yes/y or no/n |
| `--output-type, -t` | `txt` | txt/csv | 输出文件格式｜Output file format |
| `--output-dir` | `./` | Path | 输出目录｜Output directory |
| `--verbose, -v` | — |  | 增加输出详细程度｜Increase output verbosity |
| `--quiet` | — |  | 静默模式｜Quiet mode |
| `--log-file` | — | str | 日志文件路径｜Log file path |
| `--log-level` | `INFO` | DEBUG/INFO/WARNING/ERROR/CRITICAL | 日志级别｜Log level |
| `--dry-run` | — |  | 试运行模式｜Dry run mode |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | VCF文件路径(支持.gz压缩格式)｜VCF file path (supports .gz compressed format) |
| `-o, --output` | `vcf_genotype` |  | 输出文件前缀｜Output file prefix |
| `-t, --output-type` | `txt` | txt/csv | 输出文件格式｜Output file format |
| `--output-dir` | `./` |  | 输出目录｜Output directory |
| `-s, --samples` | `all` |  | 样本选择：all（所有样本）或逗号分隔的样本名称｜Sample selection: all (all samples) or comma-separated sample names |
| `--biallelic-only` | — | store_true | 只保留双等位位点｜Keep only biallelic sites |
| `--min-length` | — | int | 最小变异长度（包含），默认不限制｜Minimum variant length (inclusive), default no limit |
| `--max-length` | — | int | 最大变异长度（包含），默认不限制｜Maximum variant length (inclusive), default no limit |
| `-e, --each` | `n` | yes/y/no/n | 按染色体拆分输出文件：yes/y（是）或no/n（否）｜Split output files by chromosome: yes/y or no/n |
| `--dry-run` | — | store_true | 试运行模式，不实际执行操作｜Dry run mode, no actual operations performed |
| `-v, --verbose` | `0` | count | 增加输出详细程度 (-v, -vv, -vvv)｜Increase output verbosity (-v, -vv, -vvv) |
| `--quiet` | — | store_true | 静默模式，仅输出错误信息｜Quiet mode, only output error messages |
| `--log-file` | — | str | 日志文件路径｜Log file path |
| `--log-level` | `INFO` | DEBUG/INFO/WARNING/ERROR/CRITICAL | 日志级别｜Log level |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- Python：标准库
- 可选：[cyvcf2](https://github.com/brentp/cyvcf2)（若安装则自动启用更快的 VCF 解析后端）
- 可选：`openpyxl`（输出 Excel 格式时）

## 引用 | Citation

- Pedersen BS, Quinlan AR. cyvcf2: fast, flexible variant analysis with Python. *Bioinformatics*, 2017, 33(12): 1867-1869. doi:10.1093/bioinformatics/btx057（若启用 cyvcf2 后端）
- Danecek P et al. The variant call format and VCFtools. *Bioinformatics*, 2011, 27(15): 2156-2158. doi:10.1093/bioinformatics/btr330

## 相关链接 | References

- [项目主页](https://github.com/lixiang117423/biopytools)
- [VCF 规范](https://samtools.github.io/hts-specs/VCFv4.3.pdf)
