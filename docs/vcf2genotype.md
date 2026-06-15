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
