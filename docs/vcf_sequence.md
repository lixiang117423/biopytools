# vcf-sequence | 从基因组与 VCF 提取区间序列

**根据 VCF 中每个样本的基因型，将指定基因组区间上的参考序列替换为各样本变异，导出每样本序列 | Reconstruct per-sample sequences of a genomic region by applying VCF variants onto the reference**

## 功能概述 | Overview

本模块解决"在某一基因组区间内，每个样本的精确序列是什么"的问题：以参考基因组 FASTA 为底板，对该区间内每个样本的 SNP/INDEL 变异进行原地替换（参考序列 → 样本等位基因），得到每条样本在该区间的序列。输出可用于构建系统发育树、做候选基因序列比对、提取候选基因区域多序列比对等。

支持选择使用第二个等位基因、不包含参考序列、按最小质量过滤位点、指定或排除样本子集。输出格式可选 `tab`（表格形式，每行一条样本序列）、`fasta`（多序列 FASTA）或 `csv`。

## 快速开始 | Quick Start

```bash
# 提取 chr1:1,000,000-1,001,000 区间的样本序列（FASTA 格式）
biopytools vcf-sequence -v variants.vcf -g genome.fa \
    -c chr1 -s 1000000 -e 1001000 \
    -o seq_output --format fasta
```

## 参数说明 | Parameters

### 必需参数 | Required

| 参数 | 描述 |
|------|------|
| `-v, --vcf` | 输入 VCF 文件路径 |
| `-g, --genome` | 基因组 FASTA 文件路径 |
| `-c, --chrom` | 染色体名称 |
| `-s, --start` | 起始位置 |
| `-e, --end` | 结束位置 |

### 常用可选参数 | Common Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-o, --output-dir` | `./sequence_output` | 输出目录 |
| `--format` | `tab` | 输出格式：`tab` / `fasta` / `csv` |
| `--second-allele` | 关闭 | 使用第二个等位基因（默认使用 GT 中第一个） |
| `--no-reference` | 关闭 | 不包含参考序列 |
| `--min-qual` | - | 最小质量过滤 |
| `--samples` | - | 指定样本（逗号分隔） |
| `--exclude-samples` | - | 排除样本（逗号分隔） |

（运行 `biopytools vcf-sequence -h` 查看完整参数列表）

## 输出 | Output

```
sequence_output/
├── {chrom}_{start}_{end}_sequences.{tab|fasta|csv}   # 各样本序列（含参考序列，除非 --no-reference）
├── {chrom}_{start}_{end}_statistics.txt               # 每样本应用变异数、序列长度
└── extraction_summary.txt                             # 提取汇总
```

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--vcf, -v` | 必填 | Path | 输入VCF文件路径｜Input VCF file path |
| `--genome, -g` | 必填 | Path | 基因组FASTA文件路径｜Genome FASTA file path |
| `--chrom, -c` | 必填 | str | 染色体名称｜Chromosome name |
| `--start, -s` | 必填 | int | 起始位置｜Start position |
| `--end, -e` | 必填 | int | 结束位置｜End position |
| `--output-dir, -o` | `./sequence_output` | Path | 输出目录｜Output directory |
| `--format` | `tab` | tab/fasta/csv | 输出格式｜Output format |
| `--second-allele` | — |  | 使用第二个等位基因｜Use second allele |
| `--no-reference` | — |  | 不包含参考序列｜Do not include reference sequence |
| `--min-qual` | — | int | 最小质量过滤｜Minimum quality filter |
| `--samples` | — | str | 指定样本｜Specify samples |
| `--exclude-samples` | — | str | 排除样本｜Exclude samples |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-v, --vcf` | 必填 |  | 输入VCF文件路径｜Input VCF file path |
| `-g, --genome` | 必填 |  | 基因组FASTA文件路径｜Genome FASTA file path |
| `-c, --chrom` | 必填 |  | 染色体名称｜Chromosome name |
| `-s, --start` | 必填 | int | 起始位置｜Start position |
| `-e, --end` | 必填 | int | 结束位置｜End position |
| `-o, --output-dir` | `./sequence_output` |  | 输出目录｜Output directory |
| `--format` | `tab` | tab/fasta/csv | 输出格式｜Output format |
| `--second-allele` | — | store_true | 使用第二个等位基因｜Use second allele |
| `--no-reference` | — | store_true | 不包含参考序列｜Do not include reference sequence |
| `--min-qual` | — | int | 最小质量值过滤｜Minimum quality filter |
| `--samples` | — |  | 指定样本列表｜Specify samples |
| `--exclude-samples` | — |  | 排除样本列表｜Exclude samples |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- Python：标准库（VCF 解析与 FASTA 索引均使用纯 Python 实现）

## 引用 | Citation

本模块基于 VCF 规范与 FASTA 格式直接处理，请引用：

- Danecek P et al. The variant call format and VCFtools. *Bioinformatics*, 2011, 27(15): 2156-2158. doi:10.1093/bioinformatics/btr330
- Li W, Cowley AC, Uludag M, Gur T, McWilliam H, Squizzato S, Park YM, Buso N, Lopez R. The EMBL-EBI bioinformatics web and programmatic tools framework. *Nucleic Acids Research*, 2015, 43(W1): W580-W584. doi:10.1093/nar/gkv279

## 相关链接 | References

- [项目主页](https://github.com/lixiang117423/biopytools)
- [VCF 规范](https://samtools.github.io/hts-specs/VCFv4.3.pdf)
