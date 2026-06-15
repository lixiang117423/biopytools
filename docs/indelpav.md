# indelpav | INDEL 存在/缺失变异分析

**从 VCF 提取 INDEL 并构建样本 × INDEL 的存在/缺失（Presence/Absence）矩阵 | Build a sample-by-INDEL presence/absence matrix from a VCF file**

## 功能概述 | Overview

本模块针对群体重测序数据，将 VCF 中的 INDEL 变异转换为 PAV（Presence/Absence Variation）矩阵：对每个样本每个 INDEL 位点，根据其等位基因状态判定为"存在"或"缺失"。这种矩阵适合做 PCA、聚类、进化树等不以精确基因型为单位的下游分析，对长度较大或复杂性较高的 INDEL 尤其有用。

流程通过 BCFtools 提取 INDEL、按用户指定的长度与质量阈值过滤，再扫描 GT 字段生成 PAV 表与汇总报告。可选包含复杂变异，输出可压缩。

## 快速开始 | Quick Start

```bash
biopytools indelpav -v variants.vcf.gz -o indel_pav.txt
```

## 参数说明 | Parameters

### 必需参数 | Required

| 参数 | 描述 |
|------|------|
| `-v, --vcf` | 输入 VCF 文件路径 |

### 常用可选参数 | Common Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-o, --output` | `./indel_pav.txt` | 输出文件路径 |
| `-t, --threads` | `12` | 线程数 |
| `--min-length` | `1` | 最小 INDEL 长度（bp） |
| `--max-length` | - | 最大 INDEL 长度（bp） |
| `-q, --min-quality` | `20.0` | 最小质量分数 |
| `-d, --min-depth` | `5` | 最小深度 |
| `--max-missing` | `0.8` | 最大缺失率（0-1） |
| `--include-complex` | 关闭 | 包含复杂变异 |
| `--compress` | 关闭 | 压缩输出文件 |
| `--bcftools-path` | `bcftools` | BCFtools 可执行文件路径 |

（运行 `biopytools indelpav -h` 查看完整参数列表）

## 输出 | Output

```
indel_pav.txt              # 样本 × INDEL PAV 矩阵（1/0/-）
indel_pav_summary.txt      # 位点数、样本数、过滤统计汇总
```

## 依赖 | Dependencies

- [BCFtools](http://www.htslib.org/)（必需，用于 INDEL 提取与质控过滤）
- Python：标准库

## 引用 | Citation

本模块为基于 BCFtools 与 VCF 规范的下游分析封装，请引用：

- Li H. A statistical framework for SNP calling, mutation discovery, association mapping and population genetical parameter estimation from sequencing data. *Bioinformatics*, 2011, 27(21): 2987-2993. doi:10.1093/bioinformatics/btr509
- Danecek P et al. The variant call format and VCFtools. *Bioinformatics*, 2011, 27(15): 2156-2158. doi:10.1093/bioinformatics/btr330

## 相关链接 | References

- [项目主页](https://github.com/lixiang117423/biopytools)
- [BCFtools 文档](http://www.htslib.org/doc/bcftools.html)
