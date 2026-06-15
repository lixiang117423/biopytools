# GWAS候选基因 | GWAS Candidate Gene Finder (gwas2gene)

**根据GWAS显著SNP从GFF注释中提取窗口内候选基因, 可整合功能注释 | Extract candidate genes near significant GWAS SNPs**

## 功能概述 | Overview

gwas2gene 模块用于从GWAS分析结果中识别显著SNP(P值低于阈值), 并在用户指定的上下游窗口内从GFF3注释文件提取候选基因。适用于GWAS下游功能基因挖掘, 将统计学显著位点映射到具体的基因ID, 并可选整合外部功能注释文件(如GO/KEGG/基因描述)以辅助生物学解释。

## 快速开始 | Quick Start

```bash
# 标准用法
biopytools gwas2gene -g gwas.txt -p Pvalue -a annotation.gff3 -o candidates.tsv

# 自定义阈值和窗口
biopytools gwas2gene -g gwas.txt -p Pvalue -t 1e-6 -w 100000 -a annotation.gff3 -o candidates.tsv

# 带功能注释
biopytools gwas2gene -g gwas.txt -p Pvalue -a annotation.gff3 -f gene_func.tsv -o candidates.tsv
```

## 参数说明 | Parameters

### 必需参数 | Required

| 参数 | 描述 |
|------|------|
| `-g, --gwas` | GWAS结果文件路径 |
| `-p, --pval-col` | P值所在列名或列号(1-based) |
| `--gff` | GFF3注释文件路径 |
| `-o, --output` | 输出文件路径 |

### 常用可选参数 | Common Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-t, --threshold` | `1e-5` | P值阈值(显著性SNP判定) |
| `-w, --window` | `200000` | 上下游窗口大小(bp) |
| `-f, --func` | `None` | 功能注释文件路径(可选) |

(运行 `biopytools gwas2gene -h` 查看完整参数列表)

## 输出 | Output

输出TSV文件, 包含显著SNP信息和对应候选基因:

| 列名 | 描述 |
|------|------|
| Chr | 染色体 |
| Pos | SNP位置 |
| Pvalue | GWAS P值 |
| Gene_ID | 候选基因ID |
| Gene_Start | 基因起始位置 |
| Gene_End | 基因终止位置 |
| Strand | 链方向 |
| Distance | SNP到基因的距离 |
| Function | 基因功能(若提供`--func`) |

## 依赖 | Dependencies

- **Python库**: pandas (文件处理), 无需外部生物信息学工具

## 引用 | Citation

无外部算法依赖, 为本工具集实现的辅助脚本。GWAS分析方法请参考相应工具的引用(如PLINK/GEMMA/EMMAX)。

## 相关链接 | References

- [项目主页](https://github.com/lixiang117423/biopytools)
