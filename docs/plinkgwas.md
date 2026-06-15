# PLINK GWAS | PLINK Genome-Wide Association Study

**基于PLINK的全基因组关联分析, 含质控、PCA校正和多种显著性校正 | GWAS pipeline with QC, PCA correction and multiple testing correction**

## 功能概述 | Overview

plinkgwas 模块封装了 PLINK 进行全基因组关联分析(GWAS), 提供端到端的完整流程: 质控(缺失率/MAF/HWE) → LD pruning → PCA计算 → 关联分析(含PCA校正以处理群体结构) → 多重检验校正(Bonferroni/FDR/提示性阈值) → 结果汇总。支持质量性状(病例/对照)和数量性状两种表型类型, 以及加性、显性、隐性等遗传模型。

## 快速开始 | Quick Start

```bash
# 标准病例/对照GWAS
biopytools plink-gwas -i data.vcf.gz -p pheno.txt -o results

# 数量性状, 自定义阈值
biopytools plink-gwas -i data.vcf.gz -p pheno.txt -T quantitative -o results --maf 0.01 --suggestive-threshold 5e-8
```

## 参数说明 | Parameters

### 必需参数 | Required

| 参数 | 描述 |
|------|------|
| `-i, --vcf` | VCF文件路径 |
| `-p, --phenotype` | 表型文件路径 |

### 常用可选参数 | Common Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-T, --trait-type` | `qualitative` | 表型类型(qualitative/quantitative) |
| `-m, --genetic-model` | `additive` | 遗传模型(additive/dominant/recessive/all) |
| `-o, --output-dir` | `gwas_results` | 输出目录 |
| `--no-strat-corr` | `False` | 禁用群体结构校正 |
| `--mind` | `None` | 个体缺失率阈值 |
| `--geno` | `0.05` | SNP缺失率阈值 |
| `--maf` | `0.05` | 最小等位基因频率 |
| `--hwe` | `None` | HWE平衡p值阈值 |
| `--ld-window-size` | `50` | LD窗口大小(kb) |
| `--ld-step-size` | `5` | LD步长(SNPs) |
| `--ld-r2-threshold` | `0.2` | LD r2阈值 |
| `--pca-components` | `10` | 计算的主成分数量 |
| `--pca-use` | `5` | 关联分析中使用的主成分数量 |
| `--correction-method` | `all` | 显著性校正(bonferroni/suggestive/fdr/all) |
| `--bonferroni-alpha` | `0.05` | Bonferroni alpha水平 |
| `--suggestive-threshold` | `1e-5` | 提示性关联阈值 |
| `--fdr-alpha` | `0.05` | FDR q值阈值 |
| `-t, --threads` | `12` | 线程数 |
| `-f, --force` | `False` | 强制覆盖已存在的输出目录 |
| `--dry-run` | `False` | 模拟运行不实际执行分析 |

(运行 `biopytools plink-gwas -h` 查看完整参数列表)

## 输出 | Output

```
gwas_results/
├── qc/                      # 质控结果(filtered.{bed,bim,fam})
├── pca/                     # PCA结果
│   ├── pca.eigenvec
│   └── pca.eigenval
├── association/             # 关联分析结果
│   └── gwas.assoc.logistic  # 或 .linear (数量性状)
├── corrected_results/       # 多重检验校正后结果
│   └── significant_snps.txt
└── gwas.log                 # 运行日志
```

## 依赖 | Dependencies

- **PLINK 1.9/2.0**: GWAS计算 (https://www.cog-genomics.org/plink/)
- **Python库**: pandas, numpy, scipy (校正计算)

## 引用 | Citation

- Purcell S., et al. (2007) PLINK: a toolset for whole-genome association and population-based linkage analyses. AJHG. 81(3):559-575.
- Chang C.C., Chow C.C., Tellier L.C., et al. (2015) Second-generation PLINK: rising to the challenge of larger and richer datasets. GigaScience. 4(1):7.

## 相关链接 | References

- [项目主页](https://github.com/lixiang117423/biopytools)
