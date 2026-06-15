# VCF PCA分析 | VCF PCA Principal Component Analysis

**基于PLINK的VCF主成分分析, 含质控和可视化 | Principal component analysis with quality control and visualization**

## 功能概述 | Overview

vcf_pca 模块用于对VCF变异文件进行主成分分析(PCA), 将高维基因型数据降维到少数几个主成分, 用于揭示群体遗传结构、识别异常个体、检测批次效应。流程默认包含完整质控(MAF/缺失率/HWE), 使用 PLINK 进行高效PCA计算, 并可选生成按分组着色的可视化图表。

## 快速开始 | Quick Start

```bash
# 标准PCA分析
biopytools vcf-pca -i variants.vcf -o pca_results

# 带样本信息和可视化
biopytools vcf-pca -i variants.vcf -o pca_results -s sample_info.txt -g population -P
```

## 参数说明 | Parameters

### 必需参数 | Required

| 参数 | 描述 |
|------|------|
| `-i, --input` | 输入VCF文件 |

### 常用可选参数 | Common Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-o, --output` | `./pca_output` | 输出目录 |
| `-s, --sample-info` | `None` | 样本信息文件 |
| `-c, --components` | `10` | 主成分数量 |
| `--maf` | `0.05` | 最小等位基因频率阈值 |
| `--missing` | `0.1` | 最大缺失率阈值 |
| `--hwe` | `1e-6` | Hardy-Weinberg平衡p值阈值 |
| `--skip-qc` | `False` | 跳过质量控制过滤 |
| `-P, --plot` | `False` | 生成PCA可视化图表 |
| `-g, --group-column` | `None` | 分组列名 |
| `--plink-path` | `plink` | PLINK软件路径 |
| `--bcftools-path` | `bcftools` | BCFtools软件路径 |

(运行 `biopytools vcf-pca -h` 查看完整参数列表)

## 输出 | Output

```
pca_output/
├── filtered.{bed,bim,fam}      # 质控后的PLINK二进制文件
├── pca.eigenvec                # 主成分得分(每个样本各PC)
├── pca.eigenval                # 特征值(各PC解释方差)
├── pca_plot.pdf/png            # PCA散点图(需 -P)
└── pca.log                     # 运行日志
```

## 依赖 | Dependencies

- **PLINK 1.9/2.0**: 质控和PCA计算 (https://www.cog-genomics.org/plink/)
- **BCFtools**: VCF处理 (http://www.htslib.org/)
- **Python库**: pandas, matplotlib (用于可视化)

## 引用 | Citation

- Purcell S., et al. (2007) PLINK: a toolset for whole-genome association and population-based linkage analyses. AJHG. 81(3):559-575.
- Patterson N., Price A.L., Reich D. (2006) Population structure and eigenanalysis. PLoS Genetics. 2(12):e190.

## 相关链接 | References

- [项目主页](https://github.com/lixiang117423/biopytools)
