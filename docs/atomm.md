# ATOMM关联分析 | ATOMM Two-Organism Mixed Model Association

**基于双物种(宿主-病原)混合效应模型的交叉感染关联分析 | Cross-infection association analysis using two-organism mixed-effect model**

## 功能概述 | Overview

atomm 模块实现 ATOMM (A Two-Organism Mixed Model) 算法, 用于宿主-病原体(host-pathogen)交叉感染设计下的全基因组关联分析。传统GWAS仅处理单物种基因型-表型关系, 而 ATOMM 同时建模宿主和病原体双方的SNP效应及交互效应, 捕捉双方遗传变异对感染表型(如病情指数、繁殖适合度)的联合贡献, 适用于植物-病原、宿主-共生体等双向遗传学研究。

## 快速开始 | Quick Start

```bash
# 标准分析
biopytools atomm --host-vcf host.vcf.gz --pathogen-vcf pathogen.vcf.gz \
    --phenotype-matrix pheno_matrix.csv -o atomm_output

# 自定义MAF和编码方式
biopytools atomm --host-vcf host.vcf.gz --pathogen-vcf pathogen.vcf.gz \
    --phenotype-matrix pheno_matrix.csv -o atomm_output --maf 0.01 --encoding dosage
```

## 参数说明 | Parameters

### 必需参数 | Required

| 参数 | 描述 |
|------|------|
| `--host-vcf` | 宿主VCF文件(支持.gz) |
| `--pathogen-vcf` | 病原VCF文件(支持.gz) |
| `--phenotype-matrix` | 交叉感染表型矩阵(行=宿主, 列=病原) |

### 常用可选参数 | Common Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-o, --output-dir` | `./output` | 输出目录 |
| `--maf` | `0.05` | MAF过滤阈值 |
| `--encoding` | `auto` | 基因型编码方式(auto/haploid/dosage) |
| `--convert-maf` | `0.05` | VCF转换时MAF阈值 |
| `--missing-value` | `NA` | 表型缺失值标记 |
| `--host-snp-range` | `None` | 宿主边际检验SNP范围(起 止) |
| `--pathogen-snp-range` | `None` | 病原边际检验SNP范围(起 止) |
| `--interaction-host-range` | `None` | 交互检验宿主SNP范围 |
| `--interaction-pathogen-range` | `None` | 交互检验病原SNP范围 |
| `--tol` | `1e-6` | 优化容忍度 |
| `--maxiter` | `10000` | 最大迭代次数 |

(运行 `biopytools atomm -h` 查看完整参数列表)

## 输出 | Output

```
output/
├── host_marginal.csv          # 宿主SNP边际效应
├── pathogen_marginal.csv      # 病原SNP边际效应
├── interaction.csv            # 宿主-病原SNP交互效应
├── kinship_host.png           # 宿主亲缘关系矩阵
├── kinship_pathogen.png       # 病原亲缘关系矩阵
└── atomm.log                  # 运行日志
```

## 依赖 | Dependencies

- **Python库**: numpy, scipy (优化算法), pandas, statsmodels
- **BCFtools**: VCF预处理(可选)
- 推荐 Python 3.10+

## 输入格式说明 | Input Format

**表型矩阵 (CSV)**: 行名为宿主样本ID, 列名为病原样本ID, 单元格为对应交叉感染的表型值。缺失值用`NA`或自定义字符。

## 引用 | Citation

- Wang L., Paul P., Chen L., et al. (2018) ATOMM: a two-organism mixed model for mapping cross-kingdom genetic interactions in host-microbe associations. PNAS. 115(19):E4350-E4359.

## 相关链接 | References

- [项目主页](https://github.com/lixiang117423/biopytools)
