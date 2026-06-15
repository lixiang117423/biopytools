# vcf2gwas GWAS | vcf2gwas GWAS Analysis Tool

**封装vcf2gwas工具进行全基因组关联分析, 所有参数透传 | vcf2gwas GWAS pipeline with argument passthrough**

## 功能概述 | Overview

vcf2gwas 模块是对第三方工具 vcf2gwas 的轻量封装, 通过 `conda run` 调用 vcf2gwas conda 环境执行分析, 所有命令行参数直接透传给原始工具。vcf2gwas 提供线性混合模型(LMM)和广义线性模型(GLM)关联分析, 内置群体结构校正和亲缘关系矩阵计算, 适合中等到大规模群体的GWAS分析场景。

## 快速开始 | Quick Start

```bash
# 线性混合模型分析
biopytools vcf2gwas -v input.vcf.gz -pf phenotype.csv -p 1 -lmm

# 带协变量
biopytools vcf2gwas -v input.vcf.gz -pf pheno.csv -cf cov.csv -c 1 -p 1 -lmm -P 3

# 透传查看vcf2gwas原始帮助
biopytools vcf2gwas --help
```

## 参数说明 | Parameters

### 包装器选项 | Wrapper Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--vcf2gwas-env` | `vcf2gwas_v.0.8.9` | vcf2gwas conda环境名 |

### vcf2gwas透传参数 | Passthrough Arguments

以下参数直接透传给 vcf2gwas(请参考 vcf2gwas 官方文档):

| 参数 | 描述 |
|------|------|
| `-v` | 输入VCF文件(.vcf.gz) |
| `-pf` | 表型文件(CSV/TSV) |
| `-cf` | 协变量文件(可选) |
| `-p` | 表型列号 |
| `-c` | 协变量列号 |
| `-lmm` | 使用线性混合模型 |
| `-glm` | 使用广义线性模型 |
| `-P` | 主成分数量 |

(运行 `biopytools vcf2gwas --help` 查看 vcf2gwas 完整参数)

## 输出 | Output

输出文件由 vcf2gwas 决定, 典型包括:

- `gwas_results_*.{txt,csv}`: 关联分析结果
- `manhattan_plot_*.png`: 曼哈顿图
- `qq_plot_*.png`: QQ图
- 运行日志文件

## 依赖 | Dependencies

- **vcf2gwas**: conda环境 (https://github.com/frankraden/vcf2gwas)
- **conda**: 环境管理 (推荐 miniforge)
- **GCTA/LIMIX**: vcf2gwas内部依赖

## 安装提示 | Installation Notes

```bash
# 创建vcf2gwas conda环境(需提前安装)
conda create -n vcf2gwas_v.0.8.9 -c bioconda vcf2gwas
```

## 引用 | Citation

- Raden F., et al. vcf2gwas: an all-in-one Python based command line tool for GWAS using VCF files. (https://github.com/frankraden/vcf2gwas)

## 相关链接 | References

- [项目主页](https://github.com/lixiang117423/biopytools)
