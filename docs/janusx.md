# JanusX GWAS和基因组选择分析模块

**全基因组关联分析和基因组选择的统一工具包 | Unified Toolkit for Genome-Wide Association Study and Genomic Selection**

## 功能概述 | Overview

JanusX分析模块封装了JanusX软件，提供GWAS(全基因组关联分析)和GS(基因组选择)的完整分析流程。支持多种GWAS模型(GLM、LMM、fastLMM、FarmCPU)和GS模型(GBLUP、rrBLUP、BayesA/B/Cpi)，以及PCA分析和PostGWAS可视化功能。

## 主要特性 | Key Features

- **多种GWAS模型**: GLM、LMM、fastLMM、FarmCPU
- **多种GS模型**: GBLUP、rrBLUP、BayesA、BayesB、BayesCpi
- **PCA分析**: 主成分分析和可视化
- **PostGWAS**: Manhattan图、QQ图和SNP注释
- **多格式支持**: VCF和PLINK格式
- **高性能**: 基于Rust编译的高性能可执行文件

## 快速开始 | Quick Start

### GWAS分析 | GWAS Analysis

```bash
# 基本GWAS分析
biopytools janusx -i genotypes.vcf.gz -p phenotypes.txt -m gwas -o gwas_results

# 使用LMM模型
python -m biopytools.janusx gwas -i genotypes.vcf.gz -p phenotypes.txt -m lmm -o gwas_results

# 多模型分析
python -m biopytools.janusx gwas -i genotypes.vcf.gz -p phenotypes.txt -m lmm lm fastlmm -o gwas_results
```

### GS分析 | Genomic Selection

```bash
# 基本GS分析
python -m biopytools.janusx gs -i genotypes.vcf.gz -p phenotypes.txt -m GBLUP -o gs_results

# 多模型比较
python -m biopytools.janusx gs -i genotypes.vcf.gz -p phenotypes.txt -m GBLUP rrBLUP BayesA -o gs_results
```

### PCA分析 | PCA Analysis

```bash
# PCA分析
python -m biopytools.janusx pca -i genotypes.vcf.gz -d 5 --plot --plot3d -o pca_results
```

### PostGWAS可视化 | PostGWAS Visualization

```bash
# 生成Manhattan和QQ图
python -m biopytools.janusx postgwas -f gwas_results.lmm.tsv --threshold 1e-6 -o plots
```

## 模块参数说明 | Module Parameters

### GWAS模块 | GWAS Module

| 参数 | 简写 | 类型 | 默认值 | 描述 |
|------|------|------|--------|------|
| `--genotype` | `-i` | str | 必需 | 基因型文件(VCF或PLINK前缀) |
| `--pheno` | `-p` | str | 必需 | 表型文件 |
| `--type` | `-t` | str | vcf | 基因型类型(vcf或bfile) |
| `--models` | `-m` | list | [lmm] | GWAS模型(lm/lmm/fastlmm/farmcpu) |
| `--maf` | | float | 0.02 | 最小等位基因频率阈值 |
| `--geno` | | float | 0.05 | 缺失率阈值 |
| `--grm` | `-k` | str | 1 | 亲缘关系矩阵方法(1/2或路径) |
| `--qcov` | `-q` | int | 0 | PCA数量或Q矩阵文件 |
| `--cov` | `-c` | str | None | 协变量文件 |
| `--ncol` | `-n` | list | None | 表型列索引(零基) |
| `--plot` | | flag | False | 生成图表 |
| `--chunksize` | | int | 100000 | SNP分块大小 |
| `--threads` | `-th` | int | -1 | 线程数(-1为全部核心) |
| `--output-dir` | `-o` | str | ./janusx_gwas_output | 输出目录 |
| `--prefix` | | str | auto | 输出文件前缀 |

### GS模块 | GS Module

| 参数 | 简写 | 类型 | 默认值 | 描述 |
|------|------|------|--------|------|
| `--genotype` | `-i` | str | 必需 | 基因型文件(VCF或PLINK前缀) |
| `--pheno` | `-p` | str | 必需 | 表型文件 |
| `--type` | `-t` | str | vcf | 基因型类型(vcf或bfile) |
| `--models` | `-m` | list | [GBLUP] | GS模型(GBLUP/rrBLUP/BayesA/BayesB/BayesCpi) |
| `--maf` | | float | 0.02 | 最小等位基因频率阈值 |
| `--geno` | | float | 0.05 | 缺失率阈值 |
| `--pcd` | | flag | False | 启用PCA降维 |
| `--ncol` | `-n` | int | None | 表型列索引(零基) |
| `--cv` | | int | None | 交叉验证折数 |
| `--plot` | | flag | False | 生成图表 |
| `--output-dir` | `-o` | str | ./janusx_gs_output | 输出目录 |
| `--prefix` | | str | auto | 输出文件前缀 |

### PCA模块 | PCA Module

| 参数 | 简写 | 类型 | 默认值 | 描述 |
|------|------|------|--------|------|
| `--genotype` | `-i` | str | 必需 | 基因型文件(VCF或PLINK前缀) |
| `--type` | `-t` | str | vcf | 基因型类型(vcf或bfile) |
| `--grm` | | str | None | 预计算的GRM前缀 |
| `--pcfile` | | str | None | 预计算的PCA文件 |
| `--dim` | `-d` | int | 3 | 输出主成分数量 |
| `--plot` | | flag | False | 生成2D散点图 |
| `--plot3d` | | flag | False | 生成3D旋转GIF |
| `--group` | `-g` | str | None | 分组文件路径 |
| `--color` | | int | 1 | 调色板索引(0-6) |
| `--output-dir` | `-o` | str | ./janusx_pca_output | 输出目录 |
| `--prefix` | | str | auto | 输出文件前缀 |

### PostGWAS模块 | PostGWAS Module

| 参数 | 简写 | 类型 | 默认值 | 描述 |
|------|------|------|--------|------|
| `--files` | `-f` | list | 必需 | GWAS结果文件列表 |
| `--chr-col` | | str | #CHROM | 染色体列名 |
| `--pos-col` | | str | POS | 位置列名 |
| `--pvalue-col` | | str | p | P值列名 |
| `--threshold` | | float | None | 显著性阈值 |
| `--noplot` | | flag | False | 禁用绘图 |
| `--color` | | int | 0 | 颜色方案索引(0-6, -1为auto) |
| `--highlight` | `--hl` | str | None | 高亮区域BED文件 |
| `--format` | | str | png | 输出格式(pdf/png/svg/tif) |
| `--anno` | `-a` | str | None | 注释文件路径 |
| `--anno-broaden` | `--ab` | int | None | 注释窗口(kb) |
| `--desc-item` | | str | description | GFF描述键 |
| `--output-dir` | `-o` | str | ./janusx_postgwas_output | 输出目录 |
| `--prefix` | | str | JanusX | 输出文件前缀 |
| `--threads` | `-th` | int | -1 | 线程数(-1为全部核心) |

## 输入文件格式 | Input File Formats

### 表型文件 | Phenotype File

制表符分隔，第一列为样本ID，其余列为表型值：
```
sample1	10.5	0.85
sample2	12.3	0.92
sample3	11.8	0.78
```

### 基因型文件 | Genotype Files

- **VCF格式**: `.vcf` 或 `.vcf.gz`
- **PLINK格式**: `.bed`/`.bim`/`.fam`(使用前缀)

### 分组文件(可选)| Group File (Optional)

制表符分隔：
```
sample1	PopA	Label1
sample2	PopA	Label2
sample3	PopB	Label3
```

## 输出文件 | Output Files

### GWAS输出 | GWAS Output

```
janusx_gwas_output/
├── prefix.trait.lmm.tsv      # LMM分析结果
├── prefix.trait.lmm.svg      # 可视化图表(如果启用--plot)
└── janusx_gwas.log           # 分析日志
```

### GS输出 | GS Output

```
janusx_gs_output/
├── prefix.trait.gs.tsv       # GEBV预测结果
├── prefix.trait.gs.svg       # 预测散点图(如果启用--plot)
└── janusx_gs.log             # 分析日志
```

### PCA输出 | PCA Output

```
janusx_pca_output/
├── prefix.eigenvec           # PC坐标
├── prefix.eigenvec.id        # 样本ID
├── prefix.eigenval           # 特征值
├── prefix.eigenvec.2D.pdf    # 2D散点图(如果启用--plot)
├── prefix.eigenvec.3D.gif    # 3D旋转GIF(如果启用--plot3d)
└── janusx_pca.log            # 分析日志
```

### PostGWAS输出 | PostGWAS Output

```
janusx_postgwas_output/
├── file.manh.png             # Manhattan图
├── file.qq.png               # QQ图
├── threshold.anno.tsv        # 显著SNP注释(如果指定--threshold)
└── janusx_postgwas.log       # 分析日志
```

## 使用示例 | Usage Examples

### 示例1：完整GWAS流程 | Example 1: Complete GWAS Pipeline

```bash
# 1. PCA分析(计算群体结构)
python -m biopytools.janusx pca -i genotypes.vcf.gz -d 5 -o pca_results

# 2. GWAS分析(使用PCA作为协变量)
python -m biopytools.janusx gwas \
    -i genotypes.vcf.gz \
    -p phenotypes.txt \
    -m lmm \
    -q 3 \
    -o gwas_results

# 3. 可视化结果
python -m biopytools.janusx postgwas \
    -f gwas_results/*.lmm.tsv \
    --threshold 1e-6 \
    -o plots
```

### 示例2：基因组选择 | Example 2: Genomic Selection

```bash
# 使用GBLUP和rrBLUP进行预测
python -m biopytools.janusx gs \
    -i genotypes.vcf.gz \
    -p phenotypes.txt \
    -m GBLUP rrBLUP \
    --cv 5 \
    --plot \
    -o gs_results
```

### 示例3：使用PLINK格式 | Example 3: Using PLINK Format

```bash
# GWAS with PLINK format
python -m biopytools.janusx gwas \
    -i /path/to/plink_prefix \
    -t bfile \
    -p phenotypes.txt \
    -m lmm \
    -o gwas_results
```

### 示例4：多模型比较 | Example 4: Multi-Model Comparison

```bash
# 比较不同GWAS模型
python -m biopytools.janusx gwas \
    -i genotypes.vcf.gz \
    -p phenotypes.txt \
    -m lm lmm fastlmm farmcpu \
    -o gwas_comparison
```

## 系统要求 | System Requirements

### 依赖软件 | Dependencies

- **JanusX** v1.0.10+ (预编译二进制文件)
  - 默认路径: `/share/org/YZWL/yzwl_lixg/software/JanusX/JanusX.bin`
- **Python** 3.10+
- **Python包**:
  - click
  - pathlib
  - logging
  - subprocess

### 硬件建议 | Hardware Recommendations

- **CPU**: 多核处理器(推荐4核以上)
- **RAM**: 最少8GB(大基因组推荐32GB以上)
- **存储**: 预留基因型文件大小3倍的磁盘空间

## 注意事项 | Important Notes

1. **内存使用**: FarmCPU模型需要将全部基因型数据加载到内存，大样本量请谨慎使用
2. **模型选择**: 无群体结构时使用GLM，有群体结构或亲缘关系时使用LMM
3. **缓存文件**: GRM和PCA会自动缓存，可加速后续分析
4. **文件路径**: 确保基因型和表型文件的样本ID一致

## 常见问题 | FAQ

### Q: 如何选择GWAS模型？
**A**: 大样本无群体结构→GLM；有群体结构→LMM；快速筛选→fastLMM；高功效需全内存→FarmCPU

### Q: 如何使用预计算的GRM？
**A**: 使用`-k`参数指定GRM文件路径(不含扩展名)

### Q: 如何解读QQ图？
**A**: 点沿对角线→正常；右上偏离→显著关联；整体偏离→群体结构或技术混杂

### Q: PCA分析的输入源优先级？
**A**: grm > pcfile > genotype，按此顺序自动选择

## 引用信息 | Citation

```bibtex
@software{JanusX,
  title = {JanusX: High-Performance GWAS and Genomic Selection Suite},
  author = {Jingxian FU},
  url = {https://github.com/FJingxian/JanusX},
  year = {2025}
}
```

## 相关资源 | Related Resources

- [JanusX GitHub](https://github.com/FJingxian/JanusX)
- [GEMMA软件](https://github.com/genetics-statistics/GEMMA)
- [GCTA软件](https://yanglab.westlake.edu.cn/software/gcta)
- [rMVP软件](https://github.com/xiaolei-lab/rMVP)
