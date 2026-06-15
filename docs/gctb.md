# 🧬 GCTB 全基因组复杂性状贝叶斯分析模块

**专业的基因组复杂性状贝叶斯分析工具 | Professional Genome-wide Complex Trait Bayesian Analysis Tool**

## 📖 功能概述 | Overview

GCTB 分析模块是一个强大的全基因组复杂性状贝叶斯分析工具，基于 GCTB 软件构建，提供从 VCF 文件到最终贝叶斯分析结果的完整自动化流程。支持个体水平和汇总水平数据分析，实现了数据转换、质量控制、LD 矩阵构建和贝叶斯分析的全流程自动化，适用于各种复杂性状的遗传架构分析研究。

## ✨ 主要特性 | Key Features

- **🔄 一键式分析流程**: VCF 转换 → 质控 → LD 矩阵 → 贝叶斯分析全自动化
- **🎯 灵活的分析模式**: 支持个体水平和汇总水平两种分析模式
- **🛡️ 智能质量控制**: MAF、缺失率、HWE 检验等质控参数可配置
- **📊 多种 LD 矩阵**: 支持稀疏、分块、特征值分解三种 LD 矩阵类型
- **🧮 三种贝叶斯模型**: 支持 BayesS、BayesR、BayesC 分析
- **⚙️ 高度可配置**: 丰富的参数设置，满足不同研究需求
- **📝 详细日志记录**: 完整的分析过程日志和错误追踪
- **💾 断点续传**: 支持从任意步骤继续分析

## 🚀 快速开始 | Quick Start

### 基本用法 | Basic Usage

```bash
# 个体水平分析
biopytools gctb \
    -i variants.vcf \
    -p phenotype.txt \
    -o results

# 汇总水平分析
biopytools gctb \
    -i variants.vcf \
    -p summary_stats.ma \
    --analysis-mode summary \
    -o results
```

### 高级用法 | Advanced Usage

```bash
# 使用 BayesR 模型进行汇总水平分析
biopytools gctb \
    -i variants.vcf \
    -p summary_stats.ma \
    --analysis-mode summary \
    --bayes-type R \
    --ld-matrix-type eigen \
    --maf-threshold 0.05 \
    -o results
```

## 📋 参数说明 | Parameters

### 必需参数 | Required Parameters

| 参数 | 描述 | 示例 |
|------|------|------|
| `-i, --vcf-file` | VCF 变异文件路径 | `variants.vcf` |
| `-p, --pheno-file` | 表型文件或 GWAS 汇总统计文件 | `phenotype.txt` |

### 输出配置 | Output Configuration

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-o, --output-dir` | `./gctb_output` | 📁 输出目录路径 |

### 软件配置 | Software Configuration

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--gctb-path` | `~/miniforge3/envs/gctb/bin/gctb` | 🛠️ GCTB 软件路径 |
| `--plink-path` | `~/miniforge3/envs/Population_genetics/bin/plink` | 🛠️ PLINK 软件路径 |

### 质量控制配置 | Quality Control Configuration

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--maf-threshold` | `0.01` | 🎯 MAF 过滤阈值 |
| `--miss-threshold` | `0.1` | 🎯 缺失率过滤阈值 |
| `--hwe-p` | `1e-6` | 🎯 HWE p 值阈值 |

### 分析配置 | Analysis Configuration

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--bayes-type` | `S` | 🧮 贝叶斯模型类型 (S/R/C) |
| `--analysis-mode` | `individual` | 📊 分析模式 (individual/summary) |
| `--ld-matrix-type` | `sparse` | 📊 LD 矩阵类型 (sparse/block/eigen) |

### 高级参数 | Advanced Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--threads` | `24` | 🔧 线程数 |
| `--seed` | `None` | 🎲 随机种子 |
| `--pi` | `None` | 📊 Polygenicity 参数 |
| `--sigma-g` | `None` | 📊 遗传方差 |
| `--rho` | `None` | 📊 SNP 效应与 MAF 关系参数 |

### 步骤控制 | Step Control

| 参数 | 描述 |
|------|------|
| `--step` | 🎯 只运行指定步骤 (convert/qc/freq/ld/analysis) |

## 📁 输入文件格式 | Input File Formats

### VCF 文件 | VCF File

标准 VCF 格式的基因型数据文件：

```vcf
##fileformat=VCFv4.2
##contig=<ID=chr1,length=248956422>
#CHROM  POS    ID     REF  ALT  QUAL  FILTER  INFO  FORMAT  Sample1
chr1    10177  rs1    A    AC   100   PASS    .     GT      0/1
chr1    10352  rs2    T    TA   100   PASS    .     GT      1/1
```

### 表型文件 | Phenotype File

个体水平分析使用的表型文件：

```
FID  IID  PHENO
1    1    1.234
1    2    -0.567
```

### GWAS 汇总统计文件 | GWAS Summary Statistics File

汇总水平分析使用的 GWAS 汇总统计文件：

```
SNP    A1   A2   FREQ   BETA    SE        P        N
rs1    A    C    0.25   0.123   0.045    1.2e-5  1000
rs2    T    G    0.30   -0.089  0.038    0.002   1000
```

## 💡 使用示例 | Usage Examples

### 示例1：基本个体水平分析 | Example 1: Basic Individual-level Analysis

```bash
biopytools gctb \
    -i /path/to/genotypes.vcf \
    -p /path/to/phenotypes.txt \
    -o individual_analysis
```

### 示例2：汇总水平分析 | Example 2: Summary-level Analysis

```bash
biopytools gctb \
    -i /path/to/genotypes.vcf \
    -p /path/to/gwas_summary.ma \
    --analysis-mode summary \
    --ld-matrix-type sparse \
    -o summary_analysis
```

### 示例3：使用 BayesR 模型 | Example 3: Using BayesR Model

```bash
biopytools gctb \
    -i genotypes.vcf \
    -p phenotypes.txt \
    --bayes-type R \
    --maf-threshold 0.05 \
    -o bayesr_analysis
```

### 示例4：分步执行 | Example 4: Step-by-Step Execution

```bash
# 步骤1: VCF 转换
biopytools gctb -i genotypes.vcf -p phenotypes.txt --step convert -o results

# 步骤2: 质量控制
biopytools gctb -i genotypes.vcf -p phenotypes.txt --step qc -o results

# 步骤3: 计算频率
biopytools gctb -i genotypes.vcf -p phenotypes.txt --step freq -o results

# 步骤4: 构建 LD 矩阵
biopytools gctb -i genotypes.vcf -p phenotypes.txt --step ld -o results

# 步骤5: 运行分析
biopytools gctb -i genotypes.vcf -p phenotypes.txt --step analysis -o results
```

### 示例5：使用特征值分解 LD 矩阵 | Example 5: Using Eigen-decomposed LD Matrix

```bash
biopytools gctb \
    -i genotypes.vcf \
    -p summary_stats.ma \
    --analysis-mode summary \
    --ld-matrix-type eigen \
    --bayes-type S \
    -o eigen_ld_analysis
```

## 📊 输出结果 | Output Results

### 输出目录结构 | Output Directory Structure

```
gctb_output/
├── 01_vcf_to_plink/           # VCF 转换结果
│   ├── input.bed
│   ├── input.bim
│   └── input.fam
├── 02_quality_control/        # 质量控制结果
│   ├── input_qc.bed
│   ├── input_qc.bim
│   ├── input_qc.fam
│   └── input_qc.frq
├── 03_ld_matrix/              # LD 矩阵
│   └── ld_matrix.ldm.sparse
├── 04_gctb_analysis/          # GCTB 分析结果
│   ├── bayess.snpRes
│   ├── bayess.parRes
│   └── bayess.log
└── 99_logs/                   # 日志文件
    └── gctb_pipeline.log
```

### 关键输出文件说明 | Key Output Files Description

#### SNP 结果文件 | SNP Results File

`bayesS.snpRes` 或 `sbayesS.snpRes`：每个 SNP 的分析结果

| 列名 | 描述 |
|------|------|
| SNP | SNP ID |
| CHR | 染色体编号 |
| BP | 物理位置 |
| A1 | 等位基因1 |
| A2 | 等位基因2 |
| BETA | 效应大小估计 |
| SE | 标准误 |
| PIP | 后验包含概率 |
| AVG | 后验均值 |

#### 参数结果文件 | Parameter Results File

`bayesS.parRes` 或 `sbayesS.parRes`：模型参数估计结果

| 参数 | 描述 |
|------|------|
| PVE | 解释的遗传方差比例 |
| PG | 多基因性 |
| H2 | 遗传力 |
| RHO | 效应大小与 MAF 关系参数 |

## 🔧 系统要求 | System Requirements

### 依赖软件 | Dependencies

- **GCTB** (版本 2.5.5)
  - 安装在 gctb conda 环境中
- **PLINK2**
  - 用于 VCF 到 PLINK 格式转换
  - 安装在 Population_genetics conda 环境中
- **OpenMP**
  - 用于并行计算
- **Python** (版本 3.7+)
- **Python 包**:
  - `click` - 命令行界面
  - `pathlib` - 路径处理

## ⚠️ 注意事项 | Important Notes

1. **VCF 文件格式**: 确保 VCF 文件格式正确，包含必需的头部信息
2. **表型文件格式**: 表型文件必须包含 FID、IID 和表型值列
3. **内存使用**: 构建完整 LD 矩阵需要大量内存，建议使用稀疏或分块 LD 矩阵
4. **线程数**: 根据可用 CPU 核心数合理设置线程数
5. **质控参数**: 根据研究需求调整质控参数阈值

## 🐛 故障排除 | Troubleshooting

### 常见问题 | Common Issues

**Q: "VCF 文件读取失败" 错误**
```bash
# 检查 VCF 文件格式
bcftools view genotypes.vcf | head -20

# 确保文件存在且格式正确
```

**Q: "质控后样本数为 0"**
```bash
# 放宽质控阈值
biopytools gctb ... --maf-threshold 0.001 --miss-threshold 0.2
```

**Q: "内存不足" 错误**
```bash
# 使用稀疏 LD 矩阵而不是分块 LD 矩阵
biopytools gctb ... --ld-matrix-type sparse

# 或减少线程数
biopytools gctb ... --threads 8
```

**Q: "GCTB 运行失败"**
```bash
# 检查 GCTB 是否可用
conda activate gctb
gctb --version

# 查看日志文件
cat gctb_output/99_logs/gctb_pipeline.log
```

## 📚 相关资源 | Related Resources

- [GCTB 官方文档](https://gctbhub.cloud.edu.au/software/gctb/)
- [GCTB GitHub 仓库](https://github.com/jianzeng/GCTB)
- [GCTB 论文](https://www.nature.com/articles/s41588-018-0101-4)

## 📄 许可证 | License

本项目采用 MIT 许可证 - 详见 [LICENSE](LICENSE) 文件

## 🔬 引用信息 | Citation

如果在学术研究中使用此工具，请引用 GCTB 相关文献：

```
Zeng, J. et al. (2018).
Statistical evidence for the importance of the genetic architecture in the prediction of complex traits.
Nature Genetics, 50, 743–749.
doi:10.1038/s41588-018-0101-4
```
