# 🧬 GEMMA GWAS 全基因组关联分析模块

**基于GEMMA的高效混合线性模型GWAS分析工具 | Efficient LMM-based GWAS Analysis Tool**

## 📖 功能概述 | Overview

GEMMA GWAS模块是基于GEMMA软件的全基因组关联分析工具，使用线性混合模型（LMM）控制群体结构和亲缘关系，适用于各种动植物GWAS研究。支持多表型批量分析、自动PCA主成分分析、Kinship矩阵计算，提供高质量的关联分析结果。

## ✨ 主要特性 | Key Features

- **🎯 混合线性模型**: 基于GEMMA的LMM方法，有效控制群体结构和亲缘关系
- **📋 多表型批量分析**: 自动识别表型文件中的所有表型并批量处理
- **🧮 PCA主成分分析**: 自动计算PCA主成分作为群体结构协变量
- **🔬 Kinship矩阵**: 自动计算亲缘关系矩阵控制样本相关性
- **⚡ 三种检验方法**: 支持Wald检验、似然比检验、Score检验或全部三种
- **📊 灵活质控**: 支持PLINK和GEMMA两套质控参数
- **🧵 高性能并行**: 支持多线程加速分析
- **📈 详细输出**: 生成标准关联分析结果和统计汇总

## 🚀 快速开始 | Quick Start

### 基本用法 | Basic Usage

```bash
# 基本GWAS分析 - 自动处理所有表型
biopytools gemma-gwas \
    -i genotype.vcf.gz \
    -p phenotypes.txt \
    -o results

# 指定PCA主成分数量
biopytools gemma-gwas \
    -i genotype.vcf.gz \
    -p phenotypes.txt \
    -o results \
    --n-pca 10

# 跳过PLINK质控（如果数据已过滤）
biopytools gemma-gwas \
    -i genotype.vcf.gz \
    -p phenotypes.txt \
    -o results \
    --no-qc

# 自定义GEMMA质控参数
biopytools gemma-gwas \
    -i genotype.vcf.gz \
    -p phenotypes.txt \
    -o results \
    --maf-gemma 0.01 \
    --miss-gemma 0.05
```

### 高级用法 | Advanced Usage

```bash
# 使用不同的LMM检验方法
biopytools gemma-gwas \
    -i genotype.vcf.gz \
    -p phenotypes.txt \
    -o results \
    --lmm 1  # 1=Wald, 2=LRT, 3=Score, 4=all

# 使用standardized kinship矩阵
biopytools gemma-gwas \
    -i genotype.vcf.gz \
    -p phenotypes.txt \
    -o results \
    --gk 2

# 跳过SNP输出（加快速度）
biopytools gemma-gwas \
    -i genotype.vcf.gz \
    -p phenotypes.txt \
    -o results \
    --notsnp

# 完整参数示例
biopytools gemma-gwas \
    -i genotype.vcf.gz \
    -p phenotypes.txt \
    -o results \
    --n-pca 15 \
    --maf 0.05 \
    --geno 0.05 \
    --lmm 4 \
    --threads 24 \
    --gemma /path/to/gemma
```

## 📋 参数说明 | Parameters

### 必需参数 | Required Parameters

| 参数 | 描述 | 示例 |
|------|------|------|
| `-i, --input` | VCF基因型文件路径（支持压缩） | `-i genotype.vcf.gz` |
| `-p, --pheno` | 表型文件路径（第一列样本ID，有表头） | `-p traits.txt` |
| `-o, --outdir` | 输出目录路径 | `-o gwas_results` |

### 分析参数 | Analysis Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--n-pca` | `10` | 🧮 PCA主成分数量（作为协变量） |
| `--threads` | `12` | 🧵 PLINK并行线程数 |

### GEMMA参数 | GEMMA Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--lmm` | `4` | 🧬 LMM检验方法 (1=Wald, 2=LRT, 3=Score, 4=all) |
| `--gk` | `1` | 🔬 Kinship矩阵方法 (1=centered, 2=standardized) |
| `--miss-gemma` | `1.0` | 💧 GEMMA缺失率阈值 (1.0=不过滤) |
| `--maf-gemma` | `0.0` | 🧬 GEMMA MAF阈值 (0.0=不过滤) |
| `--notsnp` | `False` | ⏭️ 不输出每个SNP的估计值（加快速度） |
| `--gemma` | `auto` | 🔧 GEMMA程序路径 |

### PLINK质控参数 | PLINK QC Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--maf` | `0.05` | 🧬 最小等位基因频率 |
| `--geno` | `0.1` | 📊 SNP最大缺失率 |
| `--mind` | `0.1` | 👥 样本最大缺失率 |
| `--hwe` | `1e-6` | 🧪 Hardy-Weinberg平衡p值阈值 |
| `--no-qc` | `False` | ⏭️ 跳过PLINK质控 |

### 日志参数 | Logging Parameters

| 参数 | 描述 |
|------|------|
| `-v, --verbose` | 详细模式 (-v: INFO, -vv: DEBUG) |
| `--quiet` | 静默模式 (仅ERROR) |
| `--log-file` | 日志文件路径 |

## 📁 输入文件格式 | Input File Formats

### 1. VCF基因型文件 | VCF Genotype File

标准的VCF格式文件，支持gzip压缩（`.vcf.gz`）：

```vcf
##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
#CHROM  POS     ID        REF  ALT  QUAL  FILTER  INFO  FORMAT  Sample1  Sample2
chr1    100     rs001     A    G    .     .      .     GT     0/0      0/1
chr1    200     rs002     C    T    .     .      .     GT     0/1      1/1
```

**要求**：
- 必须包含GT（基因型）字段
- 可以是压缩格式（`.vcf.gz`）
- 自动处理非标准染色体命名（使用`--allow-extra-chr`）

### 2. 表型文件 | Phenotype File

纯文本格式，第一列为样本ID，第一行为表头：

```text
id        trait1   trait2   trait3
Sample1   1.23     4.56     7.89
Sample2   2.34     5.67     8.90
Sample3   3.45     6.78     9.01
```

**要求**：
- 第一列：样本ID（必须与VCF中的样本ID一致）
- 第一行：表头（表型名称）
- 分隔符：制表符（`\t`）
- 表型值：数值型（连续性状）
- 样本顺序可以与VCF不同

## 📂 输出文件结构 | Output File Structure

```
gemma_results/
├── output/                          # GEMMA分析结果目录
│   ├── kinship.cXX.txt             # 亲缘关系矩阵
│   ├── acdA_lmm.assoc.txt          # acdA表型的GWAS结果
│   ├── acdA_lmm.log.txt            # acdA表型的GEMMA日志
│   ├── acdB_lmm.assoc.txt          # acdB表型的GWAS结果
│   ├── ...
│   └── traitN_lmm.assoc.txt        # 最后一个表型的结果
│
├── genotype.*                       # PLINK格式基因型文件
│   ├── genotype.bed
│   ├── genotype.bim
│   └── genotype.fam
│
├── pca.*                            # PCA分析结果
│   ├── pca.eigenvec                # 主成分特征向量
│   └── pca.eigenval                # 特征值
│
├── covariate.txt                    # 协变量文件（PCA主成分）
├── 16.微生物和甲烷菌用于GWAS的数据_no_header.txt  # 表型文件（无表头）
│
├── plink_convert.log                # PLINK转换日志
├── plink_pca.log                    # PLINK PCA日志
├── gemma_kinship.log                # GEMMA kinship日志
├── gemma_acdA.log                   # 每个表型的GEMMA日志
├── gemma_acdB.log
├── ...
│
└── gwas_summary.txt                 # 分析汇总报告
```

## 📊 分析流程 | Analysis Workflow

```
┌─────────────────────────────────────────────────────────────────┐
│                    GEMMA GWAS 分析流程                          │
├─────────────────────────────────────────────────────────────────┤
│                                                                   │
│  1. PLINK格式转换                                                │
│     VCF → PLINK (BED/BIM/FAM)                                    │
│                                                                   │
│  2. 质量控制（可选）                                               │
│     MAF过滤、缺失率过滤、HWE检验                                  │
│                                                                   │
│  3. FAM文件表型列更新                                            │
│     从表型文件提取第一个表型值更新FAM文件                         │
│                                                                   │
│  4. 表型文件准备                                                  │
│     生成带表头和不带表头两个版本                                  │
│                                                                   │
│  5. PCA主成分分析                                                 │
│     PLINK计算前n个主成分                                          │
│     → pca.eigenvec (主成分)                                       │
│     → pca.eigenval (特征值)                                       │
│     → covariate.txt (协变量文件)                                 │
│                                                                   │
│  6. Kinship矩阵计算                                              │
│     GEMMA计算亲缘关系矩阵                                         │
│     → output/kinship.cXX.txt                                     │
│                                                                   │
│  7. GWAS关联分析（每个表型独立分析）                             │
│     GEMMA LMM分析                                                 │
│     使用K矩阵 + PCA主成分作为协变量                               │
│     → output/{trait}_lmm.assoc.txt (GWAS结果)                    │
│     → output/{trait}_lmm.log.txt (GEMMA日志)                     │
│                                                                   │
│  8. 汇总报告生成                                                  │
│     → gwas_summary.txt                                           │
│                                                                   │
└─────────────────────────────────────────────────────────────────┘
```

## 📈 输出文件格式 | Output File Formats

### 1. GWAS结果文件 (*.assoc.txt)

```
chr     rs      ps      nmiss   allele1 allele2 af       beta     se      p_wald  p_lrt   p_score
chr1    none    100     87      A       G       0.123    -0.045   0.023   0.052   0.048   0.051
chr1    none    200     87      C       T       0.089    0.067    0.019   1.2e-05  1.1e-05  1.3e-05
```

**列说明**：
- `chr`: 染色体
- `rs`: SNP ID（通常为none）
- `ps`: 物理位置
- `nmiss`: 非缺失样本数
- `allele1/allele2`: 等位基因
- `af`: 等位基因频率
- `beta`: 效应大小（回归系数）
- `se`: 标准误
- `p_wald`: Wald检验p值
- `p_lrt`: 似然比检验p值
- `p_score`: Score检验p值

### 2. 汇总报告 (gwas_summary.txt)

```text
Phenotype_Name,Total_SNPs,Significant_1e-5,Significant_1e-6,Significant_1e-7
acdA,15103,125,87,52
acdB,15103,203,145,98
```

## 🔬 原理说明 | Methodology

### 线性混合模型（LMM）

GEMMA使用线性混合模型控制群体结构和亲缘关系：

**模型公式**：
```
y = Xβ + Zu + ε
```

其中：
- `y`: 表型值向量
- `X`: 固定效应矩阵（SNP基因型 + PCA主成分）
- `β`: 固定效应系数
- `Z`: 亲缘关系设计矩阵
- `u`: 随机效应（u ~ N(0, σg²K)），K为kinship矩阵
- `ε`: 残差（ε ~ N(0, σe²I)）

### PCA主成分的作用

PCA主成分作为协变量控制群体分层（population stratification）：
- 第1-10个主成分通常能解释主要的群体结构
- 减少假阳性关联

### Kinship矩阵的作用

Kinship矩阵控制样本间的亲缘关系/相关性：
- `gk=1`: centered相关矩阵（默认）
- `gk=2`: standardized相关矩阵

## ⚠️ 注意事项 | Important Notes

### 1. 质控策略

**默认情况下**：
- **PLINK质控**：关闭（`--no-qc`需要显式指定）
- **GEMMA质控**：基本关闭（maf=0.0, miss=1.0）

**原因**：假设用户在输入前已经完成质控。如需启用，请：
```bash
# 启用PLINK质控
biopytools gemma-gwas -i input.vcf.gz -p pheno.txt -o results \
    --maf 0.05 --geno 0.1 --mind 0.1

# 启用GEMMA质控
biopytools gemma-gwas -i input.vcf.gz -p pheno.txt -o results \
    --maf-gemma 0.01 --miss-gemma 0.05
```

### 2. 样本ID一致性

VCF和表型文件中的样本ID必须完全一致（包括大小写）。工具会自动同步样本。

### 3. 内存需求

GEMMA分析需要大量内存，样本数较多时建议：
- 使用高性能服务器
- 考虑分批处理表型
- 使用`--notsnp`参数减少内存占用

### 4. 分析时间

- **Kinship矩阵计算**：约5-15分钟（取决于样本数和SNP数）
- **每个表型的GWAS**：约3-10分钟
- 建议：2201个表型可能需要数小时到数天

### 5. 错误处理

- 每个表型独立分析，一个表型失败不影响其他表型
- 检查各个`gemma_{trait}.log`文件了解失败原因
- 常见错误：
  - 样本数不匹配
  - 表型值全为缺失
  - 内存不足

## 📚 参考文献 | References

1. **GEMMA软件**: Zhou X, Stephens M. (2012) Genome-wide efficient mixed-model analysis for association studies. *Nature Genetics* 44: 821-824.

2. **线性混合模型**: Zhou X, Stephens M. (2014) Efficient multivariate linear mixed model algorithms for GWAS. *Nature Methods* 11: 1144-1146.

3. **GEMMA文档**: http://www.xzlab.org/software.html

## 💡 使用建议 | Usage Tips

### 1. 分析前准备

```bash
# 检查样本ID一致性
bcftools query -l input.vcf.gz | sort > vcf_samples.txt
cut -f1 phenotypes.txt | tail -n +2 | sort > pheno_samples.txt
diff vcf_samples.txt pheno_samples.txt
```

### 2. 小规模测试

```bash
# 先用少量表型测试
head -n 3 phenotypes.txt > test_phenotypes.txt
biopytools gemma-gwas -i input.vcf.gz -p test_phenotypes.txt -o test_results
```

### 3. 查看分析进度

```bash
# 查看当前运行的表型
tail -f gwas_results/*.log

# 统计已完成的表型数
ls gwas_results/output/*_lmm.assoc.txt | wc -l
```

### 4. 结果可视化

使用GWAS结果文件绘制曼哈顿图和QQ图（需其他工具）。

## 🆘 常见问题 | FAQ

**Q: 为什么默认GEMMA质控是关闭的？**

A: 因为GEMMA的默认参数比较严格（maf=0.01, miss=0.05），会过滤掉很多SNP。如果您在输入前已经做过质控，建议保持默认（maf=0.0, miss=1.0）。

**Q: PCA主成分数量如何选择？**

A: 通常使用前3-10个主成分。可以使用肘部图（scree plot）或解释方差比例来决定。

**Q: Kinship矩阵的gk方法选哪个？**

A: 两者都可以，`gk=1`（centered）是默认和最常用的选择。`gk=2`（standardized）在某些情况下可能更合适。

**Q: 如何加快分析速度？**

A:
1. 使用`--notsnp`跳过SNP输出
2. 减少PCA主成分数量
3. 使用`--lmm 1`只运行Wald检验（最快）
4. 增加线程数（`--threads`）

**Q: 一个表型失败会影响其他表型吗？**

A: 不会。每个表型独立分析，失败只影响该表型。

## 📞 技术支持 | Support

如有问题或建议，请联系：BioTools Development Team
