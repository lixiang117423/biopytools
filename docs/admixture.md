# 🧬 ADMIXTURE 群体结构分析模块

**高效的群体遗传结构分析工具 | Efficient Population Genetic Structure Analysis Tool**

## 📖 功能概述 | Overview

ADMIXTURE 群体结构分析模块是一个强大的群体遗传学工具，基于最大似然估计算法进行祖先成分分析，支持自动K值优化、交叉验证和并行处理，适用于各种规模的群体遗传学研究和基因组关联分析预处理。

## ✨ 主要特性 | Key Features

- **🎚️ 智能K值分析**: 自动运行指定K值范围，寻找最优群体数
- **🔄 交叉验证优化**: 内置CV折数验证，科学确定最佳K值
- **⚡ 高性能并行**: 多线程加速计算，最大化利用计算资源
- **🛡️ 全面质量控制**: MAF、缺失率、HWE平衡多维过滤
- **📊 详细统计报告**: 完整的分析日志和结果汇总
- **🔧 灵活配置**: 支持跳过预处理，保留中间文件等选项

## 🚀 快速开始 | Quick Start

### 基本用法 | Basic Usage

```bash
# 基本ADMIXTURE分析
biopytools admixture -v input.vcf -o admixture_results

# 指定K值范围
biopytools admixture -v data.vcf -o results -k 2 -K 10

# 高性能分析
biopytools admixture -v large_dataset.vcf -o results -t 16
```

### 高级用法 | Advanced Usage

```bash
# 严格质控的完整分析
biopytools admixture -v population.vcf -o strict_analysis \
    --maf 0.05 --missing 0.05 --hwe 1e-5 \
    -k 2 -K 12 -c 10 -t 24

# 跳过预处理的快速分析
biopytools admixture -v clean_data.vcf -o quick_results \
    --skip-preprocessing --keep-intermediate \
    -k 3 -K 8 -t 32
```

## 📋 参数说明 | Parameters

### 必需参数 | Required Parameters

| 参数 | 描述 | 示例 |
|------|------|------|
| `-v, --vcf` | VCF基因型文件路径 | `-v population.vcf` |

### 输出配置 | Output Configuration

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-o, --output` | `admixture_results` | 📁 输出目录路径 |
| `-i, --keep-intermediate` | `False` | 💾 保留中间处理文件 |

### K值分析配置 | K-value Analysis Configuration

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-k, --min-k` | `2` | 📉 最小K值（最少祖先群体数） |
| `-K, --max-k` | `10` | 📈 最大K值（最多祖先群体数） |
| `-c, --cv-folds` | `5` | 🔄 交叉验证折数 |

### 质量控制参数 | Quality Control Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-m, --maf` | `0.01` | 📊 最小等位基因频率阈值 |
| `-M, --missing` | `0.1` | 🗑️ 最大缺失率阈值 |
| `-H, --hwe` | `1e-6` | ⚖️ Hardy-Weinberg平衡p值阈值 |
| `-s, --skip-preprocessing` | `False` | ⏭️ 跳过VCF预处理步骤 |

### 性能配置 | Performance Configuration

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-t, --threads` | `4` | 🧵 并行线程数 |

## 📁 输入文件格式 | Input File Formats

### VCF文件要求 | VCF File Requirements

支持标准VCF格式文件（压缩或未压缩）：

```vcf
##fileformat=VCFv4.2
##contig=<ID=1,length=249250621>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample1	Sample2	Sample3
1	14370	rs6054257	G	A	29	PASS	.	GT	0/0	1/0	1/1
1	17330	.	T	A	3	q10	.	GT	0/0	0/1	0/0
1	1110696	rs6040355	A	G,T	67	PASS	.	GT	1/2	2/1	2/2
```

**支持的文件格式**:
- `.vcf` - 标准VCF格式
- `.vcf.gz` - gzip压缩VCF格式

**文件要求**:
- 包含完整的基因型信息（GT字段）
- 至少包含2个样本
- 建议预先进行基本质量过滤

## 💡 使用示例 | Usage Examples

### 示例1：基础群体结构分析 | Example 1: Basic Population Structure Analysis

```bash
# 探索2-8个祖先群体的结构
biopytools admixture \
    -v population_samples.vcf \
    -o basic_structure_analysis \
    -k 2 -K 8 \
    -t 8
```

### 示例2：高质量严格分析 | Example 2: High-Quality Strict Analysis

```bash
# 严格的质量控制和详细的K值搜索
biopytools admixture \
    -v large_cohort.vcf.gz \
    -o high_quality_analysis \
    --maf 0.05 \
    --missing 0.02 \
    --hwe 1e-5 \
    -k 2 -K 15 \
    -c 10 \
    -t 24
```

### 示例3：快速分析预处理数据 | Example 3: Fast Analysis of Pre-processed Data

```bash
# 跳过质控，直接分析已清理的数据
biopytools admixture \
    -v clean_genotypes.vcf \
    -o fast_analysis \
    --skip-preprocessing \
    --keep-intermediate \
    -k 3 -K 10 \
    -t 16
```

### 示例4：精确交叉验证分析 | Example 4: Precise Cross-Validation Analysis

```bash
# 使用更多CV折数提高K值选择精度
biopytools admixture \
    -v diverse_population.vcf \
    -o precise_cv_analysis \
    -k 2 -K 12 \
    -c 20 \
    --maf 0.02 \
    --missing 0.05 \
    -t 32
```

### 示例5：大规模数据集分析 | Example 5: Large-Scale Dataset Analysis

```bash
# 处理大规模基因组数据
biopytools admixture \
    -v genome_wide_snps.vcf.gz \
    -o large_scale_admixture \
    -k 2 -K 20 \
    -c 15 \
    --maf 0.01 \
    --missing 0.1 \
    --hwe 1e-6 \
    -t 64 \
    --keep-intermediate
```

### 关键输出文件说明 | Key Output Files Description

- **\*.Q文件**: 个体在各祖先群体中的成分比例
- **\*.P文件**: 各祖先群体的等位基因频率
- **cv_errors.txt**: 所有K值对应的交叉验证误差
- **optimal_k_report.txt**: 推荐的最优K值及其统计支持
- **ancestry_proportions.txt**: 所有个体的祖先成分表格

## 🔧 系统要求 | System Requirements

### 依赖软件 | Dependencies

- **ADMIXTURE** (版本 1.3.0+)
  - 下载地址: https://dalexander.github.io/admixture/
- **PLINK** (版本 1.9+)
  - 用于VCF预处理和格式转换
- **Python** (版本 3.7+)
- **Python包**:
  - `pandas` - 数据处理
  - `numpy` - 数值计算

### 安装依赖软件 | Installing Dependencies

```bash
# 安装ADMIXTURE
wget https://dalexander.github.io/admixture/binaries/admixture_linux-1.3.0.tar.gz
tar -xzf admixture_linux-1.3.0.tar.gz
sudo mv admixture_linux-1.3.0/admixture /usr/local/bin/

# 安装PLINK
wget https://www.cog-genomics.org/static/bin/plink190904/plink_linux_x86_64.zip
unzip plink_linux_x86_64.zip
sudo mv plink /usr/local/bin/

# 安装Python包
pip install pandas numpy matplotlib click
```

### 硬件建议 | Hardware Recommendations

- **CPU**: 多核处理器（推荐16核以上用于大数据集）
- **RAM**: 最少8GB（大数据集推荐32GB以上）
- **存储**: 至少预留数据集大小3倍的磁盘空间
- **网络**: 如需下载参考数据，建议高速网络连接

## ⚠️ 注意事项 | Important Notes

1. **数据质量**: 输入VCF文件质量直接影响ADMIXTURE分析结果
2. **K值选择**: CV误差最小的K值通常是最优选择，但需结合生物学意义
3. **收敛性**: 某些K值可能需要多次运行以确保结果收敛
4. **计算时间**: 大数据集和高K值会显著增加计算时间
5. **内存使用**: 高密度SNP数据可能需要大量内存

## 🐛 故障排除 | Troubleshooting

### 常见问题 | Common Issues

**Q: "admixture: command not found" 错误**
```bash
# 检查ADMIXTURE安装
which admixture
# 如未安装，请按上述方法安装
```

**Q: CV误差未显示明显最小值**
```bash
# 增加K值搜索范围
biopytools admixture -v input.vcf -o results -k 2 -K 15

# 或调整质控参数
biopytools admixture -v input.vcf -o results --maf 0.05
```

**Q: 内存不足错误**
```bash
# 减少线程数
biopytools admixture -v input.vcf -o results -t 4

# 或增加MAF阈值减少SNP数量
biopytools admixture -v input.vcf -o results --maf 0.05
```

**Q: 某些K值运行失败**
```bash
# 检查日志文件
cat admixture_results/logs/error.log

# 常见原因：样本数少于K值，需降低最大K值
biopytools admixture -v input.vcf -o results -K 5
```

**Q: 结果不收敛**
```bash
# ADMIXTURE可能需要多次运行，这是正常现象
# 可以重复运行相同命令，选择最佳结果
```

## 📊 结果解读指南 | Result Interpretation Guide

### CV误差分析 | Cross-Validation Error Analysis

- **最优K值**: CV误差最小的K值通常为最佳选择
- **平台期**: CV误差趋于平稳时，增加K值意义不大
- **生物学验证**: 结合已知的群体历史和地理分布验证结果

### 祖先成分解读 | Ancestry Component Interpretation

- **主要成分**: 比例>50%的成分通常代表个体主要祖先
- **混合个体**: 多个成分比例相近的个体可能来自混合群体
- **群体特异性**: 某些成分在特定地理群体中高频出现

## 📚 相关资源 | Related Resources

- [ADMIXTURE软件手册](https://dalexander.github.io/admixture/admixture-manual.pdf)
- [群体结构分析最佳实践](https://www.nature.com/articles/nrg2813)
- [交叉验证在群体遗传学中的应用](https://www.genetics.org/content/195/3/693)
- [PLINK格式说明](https://www.cog-genomics.org/plink/1.9/formats)

## 📄 许可证 | License

本项目采用MIT许可证 - 详见 [LICENSE](LICENSE) 文件

---

## 🔬 引用信息 | Citation

如果在学术研究中使用此工具，请引用相关方法学文献：

```
Alexander, D. H., Novembre, J., & Lange, K. (2009). 
Fast model-based estimation of ancestry in unrelated individuals. 
Genome research, 19(9), 1655-1664.
```