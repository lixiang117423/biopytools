# 📊 GWAS Lambda GC 计算模块

**高效的GWAS结果质量控制工具 | Efficient GWAS Result Quality Control Tool**

## 📖 功能概述 | Overview

GWAS Lambda GC 模块是一个专门用于评估GWAS分析结果质量的工具，通过计算Lambda GC值来检测群体分层和假阳性问题。该工具支持批量处理多个GWAS结果文件，提供详细的统计报告和质量评估，帮助研究人员快速识别有问题的GWAS分析结果。

## ✨ 主要特性 | Key Features

- **📈 Lambda GC计算**: 基于卡方分布的精确Lambda GC值计算
- **🔍 批量处理**: 支持批量分析多个GWAS结果文件，提高工作效率
- **⚡ 高性能处理**: 优化的数值计算，快速处理大规模数据集
- **📊 详细统计**: 提供显著位点数量、总位点数和置信度评估
- **🎯 质量分级**: 智能评估结果质量（Ideal/Acceptable/Inflated/Deflated）
- **🛡️ 健壮性**: 完善的错误处理和数据验证机制

## 🧠 Lambda GC 原理 | Lambda GC Principle

### 什么是Lambda GC？

Lambda GC（Genomic Control）是衡量GWAS分析中群体膨胀程度的指标：

- **Lambda GC = 1.0**: 无群体分层，结果理想
- **Lambda GC > 1.1**: 存在群体分层或假阳性
- **Lambda GC < 0.9**: 过度校正或模型问题

### 计算公式

```
Lambda GC = 观测卡方值中位数 / 期望卡方值中位数
          = median(χ²(1-p)) / 0.4549364
```

其中：
- `p` 为每个SNP的p值
- `0.4549364` 是自由度为1的卡方分布的理论中位数

## 🚀 快速开始 | Quick Start

### 基本用法 | Basic Usage

```bash
# 使用默认设置分析GWAS结果
biopytools gwas-lambda

# 指定自定义搜索模式
biopytools gwas-lambda -p "results/*/gwas.txt"

# 自定义输出文件和目录
biopytools gwas-lambda -p "gwas_*/*.mlm" -o "my_lambda_assessment.txt" -d "./quality_control"
```

### 高级用法 | Advanced Usage

```bash
# 调整显著性阈值
biopytools gwas-lambda --threshold 1e-6

# 指定P值列（从0开始计数）
biopytools gwas-lambda --p-column 4

# 使用不同显著性阈值和输出配置
biopytools gwas-lambda \
    -p "chr*/chr*_gwas_results.txt" \
    -o "genome_wide_lambda.txt" \
    -t 5e-8 \
    -d "./gwas_quality_control"
```

## 📋 参数说明 | Parameters

### 基本参数 | Basic Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-p, --pattern` | `"feture_*/GWAS_Result.mlm.manht_input"` | 🔍 文件搜索模式（支持通配符） |
| `-o, --output` | `"Batch_Lambda_Assessment.txt"` | 📝 输出文件名 |
| `-d, --output-dir` | `"./gwas_lambda_output"` | 📁 输出目录路径 |

### 分析参数 | Analysis Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-t, --threshold` | `1e-5` | 🎯 显著性阈值（用于统计显著位点） |
| `-c, --p-column` | `3` | 📊 P值所在列索引（0-based，默认第4列） |

## 📁 输入文件格式 | Input File Format

### GWAS结果文件要求 | GWAS Result File Requirements

支持标准制表符分隔的GWAS结果文件：

```
Chromosome	Position	Marker	Allele1	Allele2	P-value	Effect
1	14370	rs6054257	G	A	2.1e-08	0.12
1	17330	.	T	A	1.5e-06	0.34
1	1110696	rs6040355	A	G	3.2e-04	0.08
```

**文件要求**：
- 支持任意列数的GWAS结果文件
- P值列可以通过 `--p-column` 参数指定
- 支持标准GWAS软件的输出格式（TASSEL、PLINK、GEMMA等）
- 自动跳过表头行和非数值数据

### 搜索模式示例 | Search Pattern Examples

```bash
# 基本模式 - 匹配子目录中的特定文件
biopytools gwas-lambda -p "results/*/gwas_results.txt"

# 高级模式 - 匹配多种文件类型
biopytools gwas-lambda -p "gwas_*/*.mlm.manht_input"

# 染色体级别分析
biopytools gwas-lambda -p "chr*/chr*_gwas_results.txt"

# 多条件模式
biopytools gwas-lambda -p "*/*.{mlm,gwas,assoc}"
```

## 💡 使用示例 | Usage Examples

### 示例1：标准TASSEL分析结果评估 | Example 1: Standard TASSEL Results Assessment

```bash
# 评估TASSEL MLM模型结果
biopytools gwas-lambda \
    -p "trait_*/*_GWAS.mlm.manht_input" \
    -o "tasel_mlm_quality.txt" \
    -d "./tasel_quality_control"

# 输出文件将包含：
# - Lambda GC值
# - 显著位点数量（P < 1e-5）
# - 总SNP数量
# - 质量评估状态
```

### 示例2：多方法比较分析 | Example 2: Multi-Method Comparison Analysis

```bash
# 比较不同GWAS方法的结果
biopytools gwas-lambda -p "glm_results/*.txt" -o "glm_quality.txt"
biopytools gwas-lambda -p "mlm_results/*.txt" -o "mlm_quality.txt"
biopytools gwas-lambda -p "farmcpu_results/*.txt" -o "farmcpu_quality.txt"
```

### 示例3：染色体级别质量控制 | Example 3: Chromosome-Level Quality Control

```bash
# 分别评估每个染色体的GWAS结果
biopytools gwas-lambda \
    -p "chr*/chr*_gwas_results.txt" \
    -o "chromosome_wide_lambda.txt" \
    -d "./chr_quality_control"
```

### 示例4：严格显著性阈值评估 | Example 4: Strict Significance Threshold Assessment

```bash
# 使用更严格的显著性阈值
biopytools gwas-lambda \
    -p "large_scale_gwas/*.txt" \
    -o "strict_threshold_assessment.txt" \
    -t 5e-8  # 全基因组显著性阈值
```

### 示例5：自定义格式文件分析 | Example 5: Custom Format Analysis

```bash
# 处理不同列顺序的文件
biopytools gwas-lambda \
    -p "gemma_results/*.assoc.txt" \
    -o "gemma_quality.txt" \
    -c 2  # P值在第3列（索引为2）
```

## 📁 输出文件说明 | Output Files Description

### 主要输出文件 | Main Output File

生成的输出文件包含以下列：

| 列名 | 描述 |
|------|------|
| `Folder` | GWAS结果文件所在文件夹路径 |
| `Lambda_GC` | Lambda GC值（NA表示无法计算） |
| `Sig_SNPs(<1e-5)` | 显著位点数量（P值小于显著性阈值） |
| `Total_SNPs` | 总有效SNP数量 |
| `Status` | 质量评估状态 |

### 质量状态说明 | Quality Status Description

| 状态 | Lambda GC范围 | 建议 |
|------|--------------|------|
| **Ideal (完美)** | 0.95 - 1.05 | ✅ 结果理想，无需调整 |
| **Acceptable (可接受)** | 0.90 - 1.10 | ✅ 结果可接受 |
| **Inflated (膨胀)** | > 1.10 | ⚠️ 存在群体分层，需要校正 |
| **Deflated (压缩)** | 0.80 - 0.90 | ⚠️ 可能过度校正 |
| **Warning: Extreme Deflated** | < 0.80 | ❌ 模型可能失效 |
| **No Signals (无显著关联)** | NA | 💡 无显著位点，正常情况 |

### 示例输出 | Example Output

```
Folder	Lambda_GC	Sig_SNPs(<1e-5)	Total_SNPs	Status
trait_1	1.0234	45	1057153	Acceptable (可接受)
trait_2	1.2345	128	1057153	Inflated (膨胀/假阳性高)
trait_3	NA	0	1057153	No Signals (无显著关联)
trait_4	0.7568	234	1057153	Warning: Extreme Deflated (异常/模型失效)
```

## 📊 结果解读指南 | Result Interpretation Guide

### Lambda GC 值解读 | Lambda GC Value Interpretation

#### 理想范围：0.95 - 1.05
- **含义**: 无明显群体分层，分析结果可靠
- **建议**: 可以直接用于后续分析

#### 可接受范围：0.90 - 1.10
- **含义**: 轻微偏差，但结果基本可靠
- **建议**: 结合其他指标评估，通常可以接受

#### 膨胀：> 1.10
- **含义**: 存在群体分层或假阳性
- **建议**:
  - 使用MLM或混合模型
  - 增加PCA主成分
  - 检查样本质量
  - 考虑使用Genomic Control校正

#### 压缩：< 0.90
- **含义**: 可能过度校正或模型问题
- **建议**:
  - 检查协变量设置
  - 验证模型参数
  - 考虑简化模型

### 结合其他指标 | Combined with Other Metrics

1. **QQ图分析**: 观察偏差曲线形状
2. **显著位点数量**: 评估发现能力
3. **Manhattan图**: 检查峰值分布
4. **样本结构**: PCA或系统发育分析

## 🔧 系统要求 | System Requirements

### 依赖软件 | Dependencies

- **Python** (版本 3.7+)
- **Python包**:
  - `numpy` - 数值计算
  - `scipy` - 统计分析（卡方分布）
  - `pandas` - 数据处理
  - `glob` - 文件搜索（Python内置）
  - `click` - 命令行界面

### 安装依赖 | Installing Dependencies

```bash
# 安装必需的Python包
pip install numpy scipy pandas click

# 或使用conda
conda install numpy scipy pandas
```

### 硬件建议 | Hardware Recommendations

- **CPU**: 多核处理器（推荐4核以上）
- **RAM**: 最少4GB（大规模数据推荐16GB以上）
- **存储**: 预留足够空间存放结果文件
- **网络**: 不需要网络连接

## ⚠️ 注意事项 | Important Notes

1. **数据质量**: 确保GWAS分析结果文件格式正确
2. **P值范围**: P值必须在(0, 1]范围内
3. **显著性阈值**: 根据研究目标选择合适的阈值
4. **批量处理**: 大量文件处理时注意存储空间
5. **结果解读**: Lambda GC只是质量指标之一，需要结合其他方法

## 🐛 故障排除 | Troubleshooting

### 常见问题 | Common Issues

**Q: "未找到匹配文件" 错误**
```bash
# 检查文件路径和模式
ls -la your_pattern_here*

# 使用绝对路径
biopytools gwas-lambda -p "/full/path/to/files/*.txt"
```

**Q: "P值列索引错误"**
```bash
# 检查文件列格式
head -5 your_gwas_file.txt | column -t

# 尝试不同的列索引
biopytools gwas-lambda -p "*.txt" -c 4  # 第5列
biopytools gwas-lambda -p "*.txt" -c 6  # 第7列
```

**Q: Lambda GC值为无穷大**
```bash
# 这通常已修复，如果仍有问题请检查数据
# 检查P值是否有异常小的值
awk 'NR>1 && $NF<1e-300' your_file.txt
```

**Q: 内存不足错误**
```bash
# 对于极大数据集，可以分批处理
# 或使用更强大的计算资源
```

**Q: 输出文件为空**
```bash
# 检查文件权限
ls -la output_directory/

# 检查磁盘空间
df -h .
```

## 📚 相关资源 | Related Resources

### 学术文献 | Academic Papers

- [Devlin B, Roeder K. (1999) Genomic control for association studies.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1378117/)
- [The 1000 Genomes Project Consortium. (2015) A global reference for human genetic variation.](https://www.nature.com/articles/nature15393)
- [Yang J et al. (2011) Genomic inflation factors.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3116406/)

### 相关工具 | Related Tools

- [PLINK](https://www.cog-genomics.org/plink/) - GWAS分析软件
- [GEMMA](https://github.com/genetics-statistical-genomics/gemma) - 混合模型GWAS
- [TASSEL](https://www.maizegenetics.net/tassel/) - 关联分析软件
- [SAIGE](https://github.com/weizhouUMICH/SAIGE) - 大规模GWAS工具

### 数据可视化 | Data Visualization

- [qqman](https://github.com/stephenturner/qqman) - R包，用于GWAS结果可视化
- [CMplot](https://github.com/YinLiLin/R-CMplot) - R包，用于曼哈顿图和QQ图
- [LocusZoom](http://locuszoom.org/) - 在版GWAS结果可视化

## 📄 许可证 | License

本项目采用MIT许可证 - 详见 [LICENSE](LICENSE) 文件

---

## 🔬 引用信息 | Citation

如果在学术研究中使用此工具，请引用相关方法学文献：

```
Devlin B, Roeder K. (1999)
Genomic control for association studies.
Biometrics 55:997-1004.

The 1000 Genomes Project Consortium. (2015)
A global reference for human genetic variation.
Nature 526:68-74.
```