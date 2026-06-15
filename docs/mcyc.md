# 🧬 MCycDB 甲烷循环基因丰度分析模块

**基于MCycDB数据库的甲烷循环功能基因定量分析工具 | Methane Cycle Gene Abundance Analysis Tool Based on MCycDB Database**

## 📖 功能概述 | Overview

MCycDB 甲烷循环基因丰度分析模块是一个专门用于分析甲烷循环相关基因丰度的工具，基于权威的MCycDB数据库，通过Diamond比对和多种标准化方法，为宏基因组数据中甲烷代谢功能基因的定量分析提供完整解决方案。

## ✨ 主要特性 | Key Features

- **🗄️ 专业数据库**: 基于权威的MCycDB数据库，涵盖完整的甲烷循环功能基因
- **⚡ 高效比对**: 使用Diamond进行快速序列比对，支持多线程并行处理
- **📊 多种标准化**: 支持原始计数、TPM和CLR三种标准化方法
- **🔄 自动流程**: 全自动化分析流程，从数据处理到矩阵计算一键完成
- **🛡️ 质量控制**: 完善的错误处理和数据验证机制
- **📈 结果可视化**: 生成标准化的丰度矩阵，适用于下游统计分析

## 🔬 甲烷循环背景 | Methane Cycle Background

### 甲烷循环的重要性

甲烷循环是全球碳循环的重要组成部分，涉及甲烷的产生（甲烷生成）和消耗（甲烷氧化）两个主要过程：

- **甲烷生成 (Methanogenesis)**: 产甲烷古菌将有机物转化为甲烷
- **甲烷氧化 (Methanotrophy)**: 甲烷氧化细菌将甲烷转化为二氧化碳

### MCycDB数据库

MCycDB是一个专门收集和注释甲烷循环相关基因的数据库，包含：
- 产甲烷途径酶基因
- 甲烷氧化途径酶基因
- 甲烷代谢调控基因
- 辅助代谢功能基因

## 🚀 快速开始 | Quick Start

### 基本用法 | Basic Usage

```bash
# 基本甲烷循环基因丰度分析
biopytools mcyc samples.txt

# 指定输出目录和线程数
biopytools mcyc samples.txt -o ./mcyc_results -t 8

# 跳过Diamond比对（如果已有中间结果）
biopytools mcyc samples.txt --skip-diamond

# 保留临时文件用于调试
biopytools mcyc samples.txt --keep-temp
```

### 输入文件格式 | Input File Format

创建一个制表符分隔的样本文件：

```txt
sample1    /path/to/sample1_R1.fastq.gz
sample2    /path/to/sample2_R1.fastq.gz;/path/to/sample2_R2.fastq.gz
sample3    /path/to/sample3.merged.fastq.gz
```

**文件格式说明**：
- 第一列：样本名称（不能包含特殊字符）
- 第二列：FastQ文件路径
  - 单端测序：直接指定文件路径
  - 双端测序：用分号分隔R1和R2文件路径
- 支持压缩格式（.fastq.gz, .fq.gz）

## 📋 参数说明 | Parameters

### 基本参数 | Basic Parameters

| 参数 | 类型 | 描述 |
|------|------|------|
| `input-list` | 必需 | 📄 输入样本文件路径 |
| `-o, --output` | 可选 | 📁 输出目录（默认：当前目录） |

### 分析参数 | Analysis Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-t, --threads` | `4` | 🧵 并行线程数 |
| `--mcyc-base` | `auto` | 🗄️ MCycDB数据库基础路径 |
| `--skip-diamond` | `False` | ⏭️ 跳过Diamond比对步骤 |
| `--keep-temp` | `False` | 💾 保留临时文件 |

## 🧮 分析流程 | Analysis Pipeline

### 步骤1: 依赖检查 | Step 1: Dependency Check
检查必需的工具：
- **Diamond**: 序列比对工具
- **SeqKit**: 序列统计工具
- **Perl**: MCycDB脚本运行环境

### 步骤2: 数据库准备 | Step 2: Database Preparation
- 检查MCycDB数据库文件完整性
- 构建Diamond索引（如不存在）
- 创建数据库文件软链接

### 步骤3: 数据预处理 | Step 3: Data Preprocessing
- 处理输入样本文件列表
- 创建统一的FastQ文件目录
- 合并双端测序文件（如适用）

### 步骤4: Read统计 | Step 4: Read Counting
使用SeqKit统计每个样本的序列数量：
```bash
seqkit stats -j 4 samples/*.fastq.gz -T
```

### 步骤5: Diamond比对 | Step 5: Diamond Alignment
- 配置MCycDB Perl脚本
- 执行Diamond序列比对
- 生成原始计数矩阵

### 步骤6: 矩阵计算 | Step 6: Matrix Calculation
计算三种标准化矩阵：

#### 6.1 原始计数 (Raw Counts)
- 保持原始read计数
- 输出：`Matrix_1_RawCounts.txt`

#### 6.2 TPM标准化 (Transcripts Per Million)
```python
TPM = (RawCount / LibrarySize) * 1e6
```
- 消除测序深度差异
- 适用于跨样本比较
- 输出：`Matrix_2_TPM.txt`

#### 6.3 CLR变换 (Centered Log Ratio)
```python
CLR(x) = ln(x / geometric_mean(x))
```
- 添加伪计数避免log(0)
- 适用于组成数据分析
- 特别适合GWAS等统计分析
- 输出：`Matrix_3_CLR.txt`

## 💡 使用示例 | Usage Examples

### 示例1: 基础单样本分析 | Example 1: Basic Single Sample Analysis

```bash
# 创建样本文件
cat > single_sample.txt << EOF
sample1    /data/methanogens_metagenome.fastq.gz
EOF

# 运行分析
biopytools mcyc single_sample.txt -o ./results
```

### 示例2: 多样本批量分析 | Example 2: Multi-Sample Batch Analysis

```bash
# 创建多样本文件
cat > multiple_samples.txt << EOF
soil_sample1    /data/soil1_R1.fastq.gz;/data/soil1_R2.fastq.gz
soil_sample2    /data/soil2_R1.fastq.gz;/data/soil2_R2.fastq.gz
wetland_sample1 /data/wetland1.merged.fastq.gz
wetland_sample2 /data/wetland2.merged.fastq.gz
EOF

# 高性能分析
biopytools mcyc multiple_samples.txt -o ./batch_results -t 16
```

### 示例3: 大规模数据分析 | Example 3: Large-Scale Data Analysis

```bash
# 大规模数据处理（跳过临时文件清理）
biopytools mcyc large_cohort.txt \
    -o ./large_scale_results \
    -t 32 \
    --keep-temp

# 如果分析中断，可以跳过Diamond比对继续
biopytools mcyc large_cohort.txt \
    -o ./large_scale_results \
    --skip-diamond \
    --keep-temp
```

### 示例4: 自定义数据库路径 | Example 4: Custom Database Path

```bash
# 使用自定义MCycDB安装路径
biopytools mcyc samples.txt \
    --mcyc-base /path/to/custom/MCycDB \
    -o ./custom_results
```

## 📁 输出文件说明 | Output Files Description

### 主要输出文件 | Main Output Files

```
output_directory/
├── Matrix_1_RawCounts.txt    # 原始计数矩阵
├── Matrix_2_TPM.txt          # TPM标准化矩阵
├── Matrix_3_CLR.txt          # CLR变换矩阵
└── ...
```

### 矩阵文件格式 | Matrix File Format

所有矩阵文件均为制表符分隔格式：

```
GeneID    Sample1    Sample2    Sample3    ...
mcrA      123        456        789       ...
pmoA      234        567        890       ...
mtrA      345        678        901       ...
```

### 临时文件 | Temporary Files

默认情况下，临时文件会被自动清理。使用 `--keep-temp` 可以保留用于调试：

```
working_directory/
├── mcyc_staging_fastq/        # 统一FastQ文件目录
├── MCyc_Raw_Temp.txt          # Diamond原始比对结果
├── sample_info.txt            # Read统计信息
├── MCycDB_FunctionProfiler.PL # 修改的Perl脚本
└── MCycDB_*.dmnd              # 数据库软链接
```

## 🔧 系统要求 | System Requirements

### 依赖软件 | Dependencies

**必需工具**：
- **Diamond** (版本 2.0+)
  - 高性能序列比对工具
  - 下载：https://github.com/bbuchfink/diamond
- **SeqKit** (版本 2.0+)
  - 序列处理和统计工具
  - 下载：https://github.com/shenwei356/seqkit
- **Perl** (版本 5.10+)
  - 运行MCycDB脚本
- **Python** (版本 3.7+)

**Python包**：
```bash
pip install numpy pandas click
```

### 安装依赖 | Installing Dependencies

```bash
# 安装Diamond
wget https://github.com/bbuchfink/diamond/releases/download/v2.0.15/diamond-linux64.tar.gz
tar xzf diamond-linux64.tar.gz
sudo cp diamond /usr/local/bin/

# 安装SeqKit
wget https://github.com/shenwei356/seqkit/releases/download/v2.5.1/seqkit_linux_amd64.tar.gz
tar xzf seqkit_linux_amd64.tar.gz
sudo cp seqkit /usr/local/bin/

# 检查Perl
perl --version

# 安装Python包
pip install numpy pandas click
```

### MCycDB数据库 | MCycDB Database

确保MCycDB数据库已正确安装：

```bash
# 检查MCycDB文件
ls -la /share/org/YZWL/yzwl_lixg/software/MCycDB/MCycDB-main/
# 应该包含：
# - MCycDB_2021.faa
# - id2gene.map
# - MCycDB_FunctionProfiler.PL
```

### 硬件建议 | Hardware Recommendations

- **CPU**: 多核处理器（推荐8核以上）
- **RAM**: 最少16GB（大数据集推荐64GB以上）
- **存储**: 至少预留数据集大小10倍的磁盘空间
- **临时存储**: 额外10GB用于临时文件

## ⚠️ 注意事项 | Important Notes

1. **数据质量**: 确保输入FastQ文件质量良好
2. **文件格式**: 严格遵循样本文件格式要求
3. **路径规范**: 使用绝对路径避免路径问题
4. **权限设置**: 确保对数据库文件和工作目录有读写权限
5. **内存管理**: 大数据集分析时注意内存使用情况

## 🐛 故障排除 | Troubleshooting

### 常见问题 | Common Issues

**Q: "Diamond not found" 错误**
```bash
# 检查Diamond安装
which diamond
diamond version

# 重新安装或添加到PATH
export PATH=$PATH:/path/to/diamond
```

**Q: "MCycDB files not found" 错误**
```bash
# 检查数据库路径
ls -la /share/org/YZWL/yzwl_lixg/software/MCycDB/MCycDB-main/

# 使用自定义路径
biopytools mcyc samples.txt --mcyc-base /correct/path/to/MCycDB
```

**Q: 内存不足错误**
```bash
# 减少线程数
biopytools mcyc samples.txt -t 2

# 分批处理样本
# 将大样本列表拆分为小文件
```

**Q: Diamond索引构建失败**
```bash
# 手动构建索引
cd /share/org/YZWL/yzwl_lixg/software/MCycDB/MCycDB-main/
diamond makedb --in MCycDB_2021.faa --db MCycDB_2021 --quiet
```

**Q: Perl脚本执行失败**
```bash
# 检查Perl安装
perl --version

# 检查脚本权限
ls -la /share/org/YZWL/yzwl_lixg/software/MCycDB/MCycDB-main/MCycDB_FunctionProfiler.PL
chmod +x /path/to/MCycDB_FunctionProfiler.PL
```

**Q: 样本文件解析错误**
```bash
# 检查文件格式
cat samples.txt
# 确保使用制表符分隔，每行两列

# 检查文件路径
cat samples.txt | cut -f2 | xargs ls -la
```

## 📊 结果解读指南 | Result Interpretation Guide

### 矩阵选择 | Matrix Selection

**Matrix_1_RawCounts.txt** - 原始计数
- 适用于需要绝对计数的情况
- 可能受测序深度影响
- 用于简单的基因存在性分析

**Matrix_2_TPM.txt** - TPM标准化
- 消除测序深度差异
- 适用于跨样本丰度比较
- 最常用的标准化方法

**Matrix_3_CLR.txt** - CLR变换
- 适用于组成数据分析
- 推荐用于GWAS、PCA等统计分析
- 适用于多元统计分析

### 基因功能解读 | Gene Function Interpretation

MCycDB中的主要甲烷循环基因类型：

- **mcrA/mcrB/mcrG**: 甲基辅酶M还原酶（甲烷生成关键酶）
- **pmoA/pmoB/pmoC**: 颗粒状甲烷单加氧酶（甲烷氧化关键酶）
- **mxaF**: 甲烷脱氢酶
- **mtrA/mtrB/mtrC**: 甲基转移酶

### 丰度分析 | Abundance Analysis

```python
# 使用Python进行简单分析示例
import pandas as pd
import matplotlib.pyplot as plt

# 读取TPM矩阵
tpm_df = pd.read_csv("Matrix_2_TPM.txt", sep='\t', index_col=0)

# 计算样本间相关性
correlation_matrix = tpm_df.T.corr()

# 可视化关键基因丰度
key_genes = ['mcrA', 'pmoA', 'mxaF']
tpm_df.loc[key_genes].T.plot(kind='bar')
plt.title('Methane Cycle Gene Abundance')
plt.ylabel('TPM')
plt.show()
```

## 📚 相关资源 | Related Resources

### 学术文献 | Academic Papers

- [Methane cycle in ecosystems: a review](https://doi.org/10.1038/nrmicro3378)
- [MCycDB: A database of methane cycle genes](https://doi.org/10.1093/nar/gkaa893)
- [Diamond: Fast and sensitive protein alignment](https://doi.org/10.1038/nmeth.3176)
- [CLR transformation for compositional data analysis](https://doi.org/10.1016/j.chemolab.2003.11.015)

### 相关工具 | Related Tools

- [HUMAnN3](https://github.com/biobakery/humann) - 宏基因组功能分析
- [MetaPhlAn3](https://github.com/biobakery/MetaPhlAn) - 宏基因组物种分析
- [KEGG](https://www.kegg.jp/) - 代谢通路数据库
- [eggNOG](http://eggnog.embl.de/) - 功能注释数据库

### 甲烷循环研究 | Methane Cycle Research

- [Global Methane Initiative](https://www.globalmethane.org/)
- [International Methane Emissions Observatory](https://www.globalmethane.org/imeo)
- [Methane Cycle Research](https://www.nature.com/subjects/methane-cycle)

## 📄 许可证 | License

本项目采用MIT许可证 - 详见 [LICENSE](LICENSE) 文件

MCycDB数据库的使用请遵循原始数据库许可证要求。

---

## 🔬 引用信息 | Citation

如果在学术研究中使用此工具，请引用相关方法学文献：

```
Diamond 2.0:
Buchfink B, Reuter K, Drost HG. (2021)
Fast and sensitive protein alignment using DIAMOND.
Nat Methods 18:366-368.

MCycDB:
Yang L, Xie X, Liu F, et al. (2020)
MCycDB: a database of methane cycle genes for metagenomic profiling.
Nucleic Acids Res 48:D830-D837.
```