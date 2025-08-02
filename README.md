# Biopytools

**生物信息学常用工具包** | **Bioinformatics Common Toolkits**

一个集成了多种生物信息学分析工具的Python包，提供统一的API和命令行接口。支持从数据质控到高级分析的完整生物信息学分析流程。

A Python package integrating various bioinformatics analysis tools with unified API and command-line interfaces. Supports complete bioinformatics analysis workflows from data quality control to advanced analyses.

---

## ✨ 主要特性 | Key Features

- 🔧 **统一接口** | **Unified Interface**: 为多种生物信息学工具提供一致的API和命令行接口
- 🚀 **高效处理** | **Efficient Processing**: 支持多线程和批量处理，提高分析效率
- 📊 **可视化结果** | **Visualization**: 内置图表生成和分析报告功能
- 🔄 **模块化设计** | **Modular Design**: 独立的功能模块，可单独使用或组合使用
- 📋 **完整流程** | **Complete Pipeline**: 从原始数据到最终结果的端到端解决方案
- 🌳 **系统发育分析** | **Phylogenetic Analysis**: 支持VCF到系统发育树的完整分析流程
- 🧬 **基因预测** | **Gene Prediction**: 集成Augustus基因预测与训练流程

---

## 📦 安装 | Installation

### 环境要求 | Requirements
- Python 3.10+
- 推荐使用conda/mamba进行环境管理 | Recommended to use conda/mamba for environment management

### 从PyPI安装 | Install from PyPI

```bash
pip install biopytools
```

### 从源码安装 | Install from Source

```bash
# 创建虚拟环境 | Create virtual environment
conda create -n biopytools python=3.10 
conda activate biopytools

# 克隆并安装 | Clone and install
git clone https://github.com/lixiang117423/biopytools.git
cd biopytools
pip install -e .
```

### 安装依赖 | Install Dependencies

核心依赖包将自动安装：
```bash
# 核心依赖 | Core dependencies
pyfastx>=0.8.4
python-rocksdb
pandas>=1.3.0
numpy>=1.20.0

# 可选分析依赖 | Optional analysis dependencies
scikit-learn>=1.0.0
matplotlib>=3.5.0
seaborn>=0.11.0
```

### 验证安装 | Verify Installation

```bash
# 检查版本 | Check version
python -c "import biopytools; biopytools.show_version()"

# 检查可用模块 | Check available modules
python -c "import biopytools; biopytools.list_modules()"

# 检查依赖 | Check dependencies
python -c "import biopytools; biopytools.check_dependencies()"
```

---

## 🛠️ 功能模块 | Available Modules

### 数据获取与质控 | Data Acquisition & Quality Control

#### **run_ena_downloader** - ENA数据下载工具
```bash
run_ena_downloader -a PRJNA661210 -p ftp -m save
```
**功能**: 自动元数据下载、多协议支持(FTP/Aspera)、批量FASTQ文件下载

#### **run_fastp** - 高速FASTQ质控工具
```bash
run_fastp -i input_dir -o output_dir -t 16
```
**功能**: 自动适应接头去除、质量过滤和长度过滤、HTML和JSON报告生成

#### **run_rnaseq** - RNA测序分析流程
```bash
run_rnaseq -g genome.fa -a annotation.gtf -i fastq_dir -o output_dir
```
**功能**: HISAT2序列比对、StringTie转录本定量、表达矩阵合并

### 基因组分析 | Genomic Analysis

#### **run_minimap2** - 长读序列比对工具
```bash
run_minimap2 -r reference.fa -i reads.fastq -o output_dir
```
**功能**: 支持多种测序技术、高效的长读序列比对、灵活的参数配置

#### **run_repeat_masker** - 重复序列分析
```bash
run_repeat_masker -g genome.fa -s species -o output_dir
```
**功能**: RepeatMasker重复序列屏蔽、RepeatModeler从头预测、TRF串联重复检测

### 基因预测与注释 | Gene Prediction & Annotation

#### **run_augustus_prediction** - Augustus基因预测训练流程
```bash
run_augustus_prediction -s MySpecies -g genome.fa -a annotations.gff3
```
**功能**: 完整的Augustus模型训练、自动数据集拆分和验证、预测性能评估

#### **run_augustus_multi_rnaseq** - 多转录组基因预测
```bash
run_augustus_multi_rnaseq -g genome.fa -i fastq_dir -s species
```
**功能**: 结合多个转录组数据、Augustus基因结构预测、自动化流程管理

#### **run_annovar** - 变异注释工具
```bash
run_annovar -v input.vcf -d database_dir -o output_dir
```
**功能**: ANNOVAR变异功能注释、多数据库整合注释、可定制注释流程

### 变异分析 | Variant Analysis

#### **run_vcf_extractor** - VCF基因型提取
```bash
run_vcf_extractor -v input.vcf -o output_dir
```
**功能**: 高效的基因型数据提取、支持大规模VCF文件、多种输出格式

#### **parse_sequence_vcf** - VCF序列工具包
```bash
parse_sequence_vcf -v input.vcf -g genome.fa -o output
```
**功能**: VCF变异序列解析、序列坐标转换、变异效应预测

#### **run_vcf_pca** - VCF主成分分析
```bash
run_vcf_pca -v input.vcf -o output_dir
```
**功能**: 群体遗传结构分析、主成分分析可视化、样本聚类分析

#### **run_vcf_ld_heatmap** - VCF连锁不平衡热图生成器
```bash
run_vcf_ld_heatmap -i variants.vcf -o ld_heatmap.png
```
**功能**: 从VCF文件生成连锁不平衡热图、支持区域指定和样本过滤

### 系统发育分析 | Phylogenetic Analysis

#### **run_vcf_njtree** - VCF系统发育分析
```bash
run_vcf_njtree -i variants.vcf -o phylo_results
```
**功能**: VCF文件到距离矩阵转换、邻接法系统发育树构建、多种输出格式支持

### 群体遗传学分析 | Population Genetics Analysis

#### **run_admixture** - 群体结构分析
```bash
run_admixture -v input.vcf -o output_dir -k 2 -K 10
```
**功能**: ADMIXTURE群体结构推断、交叉验证最优K值选择、协变量生成和可视化

#### **run_plink_gwas** - 全基因组关联分析
```bash
run_plink_gwas -v input.vcf -p phenotype.txt -o output_dir
```
**功能**: PLINK GWAS分析流程、群体分层控制、曼哈顿图和QQ图生成

### K-mer分析 | K-mer Analysis

#### **run_kmer_analysis** - K-mer频率分析
```bash
run_kmer_analysis -i input.fasta -k 21 -o output_dir
```
**功能**: 高效的K-mer计数、频率分布统计、可视化分析结果

#### **run_kmer_pav** - K-mer存在/缺失变异分析
```bash
run_kmer_pav -i genome_list.txt -k 31 -o output_dir
```
**功能**: 基于K-mer的PAV检测、基因组比较分析、结构变异识别

### 实用工具 | Utility Tools

#### **parse_gene_info** - 基因信息解析
```bash
parse_gene_info -g annotation.gff -o output_dir
```
**功能**: GFF/GTF文件解析、基因结构信息提取、注释信息标准化

#### **parse_longest_mrna** - 最长转录本提取
```bash
parse_longest_mrna -g annotation.gtf -o longest.gtf
```
**功能**: 每个基因最长转录本选择、注释文件简化、下游分析准备

#### **parse_sample_hete** - 样本杂合度统计
```bash
parse_sample_hete -v input.vcf -o heterozygosity.txt
```
**功能**: 个体杂合度计算、群体遗传多样性评估、统计结果可视化

---

## 📚 使用指南 | Usage Guide

### 基础用法 | Basic Usage

```python
# Python API 使用示例 | Python API Usage Example
from biopytools.fastp import FastpProcessor, FastpConfig

# 创建配置 | Create configuration
config = FastpConfig(
    input_dir="./raw_data",
    output_dir="./clean_data",
    threads=16,
    quality_threshold=30
)

# 运行分析 | Run analysis
processor = FastpProcessor(config)
results = processor.run_batch_processing()
```

### 命令行用法 | Command Line Usage

```bash
# 查看帮助信息 | View help information
run_fastp --help

# 运行分析 | Run analysis
run_fastp -i raw_data -o clean_data -t 16 -q 30
```

### 流程组合示例 | Pipeline Combination Examples

#### 完整基因组分析流程 | Complete Genome Analysis Pipeline
```bash
# 1. 数据获取 | Data acquisition
run_ena_downloader -a PRJNA661210 -p ftp -m save

# 2. 数据质控 | Data quality control
run_fastp -i raw_fastq -o clean_fastq -t 16

# 3. 基因预测 | Gene prediction
run_augustus_prediction -s MyGenome -g genome.fa -a reference.gff3

# 4. 变异分析 | Variant analysis
run_vcf_pca -v variants.vcf -o pca_results
run_vcf_njtree -i variants.vcf -o phylo_tree
```

#### RNA-seq分析流程 | RNA-seq Analysis Pipeline
```bash
# 1. 质控 | Quality control
run_fastp -i raw_fastq -o clean_fastq -t 16

# 2. 比对和定量 | Alignment and quantification  
run_rnaseq -g genome.fa -a annotation.gtf -i clean_fastq -o rnaseq_results

# 3. 多转录组基因预测 | Multi-RNA-seq gene prediction
run_augustus_multi_rnaseq -g genome.fa -i clean_fastq -s species
```

---

## 📊 输出结果 | Output Results

每个模块都会生成详细的分析报告和结果文件，包括：

- **日志文件**: 详细的运行日志和错误信息
- **统计报告**: 数据处理和分析统计  
- **可视化图表**: 结果图表和质量控制图
- **标准格式输出**: 兼容下游分析的标准文件格式
- **系统发育树**: Newick格式的系统发育树文件

---

## 🔧 依赖软件 | Dependencies

### 必需软件 | Required Software
- **fastp**: FASTQ质控
- **HISAT2**: RNA-seq比对  
- **StringTie**: 转录本定量
- **minimap2**: 长读序列比对
- **PLINK**: 遗传分析
- **ADMIXTURE**: 群体结构分析
- **ANNOVAR**: 变异注释
- **Augustus**: 基因预测
- **RepeatMasker**: 重复序列分析
- **VCF2Dis**: VCF距离矩阵计算
- **BCFtools**: VCF文件处理

### Python依赖 | Python Dependencies
```
pyfastx>=0.8.4
python-rocksdb
pandas>=1.3.0
numpy>=1.20.0
scikit-learn>=1.0.0
matplotlib>=3.5.0
seaborn>=0.11.0
scikit-bio>=0.5.7
scipy
scikit-allel
```

---

## 🔄 持续集成 | Continuous Integration

项目配置了GitHub Actions自动化流程：
- 自动构建和测试
- 自动发布到PyPI
- 使用Trusted Publishing进行安全发布

---

## 📖 文档和示例 | Documentation & Examples

- **项目主页**: [GitHub Repository](https://github.com/lixiang117423/biopytools)
- **使用示例**: 查看各模块目录下的示例文件
- **API参考**: 查看代码中的docstring文档

---

## 🐛 问题反馈 | Issues & Support

如果您在使用过程中遇到问题，请通过以下方式联系我们：

- 📧 **邮箱**: lixiang117423@gmail.com
- 🐛 **GitHub Issues**: [https://github.com/lixiang117423/biopytools/issues](https://github.com/lixiang117423/biopytools/issues)
- 💬 **讨论区**: GitHub Discussions页面

提交问题时，请包含以下信息：
- 操作系统和Python版本
- 完整的错误信息
- 重现问题的最小示例
- 输入数据的基本信息

---

## 🤝 贡献指南 | Contributing

我们欢迎社区贡献！您可以通过以下方式参与项目：

### 代码贡献 | Code Contributions
1. Fork项目仓库
2. 创建功能分支
3. 提交您的更改
4. 创建Pull Request

### 文档改进 | Documentation Improvements
- 修正文档错误
- 添加使用示例  
- 翻译文档内容

### 问题报告 | Issue Reporting
- 报告bugs
- 提出功能请求
- 改进现有功能

---

## 📄 许可证 | License

本项目使用MIT许可证 - 查看 [LICENSE](LICENSE) 文件了解详情

---

## 📈 更新日志 | Changelog

### v1.19.1 (Latest)
- 🐛 替换run_vcf_filter中的PLINK，使用BCFtools

### v1.19.0
- ✨ 增加run_haplotype_extractor模块

### v1.18.3
- ✨ 修复run_popgen_analysis模块软件路径检查bug

### v1.18.0
- ✨ 增加run_vcf_njtree模块
- 🐛 修复run_admixture模块染色体编号问题

### v1.17.0
- ✨ 增加run_augustus_prediction模块

### v1.16.0
- ✨ 增加run_ena_downloader模块

[查看完整更新日志...](CHANGELOG.md)

---

## 🌟 致谢 | Acknowledgments

感谢所有为biopytools项目做出贡献的开发者和用户！

特别感谢以下开源项目：
- fastp, HISAT2, StringTie, minimap2, PLINK, ADMIXTURE, ANNOVAR, Augustus, RepeatMasker, VCF2Dis等优秀的生物信息学工具
- pandas, numpy, matplotlib, scikit-bio等Python科学计算库

---

<div align="center">

**⭐ 如果这个项目对您有帮助，请给我们一个Star！**

[🏠 项目主页](https://github.com/lixiang117423/biopytools) • [📚 文档](https://lixiang117423.github.io/article/biopytools-readme/) • [🐛 问题反馈](https://github.com/lixiang117423/biopytools/issues)

</div>