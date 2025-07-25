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

---

## 📦 安装 | Installation

### 环境要求 | Requirements
- Python 3.10
- 推荐使用conda/mamba进行环境管理 | Recommended to use conda/mamba for environment management

### 从源码安装 | Install from Source

```bash
# 创建虚拟环境 | Create virtual environment
mamba create -n biopytools python=3.10 
mamba activate biopytools

# 安装依赖 | Install dependencies
mamba install conda-forge::python-rocksdb

# 克隆并安装 | Clone and install
git clone https://github.com/lixiang117423/biopytools.git
cd biopytools
pip install -e .
```

### 安装开发依赖 | Install Development Dependencies

```bash
pip install -e ".[dev]"
```

### 验证安装 | Verify Installation

```bash
# 检查版本 | Check version
python -c "import biopytools; biopytools.show_version()"

# 检查可用模块 | Check available modules
python -c "import biopytools; biopytools.list_modules()"
```

---

## 🛠️ 功能模块 | Available Modules

### 数据质控与预处理 | Data Quality Control & Preprocessing

#### **run_fastp** - 高速FASTQ质控工具 | High-speed FASTQ Quality Control
```bash
run_fastp -i input_dir -o output_dir -t 16
```
**功能** | **Features**: 
- 自动适应接头去除 | Adaptive adapter trimming
- 质量过滤和长度过滤 | Quality and length filtering  
- HTML和JSON报告生成 | HTML and JSON report generation

#### **run_rnaseq** - RNA测序分析流程 | RNA-seq Analysis Pipeline
```bash
run_rnaseq -g genome.fa -a annotation.gtf -i fastq_dir -o output_dir
```
**功能** | **Features**:
- HISAT2序列比对 | HISAT2 alignment
- StringTie转录本定量 | StringTie transcript quantification
- 表达矩阵合并 | Expression matrix merging

### 基因组分析 | Genomic Analysis

#### **run_minimap2** - 长读序列比对工具 | Long-read Alignment Tool
```bash
run_minimap2 -r reference.fa -i reads.fastq -o output_dir
```
**功能** | **Features**:
- 支持多种测序技术 | Support various sequencing technologies
- 高效的长读序列比对 | Efficient long-read alignment
- 灵活的参数配置 | Flexible parameter configuration

#### **run_repeat_masker** - 重复序列分析 | Repeat Sequence Analysis
```bash
run_repeat_masker -g genome.fa -s species -o output_dir
```
**功能** | **Features**:
- RepeatMasker重复序列屏蔽 | RepeatMasker repeat masking
- RepeatModeler从头预测 | RepeatModeler de novo prediction
- TRF串联重复检测 | TRF tandem repeat detection

#### **run_augustus_multi_rnaseq** - 多转录组基因预测 | Multi-RNA-seq Gene Prediction
```bash
run_augustus_multi_rnaseq -g genome.fa -i fastq_dir -s species
```
**功能** | **Features**:
- 结合多个转录组数据 | Integrate multiple RNA-seq datasets
- Augustus基因结构预测 | Augustus gene structure prediction
- 自动化流程管理 | Automated pipeline management

### 变异分析 | Variant Analysis

#### **run_vcf_extractor** - VCF基因型提取 | VCF Genotype Extraction
```bash
run_vcf_extractor -v input.vcf -o output_dir
```
**功能** | **Features**:
- 高效的基因型数据提取 | Efficient genotype data extraction
- 支持大规模VCF文件 | Support large-scale VCF files
- 多种输出格式 | Multiple output formats

#### **parse_sequence_vcf** - VCF序列工具包 | VCF Sequence Toolkit
```bash
parse_sequence_vcf -v input.vcf -g genome.fa -o output
```
**功能** | **Features**:
- VCF变异序列解析 | VCF variant sequence parsing
- 序列坐标转换 | Sequence coordinate conversion
- 变异效应预测 | Variant effect prediction

#### **run_vcf_pca** - VCF主成分分析 | VCF Principal Component Analysis
```bash
run_vcf_pca -v input.vcf -o output_dir
```
**功能** | **Features**:
- 群体遗传结构分析 | Population genetic structure analysis
- 主成分分析可视化 | PCA visualization
- 样本聚类分析 | Sample clustering analysis

### 群体遗传学分析 | Population Genetics Analysis

#### **run_admixture** - 群体结构分析 | Population Structure Analysis
```bash
run_admixture -v input.vcf -o output_dir -k 2 -K 10
```
**功能** | **Features**:
- ADMIXTURE群体结构推断 | ADMIXTURE population structure inference
- 交叉验证最优K值选择 | Cross-validation for optimal K selection
- 协变量生成和可视化 | Covariate generation and visualization

#### **run_plink_gwas** - 全基因组关联分析 | Genome-wide Association Study
```bash
run_plink_gwas -v input.vcf -p phenotype.txt -o output_dir
```
**功能** | **Features**:
- PLINK GWAS分析流程 | PLINK GWAS analysis pipeline
- 群体分层控制 | Population stratification control
- 曼哈顿图和QQ图生成 | Manhattan and QQ plot generation

#### run_vcf_ld_heatmap - VCF连锁不平衡热图生成器 | VCF Linkage Disequilibrium Heatmap Generator
```
bashrun_vcf_ld_heatmap -i variants.vcf -o ld_heatmap.png
```
功能 | Features:

从VCF文件生成连锁不平衡热图 | Generate LD heatmaps from VCF files
支持区域指定和样本过滤 | Support region specification and sample filtering
可自定义图形参数和输出格式 | Customizable graphics parameters and output formats

### 基因组注释与分析 | Genome Annotation & Analysis

#### **run_annovar** - 变异注释工具 | Variant Annotation Tool
```bash
run_annovar -v input.vcf -d database_dir -o output_dir
```
**功能** | **Features**:
- ANNOVAR变异功能注释 | ANNOVAR variant functional annotation
- 多数据库整合注释 | Multi-database integrated annotation
- 可定制注释流程 | Customizable annotation workflow

#### **parse_gene_info** - 基因信息解析 | Gene Information Parsing
```bash
parse_gene_info -g annotation.gff -o output_dir
```
**功能** | **Features**:
- GFF/GTF文件解析 | GFF/GTF file parsing
- 基因结构信息提取 | Gene structure information extraction
- 注释信息标准化 | Annotation information standardization

### K-mer分析 | K-mer Analysis

#### **run_kmer_analysis** - K-mer频率分析 | K-mer Frequency Analysis
```bash
run_kmer_analysis -i input.fasta -k 21 -o output_dir
```
**功能** | **Features**:
- 高效的K-mer计数 | Efficient K-mer counting
- 频率分布统计 | Frequency distribution statistics
- 可视化分析结果 | Visualization of analysis results

#### **run_kmer_pav** - K-mer存在/缺失变异分析 | K-mer Presence/Absence Variation Analysis
```bash
run_kmer_pav -i genome_list.txt -k 31 -o output_dir
```
**功能** | **Features**:
- 基于K-mer的PAV检测 | K-mer based PAV detection
- 基因组比较分析 | Comparative genomic analysis
- 结构变异识别 | Structural variation identification

### 实用工具 | Utility Tools

#### **parse_longest_mrna** - 最长转录本提取 | Longest Transcript Extraction
```bash
parse_longest_mrna -g annotation.gtf -o longest.gtf
```
**功能** | **Features**:
- 每个基因最长转录本选择 | Select longest transcript per gene
- 注释文件简化 | Annotation file simplification
- 下游分析准备 | Preparation for downstream analysis

#### **parse_sample_hete** - 样本杂合度统计 | Sample Heterozygosity Statistics
```bash
parse_sample_hete -v input.vcf -o heterozygosity.txt
```
**功能** | **Features**:
- 个体杂合度计算 | Individual heterozygosity calculation
- 群体遗传多样性评估 | Population genetic diversity assessment
- 统计结果可视化 | Statistical result visualization

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

### 流程组合 | Pipeline Combination

```bash
# 完整的RNA-seq分析流程 | Complete RNA-seq analysis pipeline
# 1. 数据质控 | Data quality control
run_fastp -i raw_fastq -o clean_fastq -t 16

# 2. 序列比对和定量 | Alignment and quantification  
run_rnaseq -g genome.fa -a annotation.gtf -i clean_fastq -o rnaseq_results

# 3. 差异表达分析 (后续分析) | Differential expression analysis (downstream)
```

---

## 📊 输出结果 | Output Results

每个模块都会生成详细的分析报告和结果文件，包括：

Each module generates detailed analysis reports and result files, including:

- **日志文件** | **Log Files**: 详细的运行日志和错误信息 | Detailed run logs and error messages
- **统计报告** | **Statistical Reports**: 数据处理和分析统计 | Data processing and analysis statistics  
- **可视化图表** | **Visualization Charts**: 结果图表和质量控制图 | Result charts and quality control plots
- **标准格式输出** | **Standard Format Output**: 兼容下游分析的标准文件格式 | Standard file formats compatible with downstream analysis

---

## 🔧 依赖软件 | Dependencies

### 必需软件 | Required Software
- **fastp**: FASTQ质控 | FASTQ quality control
- **HISAT2**: RNA-seq比对 | RNA-seq alignment  
- **StringTie**: 转录本定量 | Transcript quantification
- **minimap2**: 长读序列比对 | Long-read alignment
- **PLINK**: 遗传分析 | Genetic analysis
- **ADMIXTURE**: 群体结构分析 | Population structure analysis
- **ANNOVAR**: 变异注释 | Variant annotation
- **Augustus**: 基因预测 | Gene prediction
- **RepeatMasker**: 重复序列分析 | Repeat sequence analysis

### Python依赖 | Python Dependencies
```
pandas>=1.0.0
numpy>=1.19.0
pyfastx>=0.8.4
python-rocksdb
scikit-learn>=1.0.0
matplotlib>=3.5.0
seaborn>=0.11.0
```

---

## 📖 文档和示例 | Documentation & Examples

- **完整文档** | **Complete Documentation**: [https://lixiang117423.github.io/article/biopytools-readme/](https://lixiang117423.github.io/article/biopytools-readme/)
- **使用示例** | **Usage Examples**: 查看各模块目录下的README文件 | Check README files in each module directory
- **API参考** | **API Reference**: 查看代码中的docstring文档 | Check docstring documentation in code

---

## 🐛 问题反馈 | Issues & Support

如果您在使用过程中遇到问题，请通过以下方式联系我们：

If you encounter any issues during usage, please contact us through:

- 📧 **邮箱** | **Email**: lixiang117423@gmail.com
- 🐛 **GitHub Issues**: [https://github.com/lixiang117423/biopytools/issues](https://github.com/lixiang117423/biopytools/issues)
- 💬 **讨论区** | **Discussions**: GitHub Discussions页面 | GitHub Discussions page

提交问题时，请包含以下信息：

When submitting issues, please include:

- 操作系统和Python版本 | Operating system and Python version
- 完整的错误信息 | Complete error messages
- 重现问题的最小示例 | Minimal example to reproduce the issue
- 输入数据的基本信息 | Basic information about input data

---

## 🤝 贡献指南 | Contributing

我们欢迎社区贡献！您可以通过以下方式参与项目：

We welcome community contributions! You can participate in the project through:

### 代码贡献 | Code Contributions
1. Fork项目仓库 | Fork the repository
2. 创建功能分支 | Create a feature branch
3. 提交您的更改 | Commit your changes
4. 创建Pull Request | Create a Pull Request

### 文档改进 | Documentation Improvements
- 修正文档错误 | Fix documentation errors
- 添加使用示例 | Add usage examples  
- 翻译文档内容 | Translate documentation

### 问题报告 | Issue Reporting
- 报告bugs | Report bugs
- 提出功能请求 | Suggest feature requests
- 改进现有功能 | Improve existing features

---

## 📄 许可证 | License

本项目使用MIT许可证 - 查看 [LICENSE](LICENSE) 文件了解详情

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## 📈 更新日志 | Changelog

### v1.18.0
- ✨ 增加run_vcf_njtree模块 | Added run_vcf_njtree module
- 🐛 修复run_admixture模块染色体编号问题 | Fix  chromosome numbering issue in run_admixture

### v1.17.0
- ✨ 增加run_augustus_prediction模块 | Added run_augustus_prediction

### v1.16.0
- ✨ 增加run_ena_downloader模块 | Added run_ena_downloader module

### v1.15.0
- ✨ 增加run_vcf_ld_heatmap模块 | Added run_vcf_ld_heatmap module

### v1.14.1
- ✨ 增加run_vcf_pca模块 | Added run_vcf_pca module
- 🐛 修复parse_sequence_vcf和run_minimap2模块的序列坐标问题 | Fixed sequence coordinate issues in parse_sequence_vcf and run_minimap2 modules

### v1.13.0
- ✨ 增加parse_sequence_vcf模块 | Added parse_sequence_vcf module

### v1.12.0  
- ✨ 增加run_repeat_masker模块 | Added run_repeat_masker module

### v1.11.0
- ✨ 增加run_minimap2模块 | Added run_minimap2 module

### v1.10.0
- ✨ 增加run_kmer_pav模块 | Added run_kmer_pav module

### v1.9.0
- ✨ 增加run_admixture模块 | Added run_admixture module

### v1.8.0
- ✨ 增加run_augustus_multi_rnaseq模块 | Added run_augustus_multi_rnaseq module

### v1.7.0
- ✨ 增加parse_sample_hete模块 | Added parse_sample_hete module

### v1.6.0
- ✨ 增加parse_longest_mrna模块 | Added parse_longest_mrna module

### v1.5.0
- ✨ 增加run_kmer_analysis模块 | Added run_kmer_analysis module

### v1.4.0
- ✨ 增加run_plink_gwas模块 | Added run_plink_gwas module

### v1.3.0
- ✨ 增加parse_gene_info模块 | Added parse_gene_info module
- ✨ 增加run_annovar模块 | Added run_annovar module  
- ✨ 增加run_vcf_extractor模块 | Added run_vcf_extractor module

### v1.1.0
- ✨ 增加run_fastp模块 | Added run_fastp module
- ✨ 增加run_rnaseq模块 | Added run_rnaseq module

---

## 🌟 致谢 | Acknowledgments

感谢所有为biopytools项目做出贡献的开发者和用户！

Thanks to all developers and users who have contributed to the biopytools project!

特别感谢以下开源项目：

Special thanks to the following open source projects:

- fastp, HISAT2, StringTie, minimap2, PLINK, ADMIXTURE, ANNOVAR, Augustus, RepeatMasker等优秀的生物信息学工具
- pandas, numpy, matplotlib等Python科学计算库

---

<div align="center">

**⭐ 如果这个项目对您有帮助，请给我们一个Star！**

**⭐ If this project helps you, please give us a Star!**

[🏠 项目主页 | Homepage](https://github.com/lixiang117423/biopytools) • [📚 文档 | Documentation](https://lixiang117423.github.io/article/biopytools-readme/) • [🐛 问题反馈 | Issues](https://github.com/lixiang117423/biopytools/issues)

</div>