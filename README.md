# BioPyTools

A Python toolkit for bioinformatics analysis and computational biology.

一个用于生物信息学分析和计算生物学的Python工具包。

## 简介 | Overview

BioPyTools 是一个专为生物信息学研究设计的Python工具包，提供了一系列常用的生物数据分析功能。

BioPyTools is a Python toolkit designed for bioinformatics research, providing a series of commonly used biological data analysis functions.

## 系统要求 | Requirements

- Python >= 3.8
- NumPy >= 1.19.0
- Pandas >= 1.2.0
- Matplotlib >= 3.3.0

## Emoji显示支持 | Emoji Display Support

本工具使用emoji来增强日志的可读性和用户体验。如果你看到乱码或方块字符，请按以下步骤设置：

This tool uses emojis to enhance log readability and user experience. If you see garbled text or square characters, please follow these setup steps:

### 终端配置 | Terminal Configuration

**Linux:**
```bash
# 安装emoji字体支持 | Install emoji font support
sudo apt install fonts-noto-color-emoji  # Ubuntu/Debian
sudo yum install google-noto-emoji-fonts # CentOS/RHEL
sudo pacman -S noto-fonts-emoji          # Arch Linux

# 确保locale支持UTF-8 | Ensure UTF-8 locale support
export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8
```

**macOS:**
```bash
# macOS通常自带emoji支持 | macOS usually has built-in emoji support
# 推荐使用iTerm2获得更好的显示效果 | Recommend iTerm2 for better display
brew install --cask iterm2
```

**Windows:**
```powershell
# 推荐使用Windows Terminal | Recommend Windows Terminal
winget install Microsoft.WindowsTerminal

# 设置终端字体为支持emoji的字体 | Set terminal font to emoji-supported font
# 如：Cascadia Code, Segoe UI Emoji
```

### SSH远程服务器 | SSH Remote Server

```bash
# 确保服务器支持UTF-8 | Ensure server supports UTF-8
export LANG=en_US.UTF-8
export LC_ALL=en_US.UTF-8

# 永久设置 | Permanent setup
echo 'export LANG=en_US.UTF-8' >> ~/.bashrc
echo 'export LC_ALL=en_US.UTF-8' >> ~/.bashrc
source ~/.bashrc
```

### 测试emoji显示 | Test Emoji Display

运行以下命令测试emoji是否正确显示 | Run the following command to test emoji display:

```bash
python -c "print('🧬🔍✅❌⚠️📊 Emoji测试 | Emoji Test')"
```

如果仍有显示问题，可以设置环境变量禁用emoji | If you still have display issues, you can disable emojis with:

```bash
export BIOPY_NO_EMOJI=1
```

## 安装方法 | Installation

### 从源码安装 | Install from source

```bash
git clone https://github.com/lixiang117423/biopytools.git
cd biopytools
pip install -e .
```

## 使用方法 | Usage

### 查看帮助 | Getting Help

```bash
biopytools -h
Usage: biopytools [OPTIONS] COMMAND [ARGS]...

  BioPyTools - 生物信息学分析工具包

  要查看特定命令的帮助，请运行：biopytools <命令> -h/--help, 如biopytools fastp -h

Options:
  -v, --version  Show the version and exit.
  -h, --help     Show this message and exit.
```

### 🧬 全流程分析示例 | End-to-End Analysis Examples

#### GPU加速分析 (推荐) | GPU-Accelerated Analysis (Recommended)
```bash
# 完整的GPU加速全流程分析
biopytools fastq2vcf-parabricks \
    -i /path/to/raw_fastq \
    -r /path/to/reference.fa \
    -p /path/to/project

# 跳过质量控制，直接分析
biopytools fastq2vcf-parabricks \
    -i /path/to/raw_fastq \
    -r /path/to/reference.fa \
    -p /path/to/project \
    --skip-qc
```

#### CPU优化分析 | CPU-Optimized Analysis
```bash
# 完整的CPU优化全流程分析
biopytools fastq2vcf-gtx \
    -i /path/to/raw_fastq \
    -r /path/to/reference.fa \
    -p /path/to/project

# 大规模数据集群模式
biopytools fastq2vcf-gtx \
    -i /path/to/raw_fastq \
    -r /path/to/reference.fa \
    -p /path/to/project \
    --threads-gtx 128 \
    --gtx-window-size 10000000
```

#### 🌾 GWAS关联分析 | GWAS Association Analysis
```bash
# 基本GWAS分析 - 自动识别并处理所有表型
biopytools tassel-gwas -i genotype.vcf.gz -p traits.txt -o gwas_results

# 群体结构校正的GWAS分析
biopytools tassel-gwas -i genotype.vcf.gz -p traits.txt -o gwas_results \
    --model MLM --pca-components 10 --maf 0.05

# 高性能并行GWAS分析
biopytools tassel-gwas -i large_dataset.vcf.gz -p traits.txt -o gwas_results \
    --model BOTH --parallel --workers 8 --memory 200g

# GWAS结果质量控制 - 计算Lambda GC值
biopytools gwas-lambda -p "gwas_results/*/*.mlm.manht_input" -o "gwas_quality.txt"

# 使用严格显著性阈值评估GWAS结果
biopytools gwas-lambda -p "gwas_results/*/*.mlm" --threshold 5e-8
```

#### 🧬 基因组组装与挂载 | Genome Assembly & Scaffolding
```bash
# 基于Hi-C数据的基因组scaffolding - Pipeline模式
biopytools haphic -a assembly.fa -b hic.bam -c 24

# 使用原始FASTQ数据（自动执行BWA比对）
biopytools haphic -a assembly.fa -1 hic_R1.fastq.gz -2 hic_R2.fastq.gz -c 24

# 高性能配置，包含组装校正
biopytools haphic -a assembly.fa -b hic.bam -c 24 \
    --threads 32 --processes 16 --correct-nrounds 2

# 单倍型分相组装
biopytools haphic -a phased_assembly.fa -b hic.bam -c 24 \
    --remove-allelic-links 2

# 断点续传（默认启用）
biopytools haphic -a assembly.fa -b hic.bam -c 24
# 中断后再次运行自动跳过已完成步骤

# 强制重新运行所有步骤
biopytools haphic -a assembly.fa -b hic.bam -c 24 --force-rerun
```

### 分步骤执行 | Step-by-Step Execution
```bash
# 执行特定步骤
biopytools fastq2vcf-parabricks -i ./raw -r ./ref.fa -p ./proj --step 3  # GPU模式第3步
biopytools fastq2vcf-gtx -i ./raw -r ./ref.fa -p ./proj --step 4          # CPU模式第4步
```

## 模块文档 | README

### 🧬 全流程分析 | End-to-End Analysis

[Fastq到VCF (Parabricks) | GPU加速全流程分析](./docs/fastq2vcf_parabricks.md) - 🚀 基于NVIDIA Parabricks的GPU加速FASTQ到VCF完整分析流程

[Fastq到VCF (GTX) | CPU优化全流程分析](./docs/fastq2vcf_gtx.md) - 🔬 基于GTX的CPU优化FASTQ到VCF完整分析流程

[Parabricks | GPU加速WGS分析](./docs/parabricks.md) - 🚀 单独的Parabricks GPU加速WGS分析工具

### 📊 群体遗传与进化 | Population Genetics & Evolution

[TASSEL GWAS | 全基因组关联分析](./docs/tassel_gwas.md) - 🌾 基于TASSEL的GWAS分析工具

[GWAS Lambda | GWAS质量控制](./docs/gwas_lambda.md) - 📊 GWAS结果Lambda GC计算和质量评估工具

[MCycDB | 甲烷循环基因丰度分析](./docs/mcyc.md) - 🧬 基于MCycDB的甲烷循环功能基因定量分析工具

[Admixture | 群体结构](./docs/admixture.md) - 🧬 ADMIXTURE群体结构分析工具

[OrthoFinder | 基于基因的泛基因组构建](./docs/orthofinder.md) - 🔍 泛基因组分析和系统发育构建

[VCF系统发育 | VCF系统发育分析](./docs/vcf_phylo.md) - 🌳 基于VCF文件的系统发育分析

[VCF转PHYLIP | VCF格式转换工具](./docs/vcf2phylip.md) - 🔄 VCF到PHYLIP格式转换

### 🧪 序列分析与比对 | Sequence Analysis & Alignment

[Blast | Blast比对](./docs/blast.md) - 🔍 BLAST序列比对分析工具

[Longest mRNA | 提取最长转录本](./docs/longest_mrna.md) - 📜 从转录组数据中提取最长转录本

[Subseq | 序列子集提取工具](./docs/subseq.md) - ✂️ 序列子集提取和分析工具

### 📝 变异检测与注释 | Variant Detection & Annotation

[GTX Joint Calling | 联合变异检测命令生成](./docs/gtx-joint.md) - 🧬 按染色体或区间生成GTX Joint Calling命令脚本

[VCF Renamer | 样品名称重命名](./docs/vcf-renamer.md) - 🏷️ VCF文件样品名称重命名工具，防止名称被软件截断

[Annovar | 变异注释](./docs/annovar.md) - 📝 ANNOVAR变异功能注释工具

### 📥 数据获取与处理 | Data Acquisition & Processing

[SRA2FASTQ | SRA转FASTQ](./docs/sra2fastq.md) - 🧬 SRA数据库数据转FASTQ格式

[从CNCB获取链接 | 批量获取测序数据下载链接](./docs/get_link_from_CNCB.md) - 📥 CNCB数据库批量数据下载工具

### 🧬 基因组组装与挂载 | Genome Assembly & Scaffolding

[HapHiC | 基因组scaffolding工具](./docs/haphic.md) - 🧬 基于Hi-C数据的基因组scaffolding工具
- **🔄 Pipeline模式**: 原生四步自动化流程，一步完成所有scaffolding步骤
- **📍 断点续传**: 自动检测进度，支持中断恢复，避免重复计算
- **📊 自动可视化**: 默认生成Hi-C接触图，便于结果评估
- **🥤 Juicebox集成**: 自动生成Juicebox兼容文件，支持手动校正
- **⚡ 高效处理**: 支持BWA自动比对、智能校正、等位基因分相

[EGAPx Batch | EGAPx批量配置生成](./docs/egapx-batch.md) - 🧬 按染色体批量生成EGAPx运行配置和脚本
- **📂 智能拆分**: 按染色体/scaffold自动拆分基因组
- **📝 批量配置**: 为每个序列生成独立的YAML和脚本
- **🚀 并行执行**: 生成任务列表和并行执行脚本
- **🏷️ 灵活命名**: 支持自定义locus标签和报告名称

[Genome Analysis | 基因组特征评估](./docs/genome_analysis.md) - 🧬 基于GenomeScope2和Smudgeplot的基因组分析工具
- **📏 基因组大小**: 估算基因组大小和重复序列比例
- **🧬 杂合度**: 评估基因组杂合度水平
- **📊 倍性分析**: 使用Smudgeplot推断倍性
- **🔢 K-mer Coverage**: 自动提取用于后续分析

[Rename Chromosomes | 染色体重命名](./docs/rename_chromosomes.md) - 🧬 FASTA序列标准化重命名工具
- **🔢 自动编号**: 前N条命名为Chr01, Chr02...
- **🧱 Scaffold重命名**: 剩余序列命名为HiC_scaffold_XX
- **⚡ 流式处理**: 基于AWK的高效处理，支持大文件
- **📝 日志记录**: 完整的运行日志和统计信息


## 许可证 | License

本项目采用 MIT 许可证 - 详见 [LICENSE](LICENSE) 文件

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 作者信息 | Author

**李详 (Xiang Li)**
- Email: lixiang117423@gmail.com
- GitHub: [@lixiang117423](https://github.com/lixiang117423)

## 致谢 | Acknowledgments

- 感谢所有为本项目做出贡献的开发者 | Thanks to all developers who contributed to this project
- 感谢开源社区的支持 | Thanks to the open source community for support

## 问题反馈 | Issues

如果遇到问题或有建议，请在GitHub上提交issue：

If you encounter problems or have suggestions, please submit an issue on GitHub:

https://github.com/lixiang117423/biopytools/issues