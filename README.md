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

[Admixture | 群体结构](./docs/admixture.md) - 🧬 ADMIXTURE群体结构分析工具

[OrthoFinder | 基于基因的泛基因组构建](./docs/orthofinder.md) - 🔍 泛基因组分析和系统发育构建

[VCF系统发育 | VCF系统发育分析](./docs/vcf_phylo.md) - 🌳 基于VCF文件的系统发育分析

[VCF转PHYLIP | VCF格式转换工具](./docs/vcf2phylip.md) - 🔄 VCF到PHYLIP格式转换

### 🧪 序列分析与比对 | Sequence Analysis & Alignment

[Blast | Blast比对](./docs/blast.md) - 🔍 BLAST序列比对分析工具

[Longest mRNA | 提取最长转录本](./docs/longest_mrna.md) - 📜 从转录组数据中提取最长转录本

[Subseq | 序列子集提取工具](./docs/subseq.md) - ✂️ 序列子集提取和分析工具

### 📝 变异检测与注释 | Variant Detection & Annotation

[Annovar | 变异注释](./docs/annovar.md) - 📝 ANNOVAR变异功能注释工具

### 📥 数据获取与处理 | Data Acquisition & Processing

[SRA2FASTQ | SRA转FASTQ](./docs/sra2fastq.md) - 🧬 SRA数据库数据转FASTQ格式

[从CNCB获取链接 | 批量获取测序数据下载链接](./docs/get_link_from_CNCB.md) - 📥 CNCB数据库批量数据下载工具

### 🧬 基因组组装与挂载 | Genome Assembly & Scaffolding

[HapHiC | 基因组scaffolding工具](./docs/haphic.md) - 🧬 基于Hi-C数据的基因组scaffolding工具


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