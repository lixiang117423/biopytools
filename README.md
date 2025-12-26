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

## 模块文档 | README

### 🧬 全流程分析 | End-to-End Analysis
[Fastq到VCF (Parabricks) | GPU加速全流程分析](./docs/fastq2vcf_parabricks.md) - 🚀 基于NVIDIA Parabricks的GPU加速FASTQ到VCF完整分析流程
[Fastq到VCF (GTX) | CPU优化全流程分析](./docs/fastq2vcf_gtx.md) - 🔬 基于GTX的CPU优化FASTQ到VCF完整分析流程
[Parabricks | GPU加速WGS分析](./docs/parabricks.md) - 🚀 单独的Parabricks GPU加速WGS分析工具

### 📊 群体遗传与进化 | Population Genetics & Evolution
[GEMMA GWAS | 全基因组关联分析](./docs/gemma-gwas.md) - 🧬 基于GEMMA的混合线性模型GWAS分析工具
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
[EGAPx Batch | EGAPx批量配置生成](./docs/egapx-batch.md) - 🧬 按染色体批量生成EGAPx运行配置和脚本
[Genome Analysis | 基因组特征评估](./docs/genome_analysis.md) - 🧬 基于GenomeScope2和Smudgeplot的基因组分析工具
[Rename Chromosomes | 染色体重命名](./docs/rename_chromosomes.md) - 🧬 FASTA序列标准化重命名工具
[Dsuite | D统计量分析](./docs/dsuite.md) - 🧬 基因渗入检测和D统计量计算工具


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