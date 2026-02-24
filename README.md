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

## 环境配置 / Environment Setup

Conda 环境配置文件位于 [`conda_env/`](conda_env/) 目录下。

Conda environment files can be found in the [`conda_env/`](conda_env/) directory.

## 安装方法 | Installation

### 从源码安装 | Install from source

```bash
git clone https://github.com/lixiang117423/biopytools.git
cd biopytools
pip install -e .

# or
pip install .
```

## 使用方法 | Usage

### 查看帮助 | Getting Help

```bash
biopytools -h
Usage: biopytools [OPTIONS] COMMAND [ARGS]...

  BioPyTools - 生物信息学分析工具包

  要查看特定命令的帮助，请运行：biopytools <命令> -h/--help, 如biopytools annovar -h

Options:
  -v, --version  Show the version and exit.
  -h, --help     Show this message and exit.
```

## 模块文档 | README

[admixture](./docs/admixture.md) - [群体结构分析软件Admixture](https://genome.cshlp.org/content/19/9/1655)

[annovar](./docs/annovar.md) - [ANNOVAR变异功能注释工具](https://academic.oup.com/nar/article/38/16/e164/1749458)

[bam-cov](./docs/bam_coverage_stats.md) - BAM文件覆盖度统计

[blast](./docs/blast_v2.md) - [序列比对工具](https://academic.oup.com/nar/article/36/suppl_2/W5/2505810?login=true)

[busco](./docs/busco.md) - [BUSCO](https://academic.oup.com/bioinformatics/article/31/19/3210/211866)

[fastp](./docs/fastp.md) - [fastq文件质控](https://onlinelibrary.wiley.com/doi/10.1002/imt2.70078)

[iseq](./docs/iseq.md) - [iSeq下载测序数据](https://github.com/BioOmics/iSeq)

[rna-seq](./docs/rnaseq.md) - [Hisat2](https://www.nature.com/articles/s41587-019-0201-4) + [StringTie2](https://link.springer.com/article/10.1186/s13059-019-1910-1)转录组流程

[sra2fastq](docs/sra2fastq.md) - SRA文件转FSATQ文件-基于[parallel-fastq-dump](https://github.com/rvalieris/parallel-fastq-dump)


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