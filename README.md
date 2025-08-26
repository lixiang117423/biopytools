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

## 安装方法 | Installation


### 从源码安装 | Install from source

```bash
git clone https://github.com/lixiang117423/biopytools.git
cd biopytools
pip install -e .
```

## 使用方法 | Useage

```bash
biopytools -h           
Usage: biopytools [OPTIONS] COMMAND [ARGS]...

  BioPyTools - 生物信息学分析工具包

  要查看特定命令的帮助，请运行：biopytools <命令> -h/--help, 如biopytools fastp -h

Options:
  -v, --version  Show the version and exit.
  -h, --help     Show this message and exit.

Commands:
  admixture        ADMIXTURE群体结构分析
  annovar          ANNOVAR变异注释工具
  bismark          Bismark甲基化分析流程
  blast            BLAST序列比对分析工具
  coverage         BAM/SAM文件覆盖度分析工具
  ena-downloader   从ENA下载测序数据样品信息和下载链接
  fastp            批量运行fastp
  geneinfo         从GFF文件中提取基因和转录本的信息
  genomesyn        基因组共线性分析
  gtx              运行GTX WGS流程
  hifiasm          运行hifiasm组装流程
  kaks             计算Ks,Ks及其比值
  kmer-count       计算K-mer的丰度和滑窗分析
  kmer-query       从FASTA/FASTQ文件中提取K-mer序列
  longest-mrna     从基因组中提取每个基因的最长转录本
  minimap2         Minimap2全基因组比对与未比对区域提取
  mrna-prediction  转录组预测分析工具
  parabricks       基于GPU加速的全基因组测序(WGS)批处理分析工具
  plink-gwas       PLINK全基因组关联分析
  popgen           群体遗传学多样性分析
  raxml            RAxML系统发育分析工具
  rnaseq           RNA-seq表达定量分析流程
  split-fasta-id   FASTA序列ID分割和清理工具
  vcf-filter       VCF文件高性能筛选和格式转换
  vcf-genotype     VCF文件基因型数据提取和格式转换
  vcf-nj-tree      VCF文件系统发育树构建和分析
  vcf-pca          VCF文件主成分分析(PCA)工具
  vcf-sample-hete  VCF基因型统计分析工具
  vcf-sequence     VCF和基因组序列变异提取工具
  vcf2phylip       VCF格式转换工具
```

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
