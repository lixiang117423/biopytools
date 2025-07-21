# Biopytools

[![Python Version](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://python.org)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![PyPI Version](https://img.shields.io/badge/pypi-0.2.0-orange.svg)](#)

**生物信息学常用工具包** | **Bioinformatics Common Toolkits**

一个集成了多种生物信息学分析工具的Python包，提供统一的API和命令行接口。

A Python package integrating various bioinformatics analysis tools with unified API and command-line interfaces.

---

## 安装 | Installation

### 从源码安装 | Install from source

建议使用conda进行安装。

```bash
mamba create -n biopytools python=3.10 
mamba activate biopytools

# 先用conda安装rocksdb
mamba install conda-forge::python-rocksdb

git clone https://github.com/lixiang117423/biopytools.git
cd biopytools
pip install -e .
```

### 安装开发依赖 | Install development dependencies

```bash
pip install -e ".[dev]"
```

---

## 功能列表 | Module List

点击[这里](https://lixiang117423.github.io/article/biopytools-readme/)查看文档。

---

## 许可证 | License

本项目使用MIT许可证 - 查看 [LICENSE](LICENSE) 文件了解详情

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## 支持与反馈 | Support & Feedback

- 📧 Email: lixiang117423@gmail.com
- 🐛 Issues: [GitHub Issues](https://github.com/lixiang117423/biopytools/issues)

---

# 更新日志 | Changelog

- v1.13.1
    - 修复parse_sequence_vcf和run_minimap2模块的序列坐标问题

- v1.13.0
    - 增加parse_sequence_vcf模块

- v1.12.0
    - 增加run_repeat_masker模块

- v1.11.0
    - 增加run_minimap2模块

- v1.10.0
    - 增加run_kmer_pav模块

- v1.9.0
    - 增加run_admixture模块

- v1.8.0
    - 增加run_augustus_multi_rnaseq模块

- v1.7.0
    - 增加parse_sample_hete模块

- v1.6.0
    - 增加parse_longest_mrna模块

- v1.5.0
    - 增加run_kmer_analysis模块

- v1.4.0
    - 增加run_plink_gwas模块

- v1.3.0
    - 增加parse_gene_info模块
    - 增加run_annovar模块
    - 增加run_vcf_extractor模块

- v1.1.0
    - 增加run_fastp模块
    - 增加run_rnaseq模块
