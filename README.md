# Biohelpers

[![Python Version](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://python.org)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![PyPI Version](https://img.shields.io/badge/pypi-0.2.0-orange.svg)](#)

**生物信息学常用工具包** | **Bioinformatics Common Toolkits**

一个集成了多种生物信息学分析工具的Python包，提供统一的API和命令行接口。

A Python package integrating various bioinformatics analysis tools with unified API and command-line interfaces.

---

## 安装 | Installation

### 从源码安装 | Install from source

```bash
git clone https://github.com/lixiang117423/biopytools.git
cd biohelpers
pip install -e .
```

### 安装开发依赖 | Install development dependencies

```bash
pip install -e ".[dev]"
```

---

## Module列表 | Module List

### run_fastp

用于处理双端测序数据的批量质控和报告生成。支持生成HTML和JSON格式的质控报告。

点击[这里](https://lixiang117423.github.io/article/run-fastp/)查看详细文档。

---

## 许可证 | License

本项目使用MIT许可证 - 查看 [LICENSE](LICENSE) 文件了解详情

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## 支持与反馈 | Support & Feedback

- 📧 Email: biohelpers@example.com
- 🐛 Issues: [GitHub Issues](https://github.com/yourusername/biohelpers/issues)
- 📖 Documentation: [Read the Docs](https://biohelpers.readthedocs.io)

---

# 更新日志 | Changelog

## v1.1.0 (2025-07-15)
- 增加run_fastp模块
