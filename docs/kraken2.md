# Kraken2宏基因组分类 | Kraken2 Metagenomic Classification

**基于Kraken2+Bracken的宏基因组/扩增子物种分类与丰度估计 | Metagenomic taxonomic classification and abundance estimation with Kraken2 + Bracken**

## 功能概述 | Overview

输入一个含配对 FASTQ（默认后缀 `_1.clean.fq.gz` / `_2.clean.fq.gz`）的目录与一个 Kraken2 数据库目录，自动批量对每个样品运行 Kraken2 物种分类，并可选用 Bracken 在指定分类级别（默认种水平 S）重新估计丰度。流程支持断点续跑、记录软件版本信息，并自动汇总分类报告。

- 自动按 R1/R2 后缀匹配样品，断点续跑已完成的步骤
- Kraken2 分类 + Bracken 丰度估计一体化，可单独关闭 Bracken
- 可配置 Kraken2 置信度阈值与 Bracken 分类级别、最小读数阈值
- 自动生成 `software_versions.yml` 与项目级日志

## 快速开始 | Quick Start

```bash
biopytools kraken2 -i ./fastq/ -d ~/database/kraken2_db -o ./kraken2_output

# 关闭 Bracken，仅做 Kraken2 分类
biopytools kraken2 -i ./fastq/ -d ~/db/kraken2_db --no-bracken

# 自定义 Bracken 分类级别（属 G）与置信度
biopytools kraken2 -i ./fastq/ -d ~/db/kraken2_db --bracken-level G --confidence 0.1
```

## 参数说明 | Parameters

### 必需参数 | Required

| 参数 | 描述 |
|------|------|
| `-i, --input-dir` | 输入 FASTQ 目录（含配对 R1/R2） |
| `-d, --db-path` | Kraken2 数据库目录（应含已构建好的 `hash.k2d` 等文件，Bracken 需 `database*.kmer_distrib`） |

### 常用可选参数 | Common Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-o, --output-dir` | `./kraken2_output` | 输出目录 |
| `-t, --threads` | `12` | 线程数 |
| `--read-len` | `150` | 读长（Bracken 用） |
| `--confidence` | `0.0` | Kraken2 置信度阈值（0-1） |
| `--bracken-level` | `S` | Bracken 分类级别：`D` / `P` / `C` / `O` / `F` / `G` / `S` / `S1` |
| `--bracken-threshold` | `10` | Bracken 最小读数阈值 |
| `--no-bracken` | False | 跳过 Bracken 分析 |
| `--r1-suffix` | `_1.clean.fq.gz` | R1 文件后缀 |
| `--r2-suffix` | `_2.clean.fq.gz` | R2 文件后缀 |

（运行 `biopytools kraken2 -h` 查看完整参数列表）

## 输出 | Output

输出目录典型结构：

- `01_kraken2/`：每样品的 Kraken2 结果
  - `<sample>.kraken`：每条 read 的分类输出
  - `<sample>.kraken_report.txt`：Kraken2 分类报告（RA 嵌套计数）
- `02_bracken/`：Bracken 丰度估计
  - `<sample>.bracken.txt`：Bracken 修正后的丰度表
  - `<sample>.bracken_report.txt`：Bracken 修正后的报告
- `03_summary/`：样品级汇总
- `00_pipeline_info/software_versions.yml`：软件版本
- `99_logs/kraken2_pipeline.log`：运行日志

## 依赖 | Dependencies

- `kraken2`（物种分类主程序）
- `bracken`（丰度估计，默认启用）
- 已构建的 Kraken2 数据库（含 Bracken 的 `database*.kmer_distrib` 文件）

推荐通过 `conda install -c bioconda kraken2 bracken` 安装；数据库构建见 Kraken2 官方文档（`kraken2-build`）。

## 引用 | Citation

- Wood DE, Lu J, Langmead B. Improved metagenomic analysis with Kraken 2. *Genome Biology*. 2019, 20(1):257.
- Lu J, Breitwieser FP, Thielen P, Salzberg SL. Bracken: estimating species abundance in metagenomics data. *PeerJ Computer Science*. 2017, 3:e104.

## 相关链接 | References

- [Kraken2 项目](https://github.com/DerrickWood/kraken2)
- [Bracken 项目](https://github.com/jenniferlu717/bracken)
- [项目主页](https://github.com/lixiang117423/biopytools)
