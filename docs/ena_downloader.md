# ENA数据下载工具 | ENA Data Downloader

**从ENA数据库批量下载FASTQ测序数据及其元数据 | Batch download FASTQ sequencing data and metadata from the ENA database**

## 功能概述 | Overview

基于ENA Portal API，根据项目编号（如 `PRJNA`、`PRJEB`、`ERP`、`SRP` 等）获取样本元数据，并生成或直接执行FASTQ下载脚本。支持FTP（wget）与Aspera两种下载协议，支持仅下载元数据模式，方便下游分析整理。

- 自动调用ENA filereport API拉取丰富的样本/Run元数据
- 生成可直接运行的bash下载脚本（默认）或立即下载
- FTP与Aspera双协议，Aspera可显著提升大批量下载速度
- 元数据格式可选 TSV / CSV / XLSX
- 自动生成下载汇总报告与运行日志

## 快速开始 | Quick Start

```bash
# 最简：根据accession在当前目录生成wget下载脚本
biopytools ena-downloader -a PRJEB12345

# 在指定目录下创建专用结果目录并生成脚本
biopytools ena-downloader -a PRJNA661210 -o ./results -d

# 仅下载元数据（不生成FASTQ下载脚本）
biopytools ena-downloader -a PRJEB12345 -M -f xlsx

# 使用Aspera高速下载（需提供密钥，权限需为600）
biopytools ena-downloader -a PRJEB12345 -p aspera -k ~/.aspera/asperaweb_id_dsa.openssh

# 直接执行下载（而不是只生成脚本）
biopytools ena-downloader -a PRJEB12345 -m run
```

## 参数说明 | Parameters

### 必需参数 | Required

| 参数 | 描述 |
|------|------|
| `-a, --accession` | ENA项目/研究编号，例如 `PRJEB12345`、`PRJNA661210`、`SRP123456` |

### 常用可选参数 | Common Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-o, --output-dir` | 当前目录 | 输出目录 |
| `-d, --create-dir` | False | 创建专用输出目录（未指定 `-o` 时默认为 `<accession>.ena.download`） |
| `-f, --metadata-format` | `tsv` | 元数据文件格式，可选 `tsv` / `csv` / `xlsx` |
| `-p, --protocol` | `ftp` | 下载协议，可选 `ftp`（使用 wget）或 `aspera`（使用 ascp） |
| `-k, --aspera-key` | 无 | Aspera私钥路径，`-p aspera` 时必需，权限必须为 `600` |
| `-m, --method` | `save` | 执行模式：`save` 只生成下载脚本；`run` 直接执行下载 |
| `-M, --metadata-only` | False | 仅下载元数据，不处理FASTQ下载 |
| `-F, --fields` | 默认字段集 | 自定义元数据字段，可多次指定或空格分隔（如 `-F fastq_ftp tax_id scientific_name`，传入 `all` 使用全部字段） |
| `-r, --max-retries` | `3` | ENA API请求最大重试次数 |

（运行 `biopytools ena-downloader -h` 查看完整参数列表）

## 输出 | Output

默认（`-m save`）将在输出目录生成：

- `<accession>.meta.<tsv|csv|xlsx>`：样本/Run元数据表
- `download_<accession>_fastq_by_wget.sh` 或 `download_<accession>_fastq_by_aspera.sh`：FASTQ下载脚本
- `download_summary.txt`：下载汇总报告（含项目信息、FASTQ文件数量、下一步操作提示）
- `ena_download.log`：运行日志

随后执行 `bash download_<accession>_fastq_by_*.sh` 即可开始实际下载。

## 依赖 | Dependencies

- Python 包：`requests`、`openpyxl`（生成xlsx时）
- FTP协议：系统需安装 `wget`
- Aspera协议：系统需安装 IBM Aspera Connect（`ascp`），并提供私钥文件（权限 `600`）

## 引用 | Citation

欧洲核酸数据库（ENA）属于 INSDC。使用其数据时请引用：
- Cummins C et al. The European Nucleotide Archive in 2021. *Nucleic Acids Research*. 2022, D1.
- ENA Portal API: https://www.ebi.ac.uk/ena/portal/api/

## 相关链接 | References

- [ENA Browser](https://www.ebi.ac.uk/ena/browser/)
- [ENA Portal API 文档](https://docs.ena-browser.org/#/api/filereport)
- [项目主页](https://github.com/lixiang117423/biopytools)
