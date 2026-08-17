# iSeq 公共测序数据下载工具

**专业的公共数据库测序数据下载工具 | Professional Public Sequencing Data Download Tool**

## 📖 功能概述 | Overview

iSeq下载工具基于iSeq软件构建，提供从GSA、SRA、ENA、DDBJ等公共数据库下载测序数据和元数据的完整流程。支持多种下载方式（包括Aspera高速下载）、自动格式转换、智能文件合并等功能，适用于各种基因组学研究的数据获取需求。

## ✨ 主要特性 | Key Features

- **🌐 多数据库支持**: GSA、SRA、ENA、DDBJ数据库一站式下载
- **⚡ 高速下载**: 支持Aspera协议，大幅提升下载速度
- **🔄 智能格式处理**: 自动下载gzip格式或转换SRA为FASTQ
- **📦 批量合并**: 支持按实验/样本/研究合并FASTQ文件
- **⚙️ 灵活配置**: 支持多种下载选项和性能参数调整
- **📊 元数据获取**: 可仅下载元数据而不下载测序数据
- **🚀 并行下载**: 支持多线程并行连接，提升下载效率
- **📝 详细日志**: 完整的下载过程日志和汇总报告

## 🚀 快速开始 | Quick Start

### 基本用法 | Basic Usage

```bash
# 使用Aspera高速下载单细胞项目数据（推荐）
biopytools iseq \
    -i PRJNA1014406 \
    -a \
    -g \
    -p 10 \
    -t 16 \
    -o ./scRNA_data
```

### 高级用法 | Advanced Usage

```bash
# 仅获取项目元数据
biopytools iseq -i PRJNA1014406 -m -o ./metadata

# 合并实验级别的FASTQ文件
biopytools iseq -i PRJNA1014406 -g -e ex -o ./merged_data

# 从SRA数据库下载（不使用Aspera）
biopytools iseq -i SRR123456 -g -d sra -o ./sra_data
```

## 📋 参数说明 | Parameters

### 必需参数 | Required Parameters

| 参数 | 描述 | 示例 |
|------|------|------|
| `-i, --accession` | 项目/样本/实验ID | `-i PRJNA1014406` |

### 路径配置 | Path Configuration

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--iseq-path` | `/share/org/YZWL/yzwl_lixg/miniforge3/envs/iseq_v.1.9.8/bin/iseq` | 🛠️ iSeq软件路径 |
| `-c, --conda-env` | `iseq_v.1.9.8` | 🐍 Conda环境名 |
| `-o, --output-dir` | `./iseq_output` | 📁 输出目录路径 |

### 下载选项 | Download Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-m, --metadata-only` | `False` | 📊 仅下载元数据 |
| `-g, --gzip` | `True` | 📦 下载gzip格式FASTQ |
| `-q, --fastq` | `False` | 🔧 转换为FASTQ格式 |
| `-e, --merge` | `None` | 🔗 合并选项 (ex/sa/st) |

### 性能参数 | Performance Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-t, --threads` | `16` | ⚙️ 线程数（用于转换和压缩） |
| `-p, --parallel` | `10` | 🔗 并行连接数 |
| `-s, --speed` | `None` | 🚄 下载速度限制 (MB/s) |

### 数据库选项 | Database Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-d, --database` | `ena` | 💾 数据库选择 (ena/sra) |
| `--protocol` | `ftp` | 🌐 协议选择 (ftp/https) |

### 高级选项 | Advanced Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-a, --use-aspera` | `False` | ⚡ 使用Aspera高速下载 |
| `--skip-md5` | `False` | ⏭️ 跳过MD5校验 |
| `--quiet` | `False` | 🔇 静默模式（不显示进度条） |

## 📚 Accession格式说明 | Accession Format

支持的accession前缀：

| 类型 | 前缀 | 示例 |
|------|------|------|
| **项目** | PRJEB, PRJNA, PRJDB, PRJC, GSE | PRJNA1014406 |
| **研究** | ERP, DRP, SRP, CRA | ERP123456 |
| **生物样本** | SAMD, SAME, SAMN, SAMC | SAMN123456 |
| **样本** | ERS, DRS, SRS, GSM | ERS123456 |
| **实验** | ERX, DRX, SRX, CRX | ERX123456 |
| **运行** | ERR, DRR, SRR, CRR | SRR123456 |

## 💡 使用示例 | Usage Examples

### 示例1：高速下载单细胞数据 | Example 1: High-speed single-cell data download

```bash
# 使用Aspera + 多线程下载
biopytools iseq \
    -i PRJNA1014406 \
    -a \
    -g \
    -p 10 \
    -t 16 \
    -o ./scRNA_fastq
```

### 示例2：仅获取元数据 | Example 2: Metadata only

```bash
# 下载项目元数据信息
biopytools iseq \
    -i PRJNA1014406 \
    -m \
    -o ./project_metadata
```

### 示例3：合并FASTQ文件 | Example 3: Merge FASTQ files

```bash
# 按实验级别合并（适用于有多个run的样本）
biopytools iseq \
    -i PRJNA1014406 \
    -g \
    -e ex \
    -o ./merged_by_experiment
```

### 示例4：从SRA数据库下载 | Example 4: Download from SRA database

```bash
# 指定从SRA数据库下载
biopytools iseq \
    -i SRR12345678 \
    -g \
    -d sra \
    -t 24 \
    -o ./sra_download
```

### 示例5：限制下载速度 | Example 5: Limit download speed

```bash
# 限制下载速度为50MB/s
biopytools iseq \
    -i PRJNA1014406 \
    -g \
    -a \
    -s 50 \
    -o ./data
```

## 📊 输出结果 | Output Results

### 输出目录结构 | Output Directory Structure

```
iseq_output/
├── PRJNA1014406/              # 项目目录|Project directory
│   ├── SRR123456_1.fastq.gz  # FASTQ文件|FASTQ files
│   ├── SRR123456_2.fastq.gz
│   └── ...
├── metadata.txt              # 元数据文件|Metadata file
├── download_summary.txt      # 下载汇总报告|Download summary report
└── iseq_download.log         # 运行日志|Run log
```

## 🔧 系统要求 | System Requirements

### 依赖软件 | Dependencies

- **iSeq** (版本 1.9.8 或更新)
  - 安装环境: `iseq_v.1.9.8`
  - 软件路径: `/share/org/YZWL/yzwl_lixg/miniforge3/envs/iseq_v.1.9.8/bin/iseq`

- **Python** (版本 3.7+)

- **iSeq依赖软件**:
  - [pigz](https://github.com/madler/pigz) (>=2.8) - 多线程压缩
  - [wget](https://www.gnu.org/software/wget/) (>=1.16) - 文件下载
  - [axel](https://github.com/axel-download-accelerator/axel) (>=2.17) - 多线程下载
  - [aspera-cli](https://github.com/IBM/aspera-cli) (=4.14.0) - Aspera高速下载
  - [sra-tools](https://github.com/ncbi/sra-tools) (>=2.11.0) - SRA工具

### 环境配置 | Environment Setup

```bash
# iSeq已安装在独立conda环境中|iSeq is installed in separate conda environment
conda activate iseq_v.1.9.8

# 验证依赖软件|Verify dependencies
pigz --version
wget --version
axel --version
ascp --version
srapath --version

# 验证iSeq版本|Verify iSeq version
iseq --version
```

## ⚠️ 注意事项 | Important Notes

1. **网络连接**: 确保服务器能够访问外部数据库的网络连接
2. **存储空间**: 测序数据通常很大，确保有足够的磁盘空间
3. **Aspera密钥**: 使用Aspera下载需要正确配置密钥文件
4. **并发限制**: 过高的并发数可能导致服务器限制连接
5. **数据完整性**: 下载完成后建议进行MD5校验

## 🐛 故障排除 | Troubleshooting

### 常见问题 | Common Issues

**Q: "srapath: command not found" 错误**

```bash
# 激活conda环境
conda activate iseq_v.1.9.8

# 或重新安装sra-tools
conda install -c bioconda sra-tools
```

**Q: Aspera下载失败**

```bash
# 检查Aspera安装
ascp --version

# 检查密钥文件
ls -la ~/.aspera/connect/etc/asperaweb_id_dsa.openssh

# 如Aspera不可用，可以不使用-a参数
biopytools iseq -i PRJNA1014406 -g -p 10 -t 16 -o ./data
```

**Q: 下载速度慢**

```bash
# 增加并行连接数
biopytools iseq ... -p 20

# 或使用Aspera
biopytools iseq ... -a
```

**Q: 磁盘空间不足**

```bash
# 仅下载元数据查看项目大小
biopytools iseq -i PRJNA1014406 -m -o ./check
```

## 📚 相关资源 | Related Resources

- [iSeq官方文档](https://github.com/BioOmics/iSeq)
- [iSeq中文教程](https://github.com/BioOmics/iSeq/blob/main/docs/ChineseTutorial.md)
- [iSeq使用示例](https://github.com/BioOmics/iSeq/blob/main/docs/Examples.md)
- [ENA数据库](https://www.ebi.ac.uk/ena)
- [SRA数据库](https://www.ncbi.nlm.nih.gov/sra)
- [GSA数据库](https://ngdc.cncb.ac.cn/gsa/)

## 📄 许可证 | License

本项目采用MIT许可证 - 详见 [LICENSE](LICENSE) 文件

iSeq软件本身遵循其原始许可证（GPL-3.0）。

---

## 🔬 引用信息 | Citation

如果在学术研究中使用iSeq工具，请引用原始文献：

```
Chao H, Li Z, Chen D, Chen M.
iSeq: An integrated tool to fetch public sequencing data.
Bioinformatics, 2024, btae641.
doi: 10.1093/bioinformatics/btae641
```

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --accession` | 必填 |  | 项目/样本/实验ID｜Project/Sample/Experiment ID (e.g., PRJNA1014406) |
| `-o, --output-dir` | `./iseq_output` | Path | 输出目录｜Output directory |
| `--iseq-path` | `~/miniforge3/envs/misc/bin/iseq` |  | iSeq软件路径｜iSeq software path |
| `-c, --conda-env` | `iseq_v.1.9.8` |  | Conda环境名｜Conda environment name |
| `-m, --metadata-only` | — |  | 仅下载元数据｜Only download metadata |
| `-g, --gzip` | — |  | 下载gzip格式FASTQ｜Download FASTQ in gzip format |
| `-q, --fastq` | — |  | 转换为FASTQ格式｜Convert to FASTQ format |
| `-e, --merge` | — | ex/sa/st | 合并选项｜Merge option (ex/sa/st) |
| `-t, --threads` | `16` | int | 线程数｜Number of threads |
| `-p, --parallel` | `10` | int | 并行连接数｜Number of parallel connections |
| `-s, --speed` | — | int | 下载速度限制(MB/s)｜Download speed limit (MB/s) |
| `-d, --database` | `ena` | ena/sra | 数据库选择｜Database selection |
| `--protocol` | `ftp` | ftp/https | 协议选择｜Protocol selection |
| `-a, --use-aspera` | — |  | 使用Aspera高速下载｜Use Aspera for high-speed download |
| `--skip-md5` | — |  | 跳过MD5校验｜Skip MD5 check |
| `--quiet` | — |  | 静默模式｜Quiet mode |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --accession` | 必填 |  | 项目/样本/实验ID｜Project/Sample/Experiment ID (e.g., PRJNA1014406) |
| `-o, --output-dir` | `./iseq_output` |  | 输出目录｜Output directory |
| `-a, --iseq-path` | `~/miniforge3/envs/misc/bin/iseq` |  | iSeq软件路径｜iSeq software path |
| `-m, --metadata-only` | — | store_true | 仅下载元数据｜Only download metadata |
| `-g, --gzip` | — | store_true | 下载gzip格式FASTQ｜Download FASTQ in gzip format |
| `-q, --fastq` | — | store_true | 转换为FASTQ格式｜Convert to FASTQ format |
| `-e, --merge` | — | ex/sa/st | 合并选项｜Merge option (ex/sa/st) |
| `-t, --threads` | `16` | int | 线程数｜Number of threads |
| `-p, --parallel` | `10` | int | 并行连接数｜Number of parallel connections |
| `-s, --speed` | — | int | 下载速度限制(MB/s)｜Download speed limit (MB/s) |
| `-d, --database` | `ena` | ena/sra | 数据库选择｜Database selection |
| `--protocol` | `ftp` | ftp/https | 协议选择｜Protocol selection |
| `--use-aspera` | — | store_true | 使用Aspera下载｜Use Aspera for download |
| `--skip-md5` | — | store_true | 跳过MD5校验｜Skip MD5 check |
| `--quiet` | — | store_true | 静默模式｜Quiet mode (suppress progress bars) |

<!-- END PARAMS:auto -->
