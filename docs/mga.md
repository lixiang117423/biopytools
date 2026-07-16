# MGA共识基因组组装 | MGA Consensus Genome Assembly

**基于MGA(consensusLJA)对HiFi reads进行共识基因组组装, 支持断点续传与dry-run | Consensus genome assembly from HiFi reads via MGA (consensusLJA), with checkpoint resume and dry-run**

## 功能概述 | Overview

mga 模块封装了 [MGA](https://github.com/SRSRASTA/MGA)(consensusLJA 流程), 用于 PacBio HiFi reads 的共识基因组组装。MGA 二进制本身不在 conda 环境中, 但运行时依赖环境内的 minimap2 / samtools / python(pysam/Bio/numpy), 因此模块用显式 `conda run -n <env> --no-capture-output` 包装调用。

特性:
- **断点续传**: 最终产物 `5_polishing/assembly.fasta` 存在则整体跳过
- **read name 预检**: 扫描 reads 首行, name 含空格则 WARNING(可能导致 MGA/LJA 报错)
- **dry-run**: 只打印命令不执行, 便于拿到完整命令

## 快速开始 | Quick Start

```bash
# 基础用法
biopytools mga -r reads.fastq.gz -o out_dir/

# 指定线程
biopytools mga -r hifi.fq.gz -o out_dir/ -t 64

# 只打印命令不执行
biopytools mga -r hifi.fq.gz -o out_dir/ --dry-run
```

## 参数说明 | Parameters

### 必需参数 | Required

| 参数 | 描述 |
|------|------|
| `-r, --reads` | HiFi reads(fasta/fastq, 可 gz) |
| `-o, --output-dir` | 输出目录 |

### 常用可选参数 | Common Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-t, --threads` | `50` | 线程数 |
| `--mga-path` | 代码默认 | MGA 二进制路径 |
| `--conda-env` | `mga` | conda 环境名(提供 minimap2/samtools/python) |
| `--dry-run` | `False` | 只打印命令不执行 |

(运行 `biopytools mga -h` 查看完整参数列表)

## 输出 | Output

- `5_polishing/assembly.fasta`: 最终共识组装(断点续传判断依据)
- `00_pipeline_info/software_versions.yml`: 软件版本(MGA/minimap2/samtools/python 依赖)与运行参数
- `99_logs/mga.log`: 运行日志

## 断点续传 | Checkpoint Resume

最终产物 `5_polishing/assembly.fasta` 已存在时整体跳过; dry-run 模式不执行命令。

## 路径配置 | Path Configuration

MGA 二进制路径按优先级: 环境变量 `MGA_PATH` > `~/.config/biopytools/config.yml` 的 `tools.mga` > 代码默认值(默认指向 `~/software/MGA/consensusLJA/bin/MGA`)。

## 依赖 | Dependencies

- **MGA / consensusLJA**: 共识组装主程序
- 运行时 conda 环境(默认 `mga`)需含: minimap2、samtools、python + pysam/Bio/numpy

## 相关链接 | References

- [项目主页](https://github.com/lixiang117423/biopytools)
