# 序列提取 | Sequence Extraction

**seqkit 封装的序列提取, 自动识别单个ID / ID文件 / BED区间 | seqkit-based sequence extraction, auto-detects single ID / ID file / BED regions**

## 功能概述 | Overview

seq_extract 模块封装了 [seqkit](https://bioinf.shenwei.me/seqkit/), 提供统一的序列提取入口。根据查询输入自动选择提取方式:

- **单个 ID**: 不是文件时, 用 `seqkit grep -p <id>` 按序列名提取
- **ID 文件**(一列): 文件首行为单列时, 用 `seqkit grep -f <id_file>` 批量提取
- **BED 文件**(≥2 列): 文件首行 tab 分隔 ≥2 列时, 用 `seqkit subseq --bed <bed>` 按区间提取

输出文件名未指定时自动推导为 `{query}.{subject}.fa`。

## 快速开始 | Quick Start

```bash
# 单个 ID 提取
biopytools seq-extract -i chr1 -s genome.fa -o chr1.fa

# ID 文件批量提取
biopytools seq-extract -i gene.id.txt -s gene.fa -o gene.genomic.fa

# BED 区间提取
biopytools seq-extract -i regions.bed -s genome.fa -o regions.fa
```

## 参数说明 | Parameters

### 必需参数 | Required

| 参数 | 描述 |
|------|------|
| `-i, --input` | 查询: 单个 ID、ID 文件(一列)或 BED 文件(≥2 列) |
| `-s, --sequence` | 目标序列 FASTA 文件 |

### 常用可选参数 | Common Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-o, --output` | 自动推导 | 输出文件(默认 `{query}.{subject}.fa`) |
| `--bed` | `False` | 强制 BED 模式(跳过自动检测) |

(运行 `biopytools seq-extract -h` 查看完整参数列表)

## 输出 | Output

- 输出 FASTA 文件(路径由 `-o` 指定或自动推导)

## 自动检测规则 | Auto-detection

| 输入形态 | 判定 | 底层命令 |
|----------|------|----------|
| 非文件路径 | 单个 ID | `seqkit grep -p <id>` |
| 文件, 首行单列 | ID 文件 | `seqkit grep -f <id_file>` |
| 文件, 首行 ≥2 列(tab) | BED | `seqkit subseq --bed <bed>` |

## 路径配置 | Path Configuration

seqkit 路径按优先级: 环境变量 `SEQKIT_PATH` > `~/.config/biopytools/config.yml` 的 `tools.seqkit` > 代码默认值(默认 `seqkit`, 即依赖 PATH)。

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 查询:单个ID、ID文件(一列)或BED文件(>=2列)｜Query: single ID, ID file (1 column), or BED file (>=2 columns) |
| `-s, --sequence` | 必填 |  | 目标序列FASTA文件｜Target sequence FASTA file |
| `-o, --output` | — |  | 输出文件(默认自动推导:{query}.{subject}.fa)｜Output file (default: auto-derived) |
| `--bed` | — |  | 强制BED模式(跳过自动检测)｜Force BED mode (skip auto-detection) |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 查询:单个ID、ID文件(一列)或BED文件(>=2列)｜Query: single ID, ID file (1 column), or BED file (>=2 columns) |
| `-s, --sequence` | 必填 |  | 目标序列FASTA文件｜Target sequence FASTA file |
| `-o, --output` | — |  | 输出文件(默认自动推导:{query}.{subject}.fa)｜Output file (default: auto-derived {query}.{subject}.fa) |
| `--bed` | — | store_true | 强制BED模式(跳过自动检测)｜Force BED mode (skip auto-detection) |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- **seqkit**: FASTA/FASTQ 序列处理 (https://bioinf.shenwei.me/seqkit/)

## 相关链接 | References

- [项目主页](https://github.com/lixiang117423/biopytools)
