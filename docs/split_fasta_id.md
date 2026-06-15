# FASTA ID 分割 | FASTA ID Splitter

**从 FASTA 序列头中按分隔符提取指定位置的字段，简化冗长的序列 ID | Extract a specified field from FASTA headers to simplify verbose sequence IDs**

## 功能概述 | Overview

`split_fasta_id` 用于处理来自公共数据库（NCBI、Ensembl、UniProt 等）的 FASTA 文件——这些文件的头行往往包含多个字段（数据库标识、物种、描述等），用空格、制表符或竖线分隔，导致下游工具（比对、注释、组装评估）看到的序列 ID 冗长且不统一。本工具按用户指定的分隔符和位置，从每个头行中只提取一个字段作为新的序列 ID，序列本身保持不变。

分隔符支持自动检测（`auto`）、空格、制表符、同时识别空格和制表符（`both`），或任意单字符（如 `|`、`,`）。NCBI 风格的 `>gi|12345|ref|NC_000001.2|人类...` 可以通过指定 `|` 分隔符和位置 `3` 提取出 `NC_000001.2`。支持保留注释、跳过空头行等选项。典型场景：基因组/蛋白序列入库前的 ID 标准化、多文件合并前的 ID 统一。

## 快速开始 | Quick Start

```bash
# 默认提取第一个字段（空格分隔）
biopytools split-fasta-id -i input.fasta -o output.fasta

# NCBI 风格：用 | 分隔，取第 4 个字段作为 ID
biopytools split-fasta-id -i ncbi.fasta -o clean.fasta -d '|' -p 3
```

## 参数说明 | Parameters

### 必需参数 | Required

| 参数 | 描述 |
|------|------|
| `-i, --input` | 输入 FASTA 文件路径 |

### 常用可选参数 | Common Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-o, --output` | `output.fasta` | 输出 FASTA 文件路径 |
| `-p, --position` | `0` | 提取位置（0 表示第一个元素）|
| `-d, --delimiter` | `auto` | 分隔符：`auto` / `space` / `tab` / `both` / 任意字符（如 `,` 或 `|`）|
| `--keep-original` | 关 | 保留原始文件作为备份 |
| `--no-skip-empty` | 关 | 不跳过空的序列名行（默认会跳过）|
| `--preserve-comments` | 关 | 保留序列名行中的注释 |

（运行 `biopytools split-fasta-id -h` 查看完整参数列表）

## 输出 | Output

```
./
├── output.fasta           # ID 已简化的 FASTA 文件（序列内容不变）
└──（--keep-original 时）
    └── input.fasta.bak    # 原始文件备份
```

示例：`>gi|12345|ref|NC_000001.2| description`（`-d '|' -p 3`）→ `>NC_000001.2`

## 依赖 | Dependencies

- Python 3.7+（纯 Python 实现，无第三方依赖）

## 引用 | Citation

- FASTA 序列格式：Pearson, W. R. & Lipman, D. J. Improved tools for biological sequence comparison. *PNAS* 85, 2444-2448 (1988).

## 相关链接 | References

- [项目主页](https://github.com/lixiang117423/biopytools)
