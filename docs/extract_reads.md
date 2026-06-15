# 基于 contig-reads 对应关系提取 FASTQ reads | Extract FASTQ Reads by Contig-Reads Mapping

**根据 contig-reads 对应关系表，从 FASTQ 文件中批量提取指定 reads | Extract specified reads from FASTQ files according to a contig-reads mapping table.**

## 功能概述 | Overview

- 输入 TSV 格式的 contig-reads 对应关系（通常来自组装软件的比对结果）
- 支持从单个或一对压缩/未压缩 FASTQ 中筛选目标 reads
- 支持自动 gzip 压缩输出（默认开启，可用 `--no-compress` 关闭）
- 适合用于从原始测序数据中回溯特定 contig 的测序 reads

## 快速开始 | Quick Start

```bash
# 默认按输出后缀自动判断是否压缩
biopytools extract-reads -m contig_reads.tsv -i input.fq.gz -o output.fq.gz

# 不压缩输出
biopytools extract-reads -m contig_reads.tsv -i input.fq -o output.fq --no-compress
```

## 参数说明 | Parameters

### 必需参数 | Required

| 参数 | 描述 |
|------|------|
| `-m, --mapping` | contig-reads 对应关系文件 (TSV 格式) |
| `-i, --input` | 输入 FASTQ 文件 (支持 gzip 压缩) |
| `-o, --output` | 输出文件路径 |

### 常用可选参数 | Common Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--no-compress` | off | 不压缩输出文件（默认当 `.gz` 结尾时自动 gzip 压缩） |

（运行 `biopytools extract-reads -h` 查看完整参数列表）

## 输入文件格式 | Input Mapping Format

`--mapping` 为 TSV 文件，每行至少包含 contig 与 read 名称。程序会同时建立 `contig -> reads` 与 `read -> contig` 两种映射，最终把所有出现过的 read 名称作为目标集合。

## 输出 | Output

- 单个 FASTQ 文件，包含所有匹配到的 reads
- 若输出路径以 `.gz` 结尾且未指定 `--no-compress`，将自动 gzip 压缩
- 日志中会输出 "需要提取的 reads 总数" 与 "提取 reads 数" 等统计

## 依赖 | Dependencies

- Python 标准库（gzip、os），无第三方依赖

## 引用 | Citation

- 本工具为组装数据回溯工具，无直接引用文献。

## 相关链接 | References

- [项目主页](https://github.com/lixiang117423/biopytools)
