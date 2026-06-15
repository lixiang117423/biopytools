# 启动子提取 | Promoter Extractor

**从基因组 FASTA 和 GFF3 注释中批量提取基因上游启动子序列 | Batch-extract upstream promoter sequences from a genome FASTA and GFF3 annotation**

## 功能概述 | Overview

`promoter_extractor` 根据基因组组装（FASTA）和注释（GFF3）文件，为每个 `gene` 特征提取其转录起始位点（TSS）上游指定长度的启动子序列。对正链基因，取 gene 起始位置上游 N bp；对负链基因，取 gene 终止位置下游 N bp 并反向互补，从而保证序列方向与基因转录方向一致。

工具支持自定义启动子长度（默认 2000 bp）、最小可接受长度（当基因距 contig 末端不足时，允许截短输出还是丢弃）、指定基因子集列表（`--gene-list`，只提取感兴趣的基因）、多线程加速以及模拟运行。同时输出 BED（基因组浏览器可视化）和统计报告（提取成功率、长度分布）。典型场景： motif 分析、顺式元件挖掘、启动子克隆、基因表达调控研究。

## 快速开始 | Quick Start

```bash
# 提取所有基因上游 2000 bp 启动子
biopytools promoter-extractor -g genes.gff3 --genome genome.fa -o promoters

# 提取 1500 bp，并只针对指定基因列表
biopytools promoter-extractor -g annotation.gff3 --genome asm.fa \
    -o my_promoters -p 1500 --gene-list target_genes.txt -t 8
```

## 参数说明 | Parameters

### 必需参数 | Required

| 参数 | 描述 |
|------|------|
| `-g, --gff` | 输入 GFF3 文件路径（需含 `gene` 特征）|
| `--genome` | 输入基因组 FASTA 文件路径 |

### 常用可选参数 | Common Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-o, --output` | `promoters` | 输出前缀 |
| `-p, --promoter-length` | `2000` | 启动子长度（bp）|
| `--min-length` | `0` | 最小可接受长度（bp），低于此值则跳过 |
| `--gene-list` | 无 | 基因 ID 列表文件（每行一个基因 ID）|
| `-t, --threads` | `1` | 线程数 |
| `--no-bed` | 关 | 不输出 BED 文件 |
| `--no-stats` | 关 | 不输出统计文件 |
| `-v, --verbose` | 关 | 详细输出（`-v`: INFO，`-vv`: DEBUG）|
| `--quiet` | 关 | 静默模式（仅 ERROR）|
| `-f, --force` | 关 | 强制覆盖已存在文件 |
| `--dry-run` | 关 | 模拟运行，不实际执行 |

（运行 `biopytools promoter-extractor -h` 查看完整参数列表）

## 输出 | Output

```
./
├── promoters.fa           # 启动子序列 FASTA（按基因转录方向）
├── promoters.bed          # 启动子在基因组上的坐标（BED6）
└── promoters_stats.txt    # 统计：提取基因数、成功率、长度分布等
```

FASTA 头为基因 ID；BED 第 6 列标识链方向，方便在 IGV / JBrowse 中可视化。

## 依赖 | Dependencies

- Python 3.7+（纯 Python 实现，使用标准库读取 FASTA/GFF）

## 引用 | Citation

- 启动子提取为通用生物信息学操作，无特定引用。GFF3 规范参考 Sequence Ontology Project。

## 相关链接 | References

- [项目主页](https://github.com/lixiang117423/biopytools)
