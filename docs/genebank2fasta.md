# GenBank 序列提取 | GenBank to FASTA Extractor

**从一批 GenBank (.gb/.gbk) 注释文件中批量提取 CDS 核酸序列与蛋白序列 | Batch-extract CDS nucleotide and protein sequences from GenBank files**

## 功能概述 | Overview

`genebank2fasta` 用于对一批 GenBank 注释文件（例如从 NCBI 下载的叶绿体、线粒体、细菌或病毒基因组）进行批量结构化提取，输出按样本/基因分组的 CDS（.fasta）和蛋白质（.fasta）序列，方便后续的多序列比对、系统发育分析、同源基因聚类等工作。

工具内部基于 Biopython 的 `SeqIO` 解析 GenBank，识别每个 CDS 的 `/gene`、`/product`、`/locus_tag`、`/protein_id` 等注释，按基因名归类整理。支持最小蛋白长度过滤、是否按样本/基因分层输出、是否保留未识别基因等选项。开启 `--phylo` 后还会构建样品 × 基因的存在/缺失矩阵，便于挑选单拷贝直系同源基因构建物种树。

## 快速开始 | Quick Start

```bash
# 基本用法：从目录中读取所有 .gb 文件并提取
biopytools genebank2fasta -i /path/to/genbank -o ./output

# 构建系统发育矩阵，按样品和基因分离输出
biopytools genebank2fasta -i ./genbank_dir -o ./results -t 24 --phylo
```

## 参数说明 | Parameters

### 必需参数 | Required

| 参数 | 描述 |
|------|------|
| `-i, --input` | 输入目录，包含若干 GenBank 文件 |

### 常用可选参数 | Common Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-o, --output-dir` | `./genbank_output` | 输出目录 |
| `-t, --threads` | `12` | 并行线程数 |
| `--min-length` | `10` | 最小蛋白长度（氨基酸）|
| `--phylo` | 关 | 创建样品 × 基因系统发育分析矩阵 |
| `--no-sample-sep` | 关 | 不按样本分离输出 |
| `--no-gene-sep` | 关 | 不按基因分离输出 |
| `--keep-unknown` | 关 | 保留无法识别基因名的记录 |

（运行 `biopytools genebank2fasta -h` 查看完整参数列表）

## 输出 | Output

```
genbank_output/
├── cds/                       # 按基因分离的 CDS 核酸 FASTA
│   ├── geneA.fasta
│   └── geneB.fasta
├── pep/                       # 按基因分离的蛋白 FASTA
│   ├── geneA.fasta
│   └── geneB.fasta
├── (按样品分离目录，可选)
└── report / phylo_matrix.txt  # 统计与系统发育矩阵（--phylo 时）
```

默认按样品与基因两层目录分离，便于后续挑取同源基因集合；若序列过短（小于 `--min-length`）会被过滤。

## 依赖 | Dependencies

- Python 3.7+
- Biopython（`Bio.SeqIO` 解析 GenBank）
- pandas（统计与矩阵）

## 引用 | Citation

- Cock, P. J. et al. Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics* 25, 1422-1423 (2009).

## 相关链接 | References

- [项目主页](https://github.com/lixiang117423/biopytools)
