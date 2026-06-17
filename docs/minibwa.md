# Minibwa 短读长比对 | Minibwa Short-read Alignment

**基于 Minibwa 的多模式短读长比对工具（标准 / Hi-C / BS-seq / 长读） | Multi-mode short-read aligner based on Minibwa (standard / Hi-C / BS-seq / long-read)**

## 功能概述 | Overview

封装 [Minibwa](https://github.com/huang-olive/minibwa) 比对工具，提供从索引构建、批量比对、比对后处理（排序、去重）、统计到覆盖度分析的完整流程。通过 `--mode` 一键切换四种典型场景，省去手工拼接参数：

- `standard`：常规短读长比对（Illumina WGS / RNA-seq reagent）
- `hic`：Hi-C 双端比对，保留配对信息用于后续挂载
- `meth`：BS-seq 重亚硫酸盐比对，适配甲基化分析
- `long`：长读长（PacBio / ONT）比对

## 快速开始 | Quick Start

```bash
# 标准短读长比对
biopytools minibwa -g genome.fa -i fastq_dir -o minibwa_out -t 16

# Hi-C 模式
biopytools minibwa -g genome.fa -i hic_fastq_dir -o hic_out --mode hic -t 24

# BS-seq 甲基化模式
biopytools minibwa -g genome.fa -i bs_fastq_dir -o meth_out --mode meth

# 长读长模式（PacBio HiFi / ONT）
biopytools minibwa -g genome.fa -i long_read_dir -o long_out --mode long --preset lr
```

## 参数说明 | Parameters

### 必需参数 | Required

| 参数 | 描述 |
|------|------|
| `-g, --genome` | 参考基因组 FASTA 文件 |
| `-i, --input` | FASTQ 输入目录（自动按 `-p` 模式匹配 R1/R2） |

### 常用可选参数 | Common Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-o, --output-dir` | `./minibwa_output` | 输出目录 |
| `-p, --pattern` | `_1.fq.gz` | R1 文件名匹配模式 |
| `--mode` | `standard` | 比对模式：`standard` / `hic` / `meth` / `long` |
| `-t, --threads` | `12` | 线程数 |
| `--preset` | `adap` | Minibwa `-x` 预设：`sr` / `lr` / `adap` |
| `-k, --min-seed` | `19` | 最小种子长度 |
| `-A, --match-score` | `2` | 匹配得分 |
| `-B, --mismatch-penalty` | `8` | 错配罚分 |
| `-O, --gap-open` | `12,23` | gap 开放罚分 |
| `-E, --gap-ext` | `2,1` | gap 延伸罚分 |
| `-R, --read-group` | - | SAM read group 行 |
| `--markdup` | False | 标记 PCR 重复 |
| `--remove-dup` | False | 移除重复（需配合 `--markdup`） |
| `--skip-coverage` | False | 跳过覆盖度分析 |
| `--min-base-quality` | `0` | 最小碱基质量 |
| `--min-mapping-quality` | `0` | 最小比对质量 |
| `--max-depth` | `0` | 最大深度限制（0=无限） |
| `--window-size` | `1000000` | 覆盖度窗口大小 (bp) |
| `--step-size` | `100000` | 覆盖度步长 (bp) |
| `--minibwa-path` | `~/software/minibwa/minibwa` | minibwa 二进制路径 |
| `--samtools-path` | `~/.local/bin/samtools` | samtools 二进制路径 |
| `--resume` | False | 断点续传，跳过已完成样本 |

（运行 `biopytools minibwa -h` 查看完整参数列表）

## 输出 | Output

```
minibwa_output/
├── index/                   # 参考基因组索引
├── bam/                     # 比对结果 BAM 文件（每个样本一个）
│   ├── sample1.sorted.bam
│   └── sample1.sorted.bam.bai
├── stats/                   # 比对统计
│   └── sample1.flagstat.txt
├── coverage/                # 覆盖度分析（除非 --skip-coverage）
│   ├── sample1.cov.tsv
│   └── sample1.depth.tsv
└── logs/                    # 运行日志
```

## 依赖 | Dependencies

- **minibwa**（[https://github.com/huang-olive/minibwa](https://github.com/huang-olive/minibwa)）
- **samtools**（用于排序、索引、flagstat、覆盖度计算）

可通过 `--minibwa-path` 和 `--samtools-path` 指定二进制路径。

## 引用 | Citation

- Huang W. Minibwa: a fast and lightweight reimplementation of BWA-MEM for short-read alignment.（详见上游仓库）
- Li H. Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. *arXiv:1303.3997*, 2013.

## 相关链接 | References

- [Minibwa 上游仓库](https://github.com/huang-olive/minibwa)
- [项目主页](https://github.com/lixiang117423/biopytools)
