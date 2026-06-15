# Wgsim基因组测序数据模拟 | Wgsim Genome Sequencing Simulation

**基于wgsim从参考基因组模拟 Illumina 双端测序 reads | Simulate Illumina paired-end reads from reference genomes using wgsim**

## 功能概述 | Overview

对输入的一个或多个参考基因组（`.fa` / `.fna` / `.fasta`）调用 wgsim 生成模拟双端 FASTQ 数据，自动完成质量行修正、gzip 压缩以及断点续跑（已存在的样本自动跳过）。

- 支持单文件或整个目录批量模拟，自动识别 `.fa`/`.fna`/`.fasta`
- 内置参数：reads数量、reads长度、测序错误率、突变率、插入片段内外距离
- 自动修正 wgsim 输出的低质量行：将所有质量字符替换为 `I`（Q40），避免下游 fastp 等工具因默认质量阈值过高而丢弃所有reads
- 自动 gzip 压缩输出为 `*_1.fq.gz` / `*_2.fq.gz`
- 完善的日志记录与断点续跑（同名输出已存在则跳过）

## 快速开始 | Quick Start

```bash
# 单文件模拟，使用默认参数（5000万对reads，150bp，错误率0.020）
biopytools wgsim -i genome.fna -o out_dir

# 目录批量模拟，并自定义reads数量与长度
biopytools wgsim -i genome_dir/ -o out_dir -N 10000000 -1 150

# 固定随机种子以便复现，并降低错误率
biopytools wgsim -i genome.fa -o out_dir -s 42 -e 0.01
```

## 参数说明 | Parameters

### 必需参数 | Required

| 参数 | 描述 |
|------|------|
| `-i, --input` | 输入基因组文件或目录（目录下自动扫描 `.fa` / `.fna` / `.fasta`） |
| `-o, --output-dir` | 输出目录（不存在时自动创建） |

### 常用可选参数 | Common Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-N, --num-reads` | `50000000` | 每个基因组模拟的 read 对数 |
| `-1, --read-length` | `150` | reads长度（bp），同时应用于R1和R2 |
| `-s, --seed` | `0` | 随机种子，相同种子结果可复现 |
| `-e, --error-rate` | `0.020` | 测序错误率（0-1，由wgsim按-e参数分配，输出后会被脚本修正为Q40） |
| `-r, --mutation-rate` | `0.001` | 参考基因组上的突变率（用于引入变异） |
| `-d, --outer-distance` | `500` | 双端reads的外部距离（fragment长度，bp） |
| `-D, --inner-distance` | `0` | 双端reads的内部距离（bp） |

（运行 `biopytools wgsim -h` 查看完整参数列表）

## 输出 | Output

对每个输入基因组 `<base>`（不含扩展名）生成：

- `<base>_1.fq.gz`：R1 端 reads（gzip 压缩）
- `<base>_2.fq.gz`：R2 端 reads（gzip 压缩）
- `wgsim.log`：运行日志（含模拟参数、每个样本成功/跳过/失败统计）

日志末尾会输出汇总：`成功=success, 跳过=skip, 失败=fail/total`。

## 依赖 | Dependencies

- 外部工具：`wgsim`（samtools 套件之一），默认查找路径 `~/miniforge3/envs/GATK_v.4.6.2.0/bin/wgsim`，可通过环境变量 `WGSIM_PATH` 覆盖
- 系统工具：`gzip`（用于压缩输出）
- 可通过 `conda install -c bioconda samtools` 或独立编译安装 wgsim

## 引用 | Citation

wgsim 由 Heng Li 编写，作为 samtools/htslib 套件的一部分发布。
- Li H. A statistical framework for SNP calling, mutation discovery, association mapping and population genetical parameter estimation from sequencing data. *Bioinformatics*. 2011, 27(21):2987-2993.
- 源码：https://github.com/lh3/wgsim

## 相关链接 | References

- [wgsim GitHub](https://github.com/lh3/wgsim)
- [samtools 项目](http://www.htslib.org/)
- [项目主页](https://github.com/lixiang117423/biopytools)
