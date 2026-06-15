# 装配质量 QV 计算 | Assembly Quality QV Calculation

**基于 Merqury 通过 k-mer 光谱分析计算基因组装配的 QV 值 | Calculate assembly QV via Merqury k-mer spectrum analysis**

## 功能概述 | Overview

QV（Quality Value）是衡量基因组装配错误率的核心指标，Q40 表示 1/10000 的错误率。本模块包装了 Merqury，利用独立于组装过程的 Illumina 或 HiFi 原始 reads 构建 k-mer 数据库，与装配结果比对，从而无参考地估计装配的准确性和完整性。

流程自动完成 meryl k-mer 数据库构建、k-mer 比对、QV 评估与 completeness 评估，并输出可读统计报告。适用于评估单倍型组装或双单倍型（hap1/hap2）组装的质量。

## 快速开始 | Quick Start

```bash
# 基本用法：Illumina reads 评估组装
biopytools assembly-qv -i fastq_dir/ -g genome.fa

# 指定 HiFi 数据评估
biopytools assembly-qv -i hifi.fq.gz -g genome.fa --data-type hifi -o ./qv_results

# 指定 k-mer 大小和线程
biopytools assembly-qv -i reads/ -g asm.fa -k 21 -t 24
```

## 参数说明 | Parameters

### 必需参数 | Required

| 参数 | 描述 |
|------|------|
| `-i, --input` | FASTQ 文件或目录（用于构建 k-mer 库的原始 reads）|
| `-g, --genome` | 待评估的基因组 FASTA 文件 |

### 常用可选参数 | Common Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-o, --output-dir` | `./assembly_qv_output` | 输出目录 |
| `-k, --kmer-size` | 自动 | K-mer 大小（不指定则自动选择）|
| `-t, --threads` | `12` | 线程数 |
| `--conda-env` | `~/miniforge3/envs/merqury_v.1.3/bin/` | Merqury conda 环境路径 |
| `--data-type` | `auto` | 数据类型（auto/illumina/hifi）|

（运行 `biopytools assembly-qv -h` 查看完整参数列表）

## 输出 | Output

```
assembly_qv_output/
├── meryl_db/                 # meryl k-mer 数据库
├── *.qv                      # QV 估计结果
├── *.completeness_stats      # 完整性统计
├── spectra-cn.*              # k-mer 拷贝数光谱图
└── assembly_qv.log           # 运行日志
```

## 依赖 | Dependencies

- Merqury（默认通过 conda 环境调用，`--conda-env` 可指定路径）
- meryl（Merqury 依赖，构建 k-mer 数据库）
- Java（Merqury 运行环境）

## 引用 | Citation

- Rhie, A., Walenz, B. P., Koren, S. & Phillippy, A. M. Merqury: reference-free quality, completeness, and phasing assessment for genome assemblies. *Genome Biology* 21, 245 (2020).

## 相关链接 | References

- [项目主页](https://github.com/lixiang117423/biopytools)
- [Merqury 官方仓库](https://github.com/marbl/merqury)
