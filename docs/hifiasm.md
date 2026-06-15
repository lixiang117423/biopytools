# HiFiasm 基因组组装 | HiFiasm Genome Assembly

**基于 PacBio HiFi reads 的基因组组装自动化流程 | Automated genome assembly pipeline based on PacBio HiFi reads**

## 功能概述 | Overview

HiFiasm 是一款主流的高质量基因组组装软件，擅长利用 PacBio HiFi 高准确度长读长数据进行单倍型分相组装。本模块对 hifiasm 进行了完整封装，支持纯 HiFi 组装、HiFi + ONT 整合组装、HiFi + Hi-C 染色体级别挂载等多种组合方案。

流程覆盖从输入数据校验、组装、purge_dup 冗余清理、GFA 转 FASTA，到 BUSCO/QUAST 质量评估、单倍型差异分析的完整链路，并自动整理日志、配置文件和结果报告。适用于二倍体、多倍体物种的参考基因组或泛基因组构建。

## 快速开始 | Quick Start

```bash
# 基本用法：纯 HiFi 组装
biopytools hifiasm -i sample.hifi.fq.gz -o ./hifiasm_results -p sample

# HiFi + Hi-C 染色体级别组装
biopytools hifiasm -i sample.hifi.fq.gz \
    --hi-c-1 hic_R1.fq.gz --hi-c-2 hic_R2.fq.gz \
    -o ./hifiasm_hic -p sample --hg-size 3g

# HiFi + ONT 整合组装并跳过评估
biopytools hifiasm -i sample.hifi.fq.gz --ont-reads ont.fq.gz \
    --skip-busco --skip-quast -t 32
```

## 参数说明 | Parameters

### 必需参数 | Required

| 参数 | 描述 |
|------|------|
| `-i, --input-reads` | HiFi 测序数据文件（FASTQ，可压缩）|

### 常用可选参数 | Common Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-o, --output-dir` | `./hifiasm_output` | 输出目录 |
| `-p, --prefix` | `sample` | 输出文件前缀 |
| `-t, --threads` | `12` | 线程数 |
| `--hg-size` | `auto` | 基因组大小估计（如 `1.4g`、`500m`）|
| `-l, --purge-level` | `3` | Purge 级别（0-3）|
| `--purge-max` | `65` | 最大 purge 覆盖度 |
| `-s, --similarity-threshold` | `0.75` | 单倍型相似性阈值 |
| `--ont-reads` | - | ONT 长读长数据文件（整合组装）|
| `--hi-c-1 / --hi-c-2` | - | Hi-C 双端数据文件（需同时提供）|
| `--extra-hifiasm-args` | `''` | 额外透传给 hifiasm 的参数 |
| `--skip-busco / --skip-quast` | 关 | 跳过对应质量评估 |
| `--busco-lineage` | `auto` | BUSCO 谱系数据集 |
| `--reference-genome` | - | 参考基因组（用于 QUAST 比对）|
| `--min-contig-length` | `1000` | 最小 contig 长度过滤 |
| `--assembly-type` | `auto` | 组装类型（auto/diploid/triploid/polyploid）|
| `--output-formats` | `both` | 输出格式（fasta/gfa/both，可多选）|
| `--memory` | `100` | 内存限制（GB）|
| `--max-runtime` | `48` | 最大运行时间（小时）|
| `--resume` | 关 | 恢复中断的分析 |
| `--dry-run` | 关 | 只打印命令不执行 |
| `--config-file` | - | YAML 配置文件路径 |

（运行 `biopytools hifiasm -h` 查看完整参数列表，包括各外部工具路径）

## 输出 | Output

```
hifiasm_output/
├── 01_assembly/        # hifiasm 原始 GFA 组装结果
├── 02_fasta/           # 转换后的单倍型 FASTA 文件
├── 03_purge/           # purge_dups 清理结果
├── 04_busco/           # BUSCO 完整性评估
├── 05_quast/           # QUAST 组装统计
├── 06_haplotype/       # 单倍型差异分析（可选）
├── logs/               # 运行日志
└── config.yaml         # 本次运行配置快照
```

## 依赖 | Dependencies

- hifiasm（默认在 PATH 中，可通过 `--hifiasm-path` 指定）
- samtools
- BUSCO（可选，`--busco-path`）
- QUAST（可选，`--quast-path`）
- purge_dups（用于冗余序列清理）

## 引用 | Citation

- Cheng, H. et al. HiFiasm: haplotype-resolved de novo assembly using PacBio HiFi reads. *Nature Methods* 18, 170-175 (2021).
- Cheng, H. et al. Haplotype-resolved de novo assembly using phased assembly graphs with hifiasm. *Nature Methods* 19, 632-635 (2022).

## 相关链接 | References

- [项目主页](https://github.com/lixiang117423/biopytools)
- [hifiasm 官方仓库](https://github.com/chhylp123/hifiasm)
