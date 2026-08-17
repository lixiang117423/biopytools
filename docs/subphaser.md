# SubPhaser 异源多倍体亚基因组分离 | SubPhaser Subgenome Phasing

**基于重复 k-mer 的异源多倍体亚基因组自动分离和命名 | Auto-phase subgenomes of allopolyploids using repetitive k-mers**

## 功能概述 | Overview

SubPhaser 是一款专为异源多倍体设计的亚基因组分离工具。它通过识别各亚基因组特异的重复 k-mer 标签，对染色体进行无偏聚类，从而将同源染色体分配到对应的亚基因组并完成系统性命名（A、B、C...）。该方法不依赖亲本参考基因组，对人工合成多倍体和自然多倍体均适用。

本模块封装了 SubPhaser 的完整流程，支持自动模式（不提供配置）和配置模式（提供亚基因组配置文件），并集成父本验证、LTR 插入时间估计、Circos 可视化、同源区块分析等扩展功能。流程默认通过 conda 环境 `SubPhaser` 调用外部工具。适用于小麦、棉花、油菜、草莓等异源多倍体基因组的亚基因组注释。

## 快速开始 | Quick Start

```bash
# 基本用法：二倍体化分离（2 个亚基因组）
biopytools subphaser -i tetraploid_genome.fa --nsg 2 -o ./subphaser_output

# 验证模式：提供双亲基因组
biopytools subphaser -i allo_genome.fa --nsg 2 \
    --parental-genomes parentA.fa parentB.fa

# 已知亚基因组分配，跳过聚类
biopytools subphaser -i genome.fa --nsg 2 --sg-assigned assignment.txt

# 指定多组输入基因组并禁用耗时步骤
biopytools subphaser -i g1.fa -i g2.fa --nsg 4 \
    --disable-ltr --disable-circos -t 32
```

## 参数说明 | Parameters

### 必需参数 | Required

| 参数 | 描述 |
|------|------|
| `-i, --genomes` | 基因组 FASTA 文件（可多个，重复 `-i` 或空格分隔）|
| `--nsg` | 亚基因组数量（>=2）|

### 常用可选参数 | Common Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-o, --output-dir` | `./subphaser_output` | 输出目录 |
| `--prefix` | - | 输出前缀 |
| `-t, --threads` | `24` | 线程数 |
| `-c, --sg-cfgs` | - | 亚基因组配置文件（不提供则自动模式）|
| `--parental-genomes` | - | 父本基因组（验证模式，需 2 个）|
| `--min-chrom-size` | `1000000` | 最小染色体长度（bp），过滤小 contigs |
| `-k, --kmer-size` | `15` | K-mer 大小 |
| `-f, --min-fold` | `2.0` | 最小倍数差异 |
| `-q, --min-freq` | `200` | 最小 k-mer 频率 |
| `--max-pval` | `0.05` | 最大 P 值 |
| `--replicates` | `1000` | Bootstrap 重复次数 |
| `--test-method` | `ttest_ind` | 统计检验方法（ttest_ind/kruskal/wilcoxon/mannwhitneyu）|
| `--disable-ltr` | 关 | 禁用 LTR 分析（大基因组耗时）|
| `--disable-circos` | 关 | 禁用 Circos 图 |
| `--disable-blocks` | 关 | 禁用同源区块分析 |
| `--just-core` | 关 | 仅运行核心 phasing |
| `--ltr-detectors` | - | LTR 检测工具（ltr_finder/ltr_harvest）|
| `--mu` | `13e-9` | 替换率/年 |
| `--window-size` | `1000000` | Circos 窗口大小（bp）|
| `--aligner` | `minimap2` | 同源区块比对工具（minimap2/unimap）|
| `--sg-assigned` | - | 已知亚基因组分配文件（跳过聚类）|
| `--target` | - | 目标染色体文件 |
| `--labels` | - | 基因组标签 |
| `--figfmt` | `pdf` | 图片格式（pdf/png）|
| `--overwrite` | 关 | 覆盖已有结果 |
| `--cleanup` | 关 | 清理临时文件 |
| `--conda-env` | `SubPhaser` | conda 环境名称 |

（运行 `biopytools subphaser -h` 查看完整参数列表）

## 输出 | Output

```
subphaser_output/
├── 01_kmers/            # 亚基因组特异 k-mer
├── 02_clustering/       # 染色体聚类结果
├── 03_phasing/          # 亚基因组分配与命名
├── 04_renamed/          # 重命名后的 FASTA
├── 05_ltr/              # LTR 插入时间分析（可选）
├── 06_circos/           # Circos 可视化（可选）
├── 07_blocks/           # 同源区块分析（可选）
├── 08_validation/       # 父本验证（可选）
├── 99_logs/             # 运行日志
└── phasing_result.txt   # 亚基因组分配结果汇总
```

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --genomes` | 必填 |  | 基因组FASTA文件(可多个)｜Genome FASTA files |
| `--nsg` | 必填 | int | 亚基因组数量(>=2)｜Number of subgenomes (>=2) |
| `-c, --sg-cfgs` | — |  | 亚基因组配置文件(可选，不提供则自动模式)｜Subgenome config files (optional, auto mode if omitted) |
| `--parental-genomes` | — |  | 父本基因组(验证模式，需2个)｜Parental genomes for validation (requires 2) |
| `-o, --output-dir` | `./subphaser_output` |  | 输出目录｜Output directory |
| `--prefix` | — |  | 输出前缀｜Output prefix |
| `-t, --threads` | `24` | int | 线程数｜Number of threads |
| `--min-chrom-size` | `1000000` | int | 最小染色体长度(bp)，过滤小contigs｜Min chromosome size (bp), filter small contigs |
| `-k, --kmer-size` | `15` | int | K-mer大小｜K-mer size |
| `-f, --min-fold` | `2.0` | float | 最小倍数差异｜Minimum fold difference |
| `-q, --min-freq` | `200` | int | 最小k-mer频率｜Minimum k-mer frequency |
| `--max-pval` | `0.05` | float | 最大P值｜Maximum P-value |
| `--replicates` | `1000` | int | Bootstrap重复次数｜Bootstrap replicates |
| `--test-method` | `ttest_ind` | ttest_ind/kruskal/wilcoxon/mannwhitneyu | 统计检验方法｜Statistical test method |
| `--disable-ltr` | — |  | 禁用LTR分析(大基因组耗时较长)｜Disable LTR analysis |
| `--disable-circos` | — |  | 禁用Circos图｜Disable Circos plot |
| `--disable-blocks` | — |  | 禁用同源区块｜Disable homologous blocks |
| `--just-core` | — |  | 仅运行核心phasing｜Only run core phasing |
| `--ltr-detectors` | — | ltr_finder/ltr_harvest | LTR检测工具｜LTR detection tools |
| `--mu` | `1.3e-08` | float | 替换率/年｜Substitution rate per year |
| `--window-size` | `1000000` | int | Circos窗口大小｜Circos window size (bp) |
| `--aligner` | `minimap2` | minimap2/unimap | 比对工具｜Aligner for homologous blocks |
| `--sg-assigned` | — |  | 已知亚基因组分配文件(跳过聚类)｜Pre-assigned subgenome file (skip clustering) |
| `--target` | — |  | 目标染色体文件｜Target chromosomes file |
| `--labels` | — |  | 基因组标签｜Genome labels |
| `--no-label` | — |  | 不添加标签前缀｜No label prefix |
| `--custom-features` | — |  | 自定义特征FASTA｜Custom feature FASTA files |
| `--figfmt` | `pdf` | pdf/png | 图片格式｜Figure format |
| `--overwrite` | — |  | 覆盖已有结果｜Overwrite existing results |
| `--cleanup` | — |  | 清理临时文件｜Clean up temporary files |
| `--conda-env` | `SubPhaser` |  | conda环境名称｜Conda environment name |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- SubPhaser（默认通过 conda 环境 `SubPhaser` 调用）
- jellyfish（k-mer 计数）
- LTR_finder / LTR_harvest（LTR 分析，可选）
- minimap2 / unmap（同源区块比对）
- Circos（可视化，可选）
- R（统计分析与绘图）

## 引用 | Citation

- Guo, R. et al. SubPhaser: a robust algorithm for detecting and phasing subgenomes in polyploids. *Genome Biology* 23, 89 (2022).

## 相关链接 | References

- [项目主页](https://github.com/lixiang117423/biopytools)
- [SubPhaser 官方仓库](https://github.com/sc-zhang/subphaser)
