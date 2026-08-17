# Samplot 结构变异可视化 | Samplot Structural Variant Visualization

**基于 Samplot 对 BAM/CRAM 中的基因组结构变异 (SV) 进行可视化，支持单区域绘图与 VCF 批量绘图 | Visualize genomic structural variants in BAM/CRAM using Samplot, supporting single-region and VCF batch plotting.**

## 功能概述 | Overview

提供两个子命令：

- `plot`：针对单个基因组区域绘制 SV 图，适合人工核查候选变异
- `vcf`：根据 VCF 批量绘制多个样本的 SV，支持多线程并行、下采样、样本筛选与 GFF3 注释叠加

两个子命令都会调用外部 samplot 可执行文件，并写入运行日志。

## 快速开始 | Quick Start

```bash
# 单区域绘图
biopytools samplot plot -b sample.bam -c chr1 -s 1000 -e 5000 -t DEL

# 批量绘制VCF中的SV
biopytools samplot vcf --vcf variants.vcf -b sample.bam -d output/ -O png -t 4
```

## 子命令 plot 参数 | plot Parameters

### 必需参数 | Required

| 参数 | 描述 |
|------|------|
| `-b, --bams` | BAM/CRAM 文件路径（可多次指定或空格分隔多个） |
| `-c, --chrom` | 染色体名称 |
| `-s, --start` | 起始位置（1-based） |
| `-e, --end` | 结束位置 |

### 常用可选参数 | Common Options (plot)

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-t, --sv-type` | 无 | SV 类型：`DEL`/`DUP`/`INV`/`BND` |
| `-o, --output-file` | 无 | 输出文件名 |
| `--output-dir` | `.` | 输出目录 |
| `-r, --reference` | 无 | 参考基因组（CRAM 必需） |
| `-d, --max-depth` | `1` | 最大正常 pair 数 |
| `-w, --window` | 自动 | 窗口大小 |
| `-z, --z` | `4` | 离均值标准差倍数 |
| `-H, --plot-height` | 自动 | 图高 |
| `-W, --plot-width` | `8` | 图宽 |
| `--dpi` | `300` | 输出 DPI |
| `--long-read` | `1000` | 长读长最小长度 |
| `--coverage-only` | off | 仅显示覆盖度 |
| `--same-yaxis-scales` | off | 多样本统一 Y 轴 |
| `-n, --titles` | 无 | 样本标题（可多次指定） |
| `--samplot-path` | `~/miniforge3/envs/samplot_v.1.3.0/bin/samplot` | samplot 可执行文件路径 |

## 子命令 vcf 参数 | vcf Parameters

### 必需参数 | Required

| 参数 | 描述 |
|------|------|
| `-b, --bams` | BAM/CRAM 文件路径（可多次指定） |
| `--vcf` | VCF 文件路径 |

### 常用可选参数 | Common Options (vcf)

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-d, --output-dir` | `samplot-out` | 输出目录 |
| `-O, --output-type` | `png` | 输出格式：`png`/`pdf`/`eps`/`jpg` |
| `-t, --threads` | `1` | 线程数 |
| `--downsample` | `1` | 下采样数 |
| `--min-bp` | `20` | 最小 SV 长度 (bp) |
| `--max-mb` | 无 | 最大 SV 长度 (MB) |
| `--sample-ids` | 无 | 指定样本 ID 列表 |
| `--plot-all` | off | 绘制所有样本 |
| `--min-call-rate` | 无 | 最小 call rate |
| `--max-hets` | 无 | 最大杂合数 |
| `--min-entries` | `6` | 每张图最少样本数 |
| `--max-entries` | `10` | 每张图最多样本数 |
| `--gff3` | 无 | GFF3 注释文件（叠加到图上） |
| `--samplot-path` | `~/miniforge3/envs/samplot_v.1.3.0/bin/samplot` | samplot 可执行文件路径 |

（运行 `biopytools samplot plot -h` / `biopytools samplot vcf -h` 查看完整参数列表）

## 输出 | Output

- `plot`：在 `--output-dir`（或 `--output-file` 指定文件名）下生成单张图片
- `vcf`：在 `--output-dir` 下按 VCF 记录批量生成图片，文件扩展名由 `--output-type` 决定
- 两个子命令都会在输出目录生成 `samplot_plot.*.log` / `samplot_vcf.*.log` 日志文件

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-b, --bams` | 必填 |  | BAM/CRAM文件路径｜BAM/CRAM file paths |
| `-c, --chrom` | 必填 |  | 染色体名称｜Chromosome name |
| `-s, --start` | 必填 | int | 起始位置｜Start position |
| `-e, --end` | 必填 | int | 结束位置｜End position |
| `-t, --sv-type` | — |  | SV类型(DEL/DUP/INV/BND)｜SV type |
| `-o, --output-file` | — |  | 输出文件名｜Output file name |
| `--output-dir` | `.` |  | 输出目录｜Output directory |
| `-r, --reference` | — |  | 参考基因组(CRAM必须)｜Reference genome (required for CRAM) |
| `-d, --max-depth` | `1` | int | 最大正常pair数｜Max normal pairs to plot |
| `-w, --window` | — | int | 窗口大小｜Window size |
| `-z, --z` | `4` | int | 标准差倍数｜Number of stdevs from mean |
| `-H, --plot-height` | — | int | 图高｜Plot height |
| `-W, --plot-width` | `8` | int | 图宽｜Plot width |
| `--dpi` | `300` | int | DPI｜Dots per inch |
| `--long-read` | `1000` | int | 长读长最小长度｜Min length for long-read |
| `--coverage-only` | `False` |  | 仅显示覆盖度｜Show only coverage |
| `--same-yaxis-scales` | `False` |  | 统一Y轴｜Use same Y-axis scales |
| `-n, --titles` | — |  | 样本标题｜Sample titles |
| `--samplot-path` | `~/miniforge3/envs/viz/bin/samplot` |  | samplot路径｜samplot binary path |
| `--vcf` | 必填 |  | VCF文件路径｜VCF file path |
| `-d, --output-dir` | `samplot-out` |  | 输出目录｜Output directory |
| `-O, --output-type` | `png` | png/pdf/eps/jpg | 输出格式｜Output format |
| `-t, --threads` | `1` | int | 线程数｜Number of threads |
| `--downsample` | `1` | int | 下采样数｜Downsample count |
| `--min-bp` | `20` | int | 最小SV长度(bp)｜Min SV length in bp |
| `--max-mb` | — | int | 最大SV长度(MB)｜Max SV length in MB |
| `--sample-ids` | — |  | 样本ID列表｜Sample ID list |
| `--plot-all` | `False` |  | 绘制所有样本｜Plot all samples |
| `--min-call-rate` | — | float | 最小call rate｜Min call rate |
| `--max-hets` | — | int | 最大杂合数｜Max heterozygotes |
| `--min-entries` | `6` | int | 最小样本数｜Min entries to plot |
| `--max-entries` | `10` | int | 最大样本数｜Max entries to plot |
| `--gff3` | — |  | GFF3注释文件｜GFF3 annotation file |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-b, --bams` | 必填 |  | [FILE] BAM/CRAM文件路径｜BAM/CRAM file paths |
| `-c, --chrom` | 必填 |  | [STR] 染色体名称｜Chromosome name |
| `-s, --start` | 必填 | int | [INT] 起始位置｜Start position |
| `-e, --end` | 必填 | int | [INT] 结束位置｜End position |
| `-t, --sv-type` | — |  | [STR] SV类型(DEL/DUP/INV/BND)｜SV type |
| `-o, --output-file` | — |  | [FILE] 输出文件名｜Output file name |
| `--output-dir` | `.` |  | [DIR] 输出目录｜Output directory |
| `-r, --reference` | — |  | [FILE] 参考基因组(CRAM必须)｜Reference genome (required for CRAM) |
| `-d, --max-depth` | `1` | int | [INT] 最大正常pair数｜Max normal pairs to plot |
| `-w, --window` | — | int | [INT] 窗口大小｜Window size |
| `-z, --z` | `4` | int | [INT] 标准差倍数｜Number of stdevs from mean |
| `-H, --plot-height` | — | int | [INT] 图高｜Plot height |
| `-W, --plot-width` | `8` | int | [INT] 图宽｜Plot width |
| `--dpi` | `300` | int | [INT] DPI｜Dots per inch |
| `--long-read` | `1000` | int | [INT] 长读长最小长度｜Min length for long-read |
| `--coverage-only` | — | store_true | [FLAG] 仅显示覆盖度｜Show only coverage |
| `--same-yaxis-scales` | — | store_true | [FLAG] 统一Y轴｜Use same Y-axis scales |
| `-n, --titles` | — |  | [STR] 样本标题列表｜Sample title list |
| `--samplot-path` | `~/miniforge3/envs/viz/bin/samplot` |  | [FILE] samplot可执行文件路径｜samplot binary path |
| `--vcf` | 必填 |  | [FILE] VCF文件路径｜VCF file path |
| `-d, --output-dir` | `samplot-out` |  | [DIR] 输出目录｜Output directory |
| `-O, --output-type` | `png` | png/pdf/eps/jpg | [STR] 输出格式｜Output format |
| `-t, --threads` | `1` | int | [INT] 线程数｜Number of threads |
| `--downsample` | `1` | int | [INT] 下采样数｜Downsample count |
| `--min-bp` | `20` | int | [INT] 最小SV长度(bp)｜Min SV length in bp |
| `--max-mb` | — | int | [INT] 最大SV长度(MB)｜Max SV length in MB |
| `--sample-ids` | — |  | [STR] 样本ID列表｜Sample ID list |
| `--plot-all` | — | store_true | [FLAG] 绘制所有样本｜Plot all samples |
| `--min-call-rate` | — | float | [FLOAT] 最小call rate｜Min call rate |
| `--max-hets` | — | int | [INT] 最大杂合数｜Max heterozygotes |
| `--min-entries` | `6` | int | [INT] 最小样本数｜Min entries to plot |
| `--max-entries` | `10` | int | [INT] 最大样本数｜Max entries to plot |
| `--gff3` | — |  | [FILE] GFF3注释文件｜GFF3 annotation file |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- samplot (默认期望版本 v1.3.0，路径可通过 `--samplot-path` 覆盖)
- BAM/CRAM 对应的 `.bai`/`.crai` 索引
- CRAM 输入需要参考基因组 (`-r`)

## 引用 | Citation

- Camp E.R. et al. Samplot: drawing attention to complex structural variation. *Bioinformatics Advances*, 2022.

## 相关链接 | References

- [项目主页](https://github.com/lixiang117423/biopytools)
