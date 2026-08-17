# 泛基因组Block构建 | Pan-genome Block Construction

**通过迭代两两比对构建泛基因组共线性块 | Build pan-genome syntenic blocks via iterative pairwise alignments.**

## 功能概述 | Overview

`pan_blocks` 用于从一组基因组出发，通过 `nucmer`/`minimap2` 完成两两全基因组比对，再基于 `bedtools` 等工具迭代地构建跨基因组的共线性 Block（泛基因组 Block），最终生成共线性可视化图与统计报告。流程被拆分为 `align`/`build`/`plot` 三个可单独运行的步骤，也可通过 `all` 一次性串行执行。同时提供 `prepare` 子命令，从 FASTA 目录自动生成基因组列表文件。

支持流程 | Pipeline:
- `prepare`: 扫描 FASTA 目录生成 `genome_list.txt`
- `align`: 两两全基因组比对 (MUMmer nucmer / minimap2)
- `build`: 基于比对构建 pan-block 矩阵
- `plot`: 绘制跨基因组共线性图 (SVG/PNG)
- `all`: 顺序执行 align + build + plot

## 快速开始 | Quick Start

```bash
# 1. 自动生成基因组列表
biopytools pan-blocks prepare -i /path/to/genomes/ -o genome_list.txt

# 2. 一键运行完整流程
biopytools pan-blocks all -i genome_list.txt -o output_dir/ -t 24

# 或者分步骤运行
biopytools pan-blocks align -i genome_list.txt -o output_dir/ -t 24
biopytools pan-blocks build  -i genome_list.txt -o output_dir/ --genome-order order.txt
biopytools pan-blocks plot   -i genome_list.txt -o output_dir/ --plot-format png
```

`genome_list.txt` 格式 | Genome list format（TAB 分隔）：

```
sampleA<TAB>/abs/path/sampleA.fa
sampleB<TAB>/abs/path/sampleB.fa
```

## 参数说明 | Parameters

### 必需参数 | Required

| 参数 | 描述 |
|------|------|
| `-i, --genome-list` | 基因组列表文件，TAB 分隔 (name<TAB>path) |
| `-i, --input-dir` (prepare) | FASTA 文件目录 |

### 常用可选参数 | Common Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-o, --output-dir` | `./pan_blocks_output` | 输出目录 |
| `-t, --threads` | `12` | 线程数 |
| `--parallel-alignments` | `4` | 同时并行执行的比对任务数 |
| `--min-alignment-length` | `10000` | delta-filter 最小比对长度 |
| `--genome-order` | - | 基因组优先级顺序文件 (每行一个名称) |
| `--chromosome` | - | 只处理指定染色体 |
| `--plot-format` | `svg` | 绘图格式：`svg` 或 `png` |
| `--nucmer` / `--delta-filter` / `--show-coords` / `--bedtools` / `--minimap2` | 自动 | 对应外部可执行文件路径，可手动覆盖 |

（运行 `biopytools pan-blocks <subcommand> -h` 查看完整参数列表）

## 子命令 | Subcommands

| 子命令 | 说明 |
|------|------|
| `prepare` | 从 FASTA 目录生成 `genome_list.txt` |
| `align`  | 仅运行两两比对 |
| `build`  | 仅构建 pan-blocks |
| `plot`   | 仅可视化 |
| `all`    | 完整流程 |

## 输出 | Output

输出目录下的典型结构（按子命令产出）：

```
output_dir/
├── 00_pipeline_info/
│   └── software_versions.yml      # 软件版本和运行参数记录
├── 01_alignments/                 # 两两比对 delta/coords 文件
├── 02_pan_blocks/                 # pan-block 矩阵及统计
├── 03_plots/                      # 共线性可视化 SVG/PNG
└── 99_logs/                       # 运行日志
```

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --genome-list` | 必填 |  | 基因组列表文件｜Genome list file (name<TAB>path) |
| `-o, --output-dir` | `./pan_blocks_output` |  | 输出目录｜Output directory |
| `-t, --threads` | `12` | int | 线程数｜Number of threads |
| `--parallel-alignments` | `4` | int | 并行比对数｜Parallel alignments |
| `--min-alignment-length` | `10000` | int | 最小比对长度｜Min alignment length |
| `--genome-order` | — |  | 基因组优先级顺序文件｜Genome priority order file |
| `--chromosome` | — |  | 指定染色体｜Specific chromosome |
| `--plot-format` | `svg` | svg/png | 绘图格式｜Plot format |
| `-i, --input-dir` | 必填 |  | FASTA文件目录｜FASTA files directory |
| `-o, --output` | `./genome_list.txt` |  | 输出文件路径｜Output file path |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --genome-list` | — |  | 基因组列表文件｜Genome list file (name<TAB>path) |
| `-o, --output-dir` | `./pan_blocks_output` |  | 输出目录｜Output directory |
| `-t, --threads` | `12` | int | 线程数｜Number of threads (default: 12) |
| `--parallel-alignments` | `4` | int | 并行比对数｜Number of parallel alignments (default: 4) |
| `--min-alignment-length` | `10000` | int | 最小比对长度｜Minimum alignment length for delta-filter (default: 10000) |
| `--genome-order` | — |  | 基因组优先级顺序文件｜Genome priority order file (one name per line) |
| `--chromosome` | — |  | 指定染色体｜Specific chromosome to process |
| `--step` | — | align/build/plot | 执行特定步骤｜Run specific step |
| `--nucmer` | — |  | nucmer可执行文件路径｜nucmer executable path |
| `--delta-filter` | — |  | delta-filter可执行文件路径｜delta-filter executable path |
| `--show-coords` | — |  | show-coords可执行文件路径｜show-coords executable path |
| `--bedtools` | — |  | bedtools可执行文件路径｜bedtools executable path |
| `--minimap2` | — |  | minimap2可执行文件路径｜minimap2 executable path |
| `--plot-format` | `svg` | svg/png | 绘图格式｜Plot format (default: svg) |
| `--plot-width` | `20` | int | 绘图宽度｜Plot width (default: 20) |
| `--plot-height` | `10` | int | 绘图高度｜Plot height (default: 10) |
| `--prepare` | — |  | 从目录自动生成genome_list.txt｜Auto-generate genome_list.txt from FASTA directory |
| `--prepare-output` | `./genome_list.txt` |  | prepare输出文件路径｜Output file path for prepare (default: ./genome_list.txt) |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- `MUMmer` (nucmer, delta-filter, show-coords)
- `bedtools`
- `minimap2`（可选）
- Python: `PyYAML`, `pandas`, `matplotlib` 等

## 引用 | Citation

- Kurtz S. et al. Versatile and open software for comparing large genomes. Genome Biol. 2004. (MUMmer)
- Quinlan AR & Hall IM. BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics. 2010.
- Li H. Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics. 2018.

## 相关链接 | References

- [项目主页](https://github.com/lixiang117423/biopytools)
