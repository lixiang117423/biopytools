# WGDI 比较基因组学 | WGDI Comparative Genomics

**基于 BLAST 的同源点图、共线性鉴定与 Ka/Ks 计算 | BLAST-based dotplot, collinearity detection, and Ka/Ks calculation.**

## 功能概述 | Overview

`wgdi` 模块封装了 WGDI (Whole-Genome Duplication Investigation) 工具，提供三个子命令，覆盖比较基因组学标准工作流：

| 子命令 | 说明 |
|------|------|
| `dotplot`     | 基于 BLAST + GFF + LENS 绘制两物种同源基因点图 |
| `collinearity` | 鉴定共线性区段（支持基因组级或染色体级模式） |
| `calks`       | 利用共线性结果对同源基因对计算 Ka/Ks |

下游 `calks` 依赖上游 `collinearity` 的输出文件，三者构成完整流水线。

## 快速开始 | Quick Start

```bash
# 1. 同源点图
biopytools wgdi dotplot \
    -b blast.txt \
    --gff1 species1.gff --gff2 species2.gff \
    --lens1 species1.lens --lens2 species2.lens \
    -o wgdi_out/

# 2. 共线性鉴定（依赖相同输入）
biopytools wgdi collinearity \
    -b blast.txt \
    --gff1 species1.gff --gff2 species2.gff \
    --lens1 species1.lens --lens2 species2.lens \
    -o wgdi_out/

# 3. 基于共线性的 Ka/Ks
biopytools wgdi calks \
    -c wgdi_out/collinearity.txt \
    --fasta1 species1.cds.fa --fasta2 species2.cds.fa \
    -o wgdi_out/
```

输入文件格式 | Input formats:
- `*.gff`: 基因注释，列为 `chr/id/start/end/order/strand`
- `*.lens`: 染色体统计，列为 `chr/length/gene_count`
- `blast.txt`: 标准BLAST m8/m8-tabular 输出

## 参数说明 | Parameters

### 必需参数 | Required

| 子命令 | 参数 | 描述 |
|------|------|------|
| `dotplot`/`collinearity` | `-b, --blast` | BLAST 结果文件 |
| `dotplot`/`collinearity` | `--gff1`, `--gff2` | 物种 1/2 的 GFF 文件 |
| `dotplot`/`collinearity` | `--lens1`, `--lens2` | 物种 1/2 的 LENS 文件 |
| `calks` | `-c, --collinearity` | 共线性结果文件 |
| `calks` | `--fasta1`, `--fasta2` | 物种 1/2 的 CDS FASTA 文件 |

### 常用可选参数 | Common Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-o, --output-dir` | `./wgdi_output` | 输出目录 |
| `-t, --threads` | `8` | 线程数 |
| `--wgdi-path` | 自动 | WGDI 软件路径 |
| `--score` | `100` | BLAST 分数阈值 |
| `--evalue` | `1e-5` | BLAST e-value 阈值 |
| `--multiple` | `1` | 每个基因保留的最佳同源数 |
| `--repeat-number` | `10`/`20` | 重复基因最大数量 |
| `--position` | `order` | 基因位置类型：`order`/`start`/`end` |
| `--comparison` (collinearity) | `genomes` | 比较模式：`genomes`/`chromosomes` |
| `--grading` (collinearity) | `50,40,25` | 评分（红,蓝,灰） |
| `--mg` (collinearity) | `40,40` | 最大 gap 值 |
| `--pvalue` (collinearity) | `1.0` | P 值阈值 |
| `--savefig` (dotplot) | `dotplot.png` | 输出图像文件 |
| `--savefile` (collinearity/calks) | `collinearity.txt`/`ks_results.txt` | 输出文件名 |

（运行 `biopytools wgdi <subcommand> -h` 查看完整参数列表）

## 输出 | Output

```
wgdi_output/
├── dotplot.png                # dotplot 子命令
├── collinearity.txt           # collinearity 子命令
├── ks_results.txt             # calks 子命令
└── 99_logs/                   # 运行日志
```

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-b, --blast` | — |  | BLAST结果文件｜BLAST result file |
| `--gff1` | — |  | 物种1 GFF文件｜Species 1 GFF file |
| `--gff2` | — |  | 物种2 GFF文件｜Species 2 GFF file |
| `--lens1` | — |  | 物种1 LENS文件｜Species 1 LENS file |
| `--lens2` | — |  | 物种2 LENS文件｜Species 2 LENS file |
| `-c, --collinearity` | — |  | 共线性分析结果文件｜Collinearity analysis result file |
| `--fasta1` | — |  | 物种1 FASTA文件｜Species 1 FASTA file |
| `--fasta2` | — |  | 物种2 FASTA文件｜Species 2 FASTA file |
| `-o, --output-dir` | `./wgdi_output` |  | 输出目录｜Output directory |
| `-t, --threads` | `8` | int | 线程数｜Number of threads |
| `--wgdi-path` | — |  | WGDI软件路径｜WGDI software path |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-b, --blast` | 必填 |  | BLAST结果文件｜BLAST result file |
| `--gff1` | 必填 |  | 物种1 GFF文件｜Species 1 GFF file |
| `--gff2` | 必填 |  | 物种2 GFF文件｜Species 2 GFF file |
| `--lens1` | 必填 |  | 物种1 LENS文件｜Species 1 LENS file |
| `--lens2` | 必填 |  | 物种2 LENS文件｜Species 2 LENS file |
| `--genome1-name` | `Genome1` |  | 物种1名称｜Species 1 name (default: Genome1) |
| `--genome2-name` | `Genome2` |  | 物种2名称｜Species 2 name (default: Genome2) |
| `--blast-reverse` | — | store_true | BLAST结果前两列互换｜Swap first two columns of BLAST result |
| `--multiple` | `1` | int | 最佳同源基因数｜Number of best homologous genes (default: 1) |
| `--score` | `100` | int | BLAST分数阈值｜BLAST score threshold (default: 100) |
| `--evalue` | `1e-05` | float | BLAST e-value阈值｜BLAST e-value threshold (default: 1e-5) |
| `--repeat-number` | `10` | int | 重复基因最大数量｜Max number of repetitive genes (default: 10) |
| `--position` | `order` | order/start/end | 基因位置类型｜Gene position type (default: order) |
| `--ancestor-left` | — |  | 左侧物种祖先染色体区域｜Ancestral chromosome region for left species |
| `--ancestor-top` | — |  | 顶部物种祖先染色体区域｜Ancestral chromosome region for top species |
| `--markersize` | `0.5` | float | 点大小｜Marker size (default: 0.5) |
| `--figsize` | `10,10` |  | 图像尺寸｜Figure size (default: 10,10) |
| `--savefig` | `dotplot.png` |  | 输出图像文件｜Output image file (default: dotplot.png) |
| `-o, --output-dir` | `./wgdi_output` |  | 输出目录｜Output directory (default: ./wgdi_output) |
| `-t, --threads` | `8` | int | 线程数｜Number of threads (default: 8) |
| `--comparison` | `genomes` | genomes/chromosomes | 比较模式｜Comparison mode (default: genomes) |
| `--process` | `8` | int | 进程数｜Number of processes (default: 8) |
| `--grading` | `50,40,25` |  | 评分参数(红,蓝,灰)｜Scoring parameters for red,blue,gray (default: 50,40,25) |
| `--mg` | `40,40` |  | 最大gap值｜Max gap values (default: 40,40) |
| `--pvalue` | `1.0` | float | P值阈值｜P-value threshold (default: 1.0) |
| `--savefile` | `collinearity.txt` |  | 输出文件｜Output file (default: collinearity.txt) |
| `-c, --collinearity` | 必填 |  | 共线性分析结果文件｜Collinearity analysis result file |
| `--fasta1` | 必填 |  | 物种1 FASTA文件｜Species 1 FASTA file |
| `--fasta2` | 必填 |  | 物种2 FASTA文件｜Species 2 FASTA file |
| `--wgdi-path` | — |  | WGDI软件路径｜WGDI software path |
| `-v, --version` | — | version |  |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- `WGDI` (Python 包，可通过 `--wgdi-path` 指定)
- `KaKs_Calculator`（calks 内部使用）
- `matplotlib`, `numpy`, `pandas`（WGDI 依赖）

## 引用 | Citation

- Sun J. et al. WGDI: A user-friendly toolkit for evolutionary analyses of whole-genome duplications. Molecular Plant. 2024. (WGDI)
- Wang DP. et al. Gamma-MYN: a new algorithm for estimating Ka and Ks. Biology Direct. 2009.

## 相关链接 | References

- [WGDI GitHub](https://github.com/SunPengChuan/WGDI)
- [项目主页](https://github.com/lixiang117423/biopytools)
