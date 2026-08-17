# Hi-C 全基因组热图 | Hi-C Whole Genome Heatmap

**使用 HiC-Pro 比对生成 contact 矩阵，PlotHiC 绘制全基因组热图 | HiC-Pro alignment + matrix generation with PlotHiC visualization**

## 功能概述 | Overview

本模块串联了 HiC-Pro 和 PlotHiC 两个工具，提供从原始 Hi-C 双端 reads 到全基因组接触矩阵热图的完整流程。HiC-Pro 负责酶切位点消化、双端比对、配对分类、binning 矩阵生成与 ICE 校正；PlotHiC 负责将校正后的矩阵渲染为 publication-ready 的全基因组热图。

流程支持 MboI、HindIII、NcoI、EcoRI、BamHI 等多种限制性内切酶，可同时输出多个分辨率的 contact map。适用于基因组组装挂载质量可视化、染色体三维结构分析、染色质 compartments 识别等场景。HiC-Pro 可直接运行或通过 singularity 镜像调用。

## 快速开始 | Quick Start

```bash
# 基本用法
biopytools hic-heatmap -i genome.fa -g EcA \
    -1 R1.fq.gz -2 R2.fq.gz -o ./hic_output

# 使用 HindIII 酶切，指定分辨率
biopytools hic-heatmap -i genome.fa -g hg19 \
    -1 R1.fq.gz -2 R2.fq.gz \
    --restriction-enzyme HindIII --resolution 50000 -t 32

# 通过 singularity 镜像运行 HiC-Pro
biopytools hic-heatmap -i asm.fa -g sample \
    -1 R1.fq.gz -2 R2.fq.gz \
    --hicpro-sif /path/to/hicpro.sif
```

## 参数说明 | Parameters

### 必需参数 | Required

| 参数 | 描述 |
|------|------|
| `-i, --genome` | 基因组 FASTA 文件 |
| `-g, --genome-id` | 基因组 ID（用于输出命名，如 hg19、mm10）|
| `-1, --fastq-r1` | R1 测序文件 |
| `-2, --fastq-r2` | R2 测序文件 |

### 常用可选参数 | Common Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-o, --output-dir` | `./hic_output` | 输出目录 |
| `-t, --threads` | `64` | 线程数 |
| `--max-memory` | `200` | HiC-Pro 最大内存（GB）|
| `--restriction-enzyme` | `MboI` | 限制性内切酶（MboI/HindIII/NcoI/EcoRI/BamHI）|
| `--bowtie2-idx` | 自动生成 | Bowtie2 索引路径 |
| `--bin-sizes` | `20000 40000 150000 500000 1000000` | Contact map bin 大小（空格分隔）|
| `--resolution` | `100000` | 热图分辨率（bp）|
| `--color-map` | `YlOrRd` | 颜色方案（PlotHiC）|
| `--dpi` | `300` | 图像 DPI |
| `--format` | `pdf` | 输出格式（pdf/png/svg 等）|
| `--bar-max` | `1` | 颜色条最大值（log 变换后）|
| `--hicpro-sif` | 空 | HiC-Pro singularity 镜像路径 |
| `--plothic-path` | `~/miniforge3/envs/plothic_v.1.0.0/bin/plothic` | PlotHiC 可执行文件路径 |
| `--force` | 关 | 强制重新运行 |
| `--verbose / --quiet` | 关 | 详细日志 / 仅错误 |

（运行 `biopytools hic-heatmap -h` 查看完整参数列表）

## 输出 | Output

```
hic_output/
├── bowtie_results/        # 双端比对 BAM
├── hic_results/
│   ├── matrix/            # 各分辨率的 contact 矩阵（含 ICE 校正）
│   └── mapC/              # PlotHiC 渲染的热图（pdf/png）
├── annotation/            # bed 格式 bin 注释
└── hic_heatmap.log        # 运行日志
```

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --genome` | 必填 |  | 基因组FASTA文件｜Genome FASTA file |
| `-g, --genome-id` | 必填 |  | 基因组ID（用于输出文件命名，如hg19, mm10）｜Genome ID (for output file naming, e.g., hg19, mm10) |
| `-o, --output-dir` | `./hic_output` |  | 输出目录｜Output directory |
| `-1, --fastq-r1` | 必填 |  | R1测序文件｜R1 sequencing file |
| `-2, --fastq-r2` | 必填 |  | R2测序文件｜R2 sequencing file |
| `-t, --threads` | `64` |  | 线程数｜Threads |
| `--max-memory` | `200` | int | HiC-Pro最大内存限制（GB）｜Maximum memory limit for HiC-Pro in GB |
| `--restriction-enzyme` | `MboI` |  | 限制性内切酶｜Restriction enzyme (default: MboI). Options: MboI, HindIII, NcoI, EcoRI, BamHI |
| `--bowtie2-idx` | — |  | Bowtie2索引路径（默认自动生成）｜Bowtie2 index path (auto-generated if not specified) |
| `--bin-sizes` | `20000 40000 150000 500000 1000000` |  | Contact map bin大小（空格分隔）｜Contact map bin sizes, space-separated |
| `--resolution` | `100000` | int | 热图分辨率（bp）｜Heatmap resolution in bp (default: 100000, 100kb) |
| `--color-map` | `YlOrRd` |  | 颜色方案｜Color scheme (PlotHiC default: YlOrRd) |
| `--dpi` | `300` | int | 图像分辨率｜Image DPI |
| `--format` | `pdf` |  | 输出格式｜Output format (pdf, png, svg, etc.) |
| `--bar-max` | `1` | int | 颜色条最大值｜Color bar maximum value (after log transform) |
| `--hicpro-sif` | `` |  | HiCPro singularity镜像路径（留空则直接使用HiC-Pro）｜HiCPro singularity image path (leave empty to use HiC-Pro directly) |
| `--plothic-path` | `~/miniforge3/envs/plothic_v.1.0.0/bin/plothic` |  | PlotHiC可执行文件路径｜PlotHiC executable path |
| `--force` | — |  | 强制重新运行｜Force rerun all steps |
| `--verbose` | — |  | 显示详细日志｜Verbose logging |
| `--quiet` | — |  | 仅显示错误｜Errors only |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --genome` | 必填 |  | 基因组FASTA文件｜Genome FASTA file |
| `-g, --genome-id` | 必填 |  | 基因组ID（用于输出文件命名，如hg19, mm10）｜Genome ID (for output file naming, e.g., hg19, mm10) |
| `-o, --output-dir` | `./hic_output` |  | 输出目录｜Output directory (default: ./hic_output) |
| `-1, --fastq-r1` | 必填 |  | R1测序文件｜R1 sequencing file |
| `-2, --fastq-r2` | 必填 |  | R2测序文件｜R2 sequencing file |
| `-t, --threads` | `64` | int | 线程数｜Threads (default: 64) |
| `--max-memory` | `200` | int | HiC-Pro最大内存限制（GB）｜Maximum memory limit for HiC-Pro in GB (default: 200) |
| `--restriction-enzyme` | `MboI` |  | 限制性内切酶名称｜Restriction enzyme name (default: MboI). Options: MboI, HindIII, NcoI, EcoRI, BamHI |
| `--bowtie2-idx` | — |  | Bowtie2索引路径（默认自动生成）｜Bowtie2 index path (auto-generated if not specified) |
| `--bin-sizes` | `20000 40000 150000 500000 1000000` |  | Contact map bin大小（空格分隔）｜Contact map bin sizes, space-separated (default: "20000 40000 150000 500000 1000000") |
| `--resolution` | `100000` | int | 热图分辨率（bp）｜Heatmap resolution in bp (default: 100000, 100kb) |
| `--color-map` | `YlOrRd` |  | 热图颜色方案｜Heatmap color scheme (default: YlOrRd, PlotHiC default) |
| `--dpi` | `300` | int | 图像分辨率｜Image resolution in DPI (default: 300) |
| `--format` | `pdf` |  | 输出格式｜Output format (default: pdf, options: pdf, png, svg, etc.) |
| `--bar-max` | `1` | int | 颜色条最大值（log变换后）｜Color bar maximum value after log transform (default: 1) |
| `--hicpro-sif` | `~/software/singularity/hicpro_latest.sif` |  | HiCPro singularity镜像路径｜HiCPro singularity image path |
| `--singularity-exec` | `~/miniforge3/envs/singularity_v.3.8.7/bin/singularity` |  | Singularity可执行文件路径｜Singularity executable path |
| `--plothic-path` | `~/miniforge3/envs/plothic_v.1.0.0/bin/plothic` |  | PlotHiC可执行文件路径｜PlotHiC executable path |
| `--force` | — | store_true | 强制重新运行所有步骤｜Force rerun all steps |
| `--verbose` | — | store_true | 显示详细日志｜Show verbose logs |
| `--quiet` | — | store_true | 仅显示错误日志｜Show error logs only |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- HiC-Pro（直接调用或通过 `--hicpro-sif` singularity 镜像）
- Bowtie2（HiC-Pro 内部比对）
- Samtools、Unix 工具集
- PlotHiC（默认通过 conda 环境调用，`--plothic-path` 可指定）
- Python + R（HiC-Pro 内部依赖）

## 引用 | Citation

- Servant, N. et al. HiC-Pro: an optimized and flexible pipeline for Hi-C data processing. *Genome Biology* 16, 259 (2015).
- Wolff, J. et al. Galaxy HiCExplorer 3: a web server for reproducible Hi-C, capture Hi-C and single-cell Hi-C data analysis, quality control and visualization. *Nucleic Acids Research* (2020).

## 相关链接 | References

- [项目主页](https://github.com/lixiang117423/biopytools)
- [HiC-Pro 官方仓库](https://github.com/nservant/HiC-Pro)
