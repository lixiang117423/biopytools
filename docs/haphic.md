# HapHiC 基因组 Scaffolding | HapHiC Genome Scaffolding

一句话理解：**用 Hi-C 数据把基因组碎片自动聚类、排序、定向，串成染色体级别的组装**，并支持等位基因感知（能区分二倍体/多倍体的两套染色体），解决「组装还是碎片、且分不清同源两套」的问题。

## 功能概述 | Overview

- 封装 HapHiC：快速、不依赖参考基因组的等位基因感知 scaffolding 工具
- 支持 BAM 输入（按 read 名排序）或 FASTQ 输入（自动 BWA 比对 + samblaster 去重 + filter_bam 过滤）
- 完整流程：聚类 → 重分配 → 排序 → 构建，外加组装校正、可视化、Juicebox 文件生成
- 内建组装校正（默认 2 轮），校正后自动重跑一遍 HapHiC 提升组装质量
- 断点续传默认开启，`--force-rerun` 可强制全量重跑
- 输出 AGP、FASTA、接触图、.hic/.assembly（进 Juicebox）等多格式

## 快速开始 | Quick Start

```bash
biopytools haphic -i assembly.fa -b hic.bam -c 24
```

最小输入：一个组装 FASTA、一个按 read 名排序的 Hi-C BAM、以及染色体数（-c）。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| scaffolding | 把碎片按正确顺序和方向串成完整染色体 |
| 等位基因感知 | 能认出「两条同源染色体各自是谁」，而不是把它们糊在一起 |
| 聚类（cluster） | 按接触信号把碎片分进不同的「团」（大致对应染色体） |
| 重分配（reassign） | 把聚类时放错位置的碎片重新归位 |
| 排序 / 定向 | 确定团内碎片的先后顺序和正反方向 |
| MAPQ | 比对质量分，越高越可信；阈值以下的对读段被丢弃 |
| 编辑距离 | 比对时允许的错配/插入删除次数上限 |
| 酶切位点（RE） | Hi-C 建库用的限制酶识别序列（默认 GATC） |
| bin | 把染色体切成的等长格子，用来统计接触 |
| 膨胀值（inflation） | MCL 聚类的「松紧」参数，越大团分得越细 |
| AGP | 记录每条染色体由哪些碎片按什么顺序方向拼成的清单 |
| 组装校正 | 用覆盖度找出错误连接并修正的过程 |

## 输入 | Input

### 组装 FASTA

`-i` / `--input`，FASTA 格式（.fa / .fasta）。

### Hi-C 数据（二选一）

- **BAM 模式**：`-b` / `--bam`，按 read 名排序（非坐标排序）的 Hi-C BAM。程序会自动检查是否已有索引。
- **FASTQ 模式**：`-1` / `--hic1` + `-2` / `--hic2`，程序自动做 BWA 比对、samblaster 去重、filter_bam 过滤。

两种模式不能混用，FASTQ 模式必须同时给 hic1 和 hic2。

```text
# BAM 输入
biopytools haphic -i assembly.fa -b hic_sorted_by_name.bam -c 24

# FASTQ 输入（自动比对）
biopytools haphic -i assembly.fa -1 hic_R1.fq.gz -2 hic_R2.fq.gz -c 24
```

## 参数说明 | Parameters

### 主要输入 | Main input

**通俗理解|In plain words:** 组装文件（-i）是必填；Hi-C 数据要么给 BAM（-b，需按 read 名排序），要么给 FASTQ 一对（-1/-2，自动比对）。染色体数（-c）按物种实际数填，它决定聚类目标组数。

### 输出配置 | Output

**通俗理解|In plain words:** `--output-dir` 默认是组装文件所在目录，`--prefix` 默认取组装文件名。想换输出位置或统一文件前缀时再指定。`--force-rerun` 会清空重跑（关闭断点续传），一般只在换参数或结果异常时用。

### Hi-C 数据处理 | Hi-C processing

**通俗理解|In plain words:** 这组参数控制「读段怎么过滤才可信」。MAPQ 阈值、编辑距离、RE 位点过滤阈值（`--re-site-cutoff`、`--min-RE-sites`）都是上游 HapHiC 的默认值，一般不用动。`--aln-format` 自动识别比对格式，BAM 输入时保持 auto 即可。

### 聚类参数 | Clustering

**通俗理解|In plain words:** 控制 MCL 聚类的松紧。膨胀值范围（`--min-inflation` 到 `--max-inflation`）是扫描区间，步长 `--inflation-step`；膨胀值越大分得越细（团越多）。`--Nx` 是过滤低覆盖的阈值，`--min-group-len` 是最小团长度（Mbp），`--flank` 是邻接矩阵侧翼区（kbp），`--bin-size-kbp` 是聚类分箱（-1=自动）。这些基本都用默认值。

### 性能参数 | Performance

**通俗理解|In plain words:** `--processes` 是并行进程数、`-t/--threads` 是线程数、`--memory-limit` 是内存上限。按机器资源调大 processes/threads 可提速，内存上限建议略低于机器实际可用内存。

### 组装修正 | Assembly correction

**通俗理解|In plain words:** 用覆盖度找出「连错的地方」并修正，`--correct-nrounds` 是修正轮数（默认 2，0=关闭），`--correct-min-coverage` 是最小覆盖度。其余参数（`--median-cov-ratio` 等）是覆盖度异常判定的内部阈值，一般不用动。修正后程序会自动重跑一遍 HapHiC。

### ALLHiC 优化 | ALLHiC optimization

**通俗理解|In plain words:** HapHiC 内部用 ALLHiC 遗传算法做排序优化。`--skip-allhic` 跳过整个优化、`--skip-ga` 只跳过遗传算法（更快但可能差一点）、`--skip-fast-sort` 跳过快速排序。其余开关（`--sort-by-input`、`--no-additional-rescue` 等）都是细调，默认即可。

### 单倍型分相 | Haplotype phasing

**通俗理解|In plain words:** 有相位信息（GFA）时用 `--gfa-files` 传入，`--phasing-weight` 是相位权重（0~1）。`--remove-allelic-links` 移除等位基因之间的干扰链接。没有 GFA 就不需要动这组。

### 可视化与 Juicebox | Visualization and Juicebox

**通俗理解|In plain words:** `--generate-plots` 生成接触图（bin 大小 `--bin-size`、最小 scaffold 长度 `--min-len`）。默认会生成 Juicebox 兼容文件（.hic/.assembly），`--no-juicebox` 关闭。`--bin-size` 默认 500（bp，用于接触图），和聚类 bin 不是一回事。

### 工具路径 | Tool paths

**通俗理解|In plain words:** 各软件默认取 conda 环境里的路径（haphic/bwa/samtools/samblaster/filter_bam/matlock/3d-dna），一般不用改；只有安装位置不同时才需要指定。

## 分析流程 | Pipeline

```text
输入组装FASTA + Hi-C数据
  -> [FASTQ模式] BWA比对 -> samblaster去重 -> filter_bam过滤
  -> Step 1: 聚类（HapHiC pipeline cluster，01_cluster/）
  -> Step 2: 重分配（02_reassign/）
  -> Step 3: 排序（03_sort/）
  -> Step 4: 构建（04_build/，产出 AGP + FASTA）
  -> [可选] 组装校正（corrected_asm.fa -> 对校正后组装重跑HapHiC，07_corrected_contig_haphic/）
  -> [可选] 生成接触图（05_plots/）
  -> [可选] 生成Juicebox文件（.hic/.assembly，06_juicebox/）
```

## 输出 | Output

```text
<output_dir>/
├── 00_pipeline_info/
│   └── software_versions.yml        # 软件版本信息
├── 00_mapping/                      # FASTQ模式：BWA比对产物
│   ├── HiC.bam
│   └── HiC.filtered.bam
├── 01_cluster/
│   ├── corrected_asm.fa             # 组装修正后的组装
│   └── hic.filtered.bam
├── 02_reassign/
│   └── final_groups/final_clusters.txt  # 最终聚类结果
├── 03_sort/                         # 各组的 .tour 排序方案
├── 04_build/
│   ├── <prefix>.fa                  # 最终染色体级组装（核心结果）
│   ├── <prefix>.agp                 # 挂载清单
│   └── <prefix>.raw.agp             # 原始挂载清单
├── 05_plots/                        # 接触图（PDF/PNG）
├── 06_juicebox/
│   ├── <prefix>.hic                 # Juicebox 接触文件
│   └── <prefix>.assembly            # Juicebox assembly 文件
├── 07_corrected_contig_haphic/      # 校正后重跑的完整流程（若启用校正）
└── 99_logs/
    ├── <prefix>_haphic.log          # 运行日志
    └── <prefix>_haphic_report.txt   # 运行报告
```

关键文件说明：

- **<prefix>.fa**：最终染色体级组装，最核心的结果。
- **<prefix>.agp**：挂载清单，记录每条染色体由哪些碎片按什么顺序方向拼成。
- **<prefix>.hic / <prefix>.assembly**：载入 Juicebox 做可视化浏览与手工校正。
- **corrected_asm.fa**：组装校正后的组装，比原始组装错误连接更少。

## 结果解读 | Interpreting Results

**通俗理解|In plain words:** 看 scaffolding 好坏，看三点：染色体数对不对、每条染色体内部是否连续（N 少不少）、接触图是否呈干净的方块。

- **<prefix>.fa**：序列条数应接近 -c 指定的染色体数；序列数远多于染色体数说明有碎片没挂上去。
- **<prefix>.agp**：每条染色体下碎片越多越碎，说明连续性差；碎片又长又少说明挂载干净。
- **05_plots 接触图**：理想是沿对角线出现若干整齐方块（每条染色体一个），方块外干净；噪点多说明有错配。
- **corrected_asm.fa**：若校正轮数 >0，校正后重跑的结果在 07_corrected_contig_haphic/，可对比校正前后 N 数量和挂载连续性。
- **运行报告**：99_logs/<prefix>_haphic_report.txt 汇总了配置、输出文件与大小。

## 参数选择建议 | Parameter Guidance

- **-c 染色体数**：按物种实际数填，务必准确，直接决定聚类组数。
- **输入类型**：有现成的按 read 名排序 BAM 用 -b（省去比对）；只有 FASTQ 用 -1/-2（自动比对）。
- **-t / --processes**：机器核多可调大加速；内存上限按机器实际给。
- **--correct-nrounds**：追求高连续性保持默认 2；想快速出结果可设 0 关闭校正。
- **--force-rerun**：换参数或数据重跑前用，否则断点续传会跳过已完成步骤。

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input, -i` | 必填 | Path | 基因组组装文件(FASTA)｜Genome assembly file path (FASTA format) |
| `--bam, -b` | — | Path | Hi-C BAM文件(按read名排序)｜Hi-C BAM file path (sorted by read name) |
| `--hic1, -1` | — | Path | Hi-C Read1文件(FASTQ)｜Hi-C Read1 file path (FASTQ format) |
| `--hic2, -2` | — | Path | Hi-C Read2文件(FASTQ)｜Hi-C Read2 file path (FASTQ format) |
| `--chr-number, -c` | `12` | int | 染色体数量｜Number of chromosomes |
| `--output-dir, -o` | `.` | Path | 输出目录路径｜Output directory path |
| `--prefix` | — |  | 输出文件前缀｜Output file prefix |
| `--force-rerun` | — |  | 强制重新运行所有步骤｜Force rerun all steps (disable resume mode) |
| `--mapq-threshold` | `1` | int | MAPQ阈值｜MAPQ threshold |
| `--edit-distance` | `3` | int | 编辑距离阈值｜Edit distance threshold |
| `--re-site-cutoff` | `5` | int | Step1 RE位点过滤阈值｜Step1 RE site filtering threshold |
| `--min-RE-sites` | `25` | int | Step2重分配最小RE位点数｜Step2 reassignment min RE sites |
| `--aln-format` | `auto` | auto/bam/pairs | 比对文件格式｜Alignment file format |
| `--min-inflation` | `1.1` | float | 最小膨胀值｜Min inflation |
| `--max-inflation` | `3.0` | float | 最大膨胀值｜Max inflation |
| `--inflation-step` | `0.1` | float | 膨胀值步长｜Inflation step |
| `--Nx` | `80` | int | Nx参数｜Nx parameter |
| `--min-group-len` | `5.0` | float | 最小组长度(Mbp)｜Min group length (Mbp) |
| `--flank` | `500` | int | 邻接矩阵侧翼区域(kbp)｜Adjacency matrix flank region (kbp) |
| `--bin-size-kbp` | `-1` | int | 聚类分箱大小(kbp),-1=自动｜Clustering bin size (kbp), -1=auto |
| `--processes` | `8` | int | 并行进程数｜Number of parallel processes |
| `--threads, -t` | `12` | int | 线程数｜Number of threads |
| `--memory-limit` | `100G` |  | 内存限制｜Memory limit (e.g., 64G, 300G) |
| `--correct-nrounds` | `2` | int | 组装修正轮数(0=禁用)｜Assembly correction rounds (0=disabled) |
| `--correct-min-coverage` | `10.0` | float | 修正最小覆盖度｜Correction min coverage |
| `--median-cov-ratio` | `0.2` | float | 覆盖率截断乘数｜Coverage cutoff multiplier |
| `--region-len-ratio` | `0.1` | float | 高覆盖区域长度比｜High-coverage region length ratio |
| `--min-region-cutoff` | `5000` | int | 高覆盖区域最小长度(bp)｜Min high-coverage region length (bp) |
| `--skip-fast-sort` | — |  | 跳过快速排序｜Skip fast sorting |
| `--skip-allhic` | — |  | 跳过ALLHiC优化｜Skip ALLHiC optimization |
| `--skip-ga` | — |  | 跳过ALLHiC遗传算法｜Skip ALLHiC genetic algorithm |
| `--sort-by-input` | — |  | 按输入顺序排序｜Sort output by input order |
| `--no-additional-rescue` | — |  | 跳过额外救援轮｜Skip additional rescue round |
| `--remove-concentrated-links` | — |  | 移除高密度集中链接｜Remove concentrated links |
| `--normalize-by-nlinks` | — |  | 按链接数归一化｜Normalize by number of links |
| `--dense-matrix` | — |  | 使用稠密矩阵｜Use dense matrix |
| `--remove-allelic-links` | — | int | 移除等位基因连锁数｜Remove allelic links count |
| `--phasing-weight` | `1.0` | float | 相位权重｜Phasing weight |
| `--gfa-files` | — |  | GFA文件路径(逗号分隔)｜GFA files path (comma-separated) |
| `--generate-plots` | — |  | 生成可视化图表｜Generate visualization plots |
| `--bin-size` | `500` | int | 接触图bin大小｜Contact map bin size |
| `--min-len` | `1.0` | float | 最小scaffold长度｜Min scaffold length |
| `--separate-plots` | — |  | 生成单独的图表｜Generate separate plots |
| `--haphic-bin` | `~/miniforge3/envs/hic/bin/haphic` |  | HapHiC可执行文件路径｜HapHiC executable path |
| `--bwa-bin` | `~/miniforge3/envs/align/bin/bwa` |  | BWA可执行文件路径｜BWA executable path |
| `--samtools-bin` | `~/miniforge3/envs/align/bin/samtools` |  | Samtools可执行文件路径｜Samtools executable path |
| `--samblaster-bin` | `~/miniforge3/envs/hic/bin/samblaster` |  | Samblaster可执行文件路径｜Samblaster executable path |
| `--haphic-filter-bam-bin` | `~/miniforge3/envs/hic/bin/filter_bam` |  | HapHiC filter_bam工具路径｜HapHiC filter_bam tool path |
| `--use-samblaster/--no-use-samblaster` | `True` |  | 使用samblaster去重｜Use samblaster deduplication |
| `--use-haphic-filter/--no-use-haphic-filter` | `True` |  | 使用HapHiC过滤｜Use HapHiC filtering |
| `--generate-juicebox/--no-generate-juicebox` | `True` |  | 生成Juicebox兼容文件｜Generate Juicebox compatible files |
| `--matlock-bin` | `~/miniforge3/envs/hic/bin/matlock` |  | Matlock可执行文件路径｜Matlock executable path |
| `--three-d-dna-dir` | `~/software/3d-dna` |  | 3D-DNA目录路径｜3D-DNA directory path |
| `--agp2assembly-script` | `~/software/3d-dna/utils/agp2assembly.py` |  | agp2assembly脚本路径｜agp2assembly script path |
| `--asm-visualizer-script` | `~/software/3d-dna/visualize/run-assembly-visualizer.sh` |  | asm-visualizer脚本路径｜asm-visualizer script path |
| `--RE` | `GATC` |  | 限制性酶切位点｜Restriction enzyme sites |
| `--quick-view` | — |  | 快速查看模式｜Quick view mode |
| `--no-agp` | — |  | 不输出AGP文件｜Don't output AGP file |
| `--no-fasta` | — |  | 不输出FASTA文件｜Don't output FASTA file |
| `--no-generate-plots` | — |  | 不生成可视化图表｜Don't generate visualization plots |
| `--no-juicebox` | — |  | 不生成Juicebox脚本｜Don't generate Juicebox script |
| `--verbose, -v` | — |  | 详细输出模式｜Verbose output mode |
| `--log-file` | — | Path | 日志文件路径｜Log file path |
| `--dry-run` | — |  | 测试模式,不执行实际命令｜Test mode, do not execute actual commands |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- haphic（默认 `~/miniforge3/envs/hic/bin/haphic`，HAPHIC_PATH 覆盖）
- bwa（默认 `~/miniforge3/envs/align/bin/bwa`，仅 FASTQ 模式，BWA_PATH 覆盖）
- samtools（默认 `~/miniforge3/envs/align/bin/samtools`，SAMTOOLS_PATH 覆盖）
- samblaster（默认 `~/miniforge3/envs/hic/bin/samblaster`）
- filter_bam（默认 `~/miniforge3/envs/hic/bin/filter_bam`，HapHiC 自带工具）
- matlock（默认 `~/miniforge3/envs/hic/bin/matlock`，生成 .hic 用）
- 3D-DNA 脚本（默认 `~/software/3d-dna`，含 agp2assembly.py 与 run-assembly-visualizer.sh）
- 均通过 conda run 自动检测并包装对应环境

## 常见问题 | FAQ

**Q1：断点续传怎么生效？换数据会不会误跳过？**
默认开启，按各步骤特征文件判断（如 01_cluster 的 corrected_asm.fa + paired_links.clm、04_build 的 <prefix>.fa/.agp）。换数据或参数想重跑，加 `--force-rerun` 清空重来。

**Q2：BAM 输入报格式错误？**
HapHiC 需要**按 read 名排序**的 BAM（不是坐标排序）。用 `samtools sort -n` 重新排序后再喂入。

**Q3：FASTQ 模式比 BAM 慢很多？**
正常。FASTQ 模式要先跑 BWA 比对 + 去重 + 过滤。已有 BAM 时直接用 -b 更省时间。

**Q4：.hic 文件没生成？**
生成 .hic 依赖 3D-DNA 的 run-assembly-visualizer.sh 与 Java。若脚本缺失或 Java 环境异常会警告并跳过，.assembly 和 AGP/FASTA 仍会生成。

**Q5：想关闭自动的组装校正？**
用 `--correct-nrounds 0` 禁用；这样不会产生 07_corrected_contig_haphic/ 目录，直接以原始组装出结果。

