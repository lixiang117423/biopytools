# ALLHiC 染色体挂载 | ALLHiC Genome Scaffolding

一句话理解：**用 Hi-C 测序数据，把基因组里一堆零散的 contig（片段）自动排好顺序、连成一条条完整的染色体**，解决「拼出了基因组但还是一条条碎片」的问题。

## 功能概述 | Overview

- 基于 ALLHiC v5.4（Asmkit 版），用 Hi-C 读段把 draft 组装挂载到染色体级别
- 8 步完整流程：数据准备 → 比对 → 等位基因检测 → 修剪 → 分区 → 提取矩阵 → 拯救 → 优化 → 构建 → 绘图 → JBAT 生成
- 内建等位基因检测与「拯救」（rescue），能把分区阶段没定位的 contig 捡回来
- 自动构建 BWA / samtools 索引，自动用 conda run 包装环境里的软件
- 断点续传：每一步按输出文件存在性判断，已完成的步骤自动跳过
- 产出可直接进 Juicebox Assembly Tools 手工校正的 JBAT 文件（.assembly / .hic / .links）

## 快速开始 | Quick Start

```bash
biopytools allhic -r genome.fa -1 hic_R1.fq.gz -2 hic_R2.fq.gz -k 12
```

最小输入：一个参考基因组 FASTA、一对 Hi-C 读段文件、以及染色体数目（-k）。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| contig | 组装出来的「碎片」，还没有被排到染色体上 |
| scaffolding / 挂载 | 把碎片按正确的顺序和方向串成完整染色体的过程 |
| Hi-C 读段 | 一种能反映「DNA 上哪些位置在空间上靠得近」的测序数据，近的片段之间读段多 |
| 比对（mapping） | 把 Hi-C 读段「贴」回参考基因组，找到每条读段来自哪里 |
| 酶切位点（RE） | Hi-C 建库时用的限制酶识别序列（默认 GATC），决定接触矩阵怎么切 bin |
| bin | 把染色体切成等长的「格子」，统计格子和格子之间的接触次数 |
| 接触矩阵 | 一张「格子 × 格子」的表，格子越近、相互作用越强 |
| 嵌合 contig | 一段 contig 里混了本不该连在一起的两段，需要修剪掉 |
| 分区（partition） | 把 contig 分到不同的染色体组 |
| 拯救（rescue） | 把上一步没被分进任何组的 contig 再捡回来归位 |
| AGP 文件 | 记录「每条染色体由哪些 contig 按什么顺序、什么方向拼成」的清单文件 |
| JBAT | Juicebox Assembly Tools 的输入格式，供人工在图形界面里微调组装 |

## 输入 | Input

### 参考基因组

FASTA 格式的 draft 组装（.fa / .fasta）。会在工作目录里被链接为 `draft.asm.fasta` 并构建 BWA 与 samtools 索引。

### Hi-C 读段

两个配对文件，分别用 `-1` / `-2` 指定，支持 .fastq.gz 压缩格式。内部会链接为 `reads_R1.fastq.gz` / `reads_R2.fastq.gz`。

```text
biopytools allhic -r draft.asm.fasta -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz -k 12
```

## 参数说明 | Parameters

### 必需参数 | Required

**通俗理解|In plain words:** 这 4 个参数告诉程序「组装在哪、Hi-C 数据在哪、最终想拼成多少条染色体」。染色体数（-k）按你物种的真实染色体数目填，它决定分区阶段把 contig 分成几组。

### 基本参数 | Basic

**通俗理解|In plain words:** 酶切位点默认 GATC（即 MboI/DpnII 建库），除非你建库用了别的酶，否则一般不用动；线程数按机器资源调大即可；工作目录是整套流程所有中间产物的存放地。

### 比对过滤 | Mapping filter

**通俗理解|In plain words:** `--mapq-step1` 是第一步比对时的最低比对质量门槛。值越大保留的读段越「干净」但越少。默认 1 对 Hi-C 很合适，一般不用动。

### 绘图参数 | Plotting

**通俗理解|In plain words:** `--bin-size` 是接触图里每个格子的长度（500k 即 50 万 bp），`--min-bin-size` 是绘热图时允许的最细格子。想看图更精细可调小，但耗时和内存会上升，一般用默认即可。

### 步骤开关 | Step toggles

**通俗理解|In plain words:** 一组 `--skip-*` 开关，用来跳过某一步。断点续传会自动跳过已完成步骤，这些开关主要用于「某一步出问题想绕过去」或「上一步已手工做完」的场景，正常从头跑不需要用。

## 分析流程 | Pipeline

```text
输入基因组FASTA + Hi-C reads (R1/R2)
  -> Step 0:   准备数据（bwa index + samtools faidx，链接输入）
  -> Step 1:   BWA比对（bwa mem -5SPM -> samtools sort -n 得到 sample.clean.bam）
  -> Step 1.5: 等位基因检测（allhic extract + minimap2 自比对 + allhic alleles）
  -> Step 2:   修剪嵌合contigs（ALLHiC_prune 得到 prunning.bam）
  -> Step 3:   分区（ALLHiC_partition 按 -k 分到各染色体组）
  -> Step 3.5: 提取接触矩阵（allhic extract 得到 sample.clean.clm）
  -> Step 4:   拯救未定位contigs（ALLHiC_rescue）
  -> Step 5:   优化contig顺序（allhic optimize，逐组出 groupN.tour）
  -> Step 6:   构建最终组装（ALLHiC_build 得到 groups.asm.fasta + groups.agp）
  -> Step 7:   绘制接触图谱（ALLHiC_plot）
  -> Step 8:   生成JBAT（asmkit 得到 groups.assembly / out.links / groups.hic）
```

## 输出 | Output

```text
allhic_output/
├── 00_data/                  # 输入数据的链接 + BWA/samtools 索引
├── 01_mapping/
│   └── sample.clean.bam      # 比对+排序后的 Hi-C BAM（后续所有步骤的输入）
├── 02_allele_table/
│   └── alleles.table         # 等位基因表
├── 03_pruning/
│   └── prunning.bam          # 修剪嵌合contig后的 BAM
├── 04_partition/
│   └── prunning.clusters.txt # contig 到染色体组的分区结果
├── 05_extract_matrix/
│   ├── sample.clean.clm      # 接触矩阵
│   └── sample.clean.counts_GATC.txt  # 各 contig 的酶切 counts
├── 06_rescue/
│   └── prunning.counts_GATC.<k>g1.txt  # 拯救后的分组计数
├── 07_optimize/
│   └── groupN.tour           # 每个染色体组的 contig 顺序方案
├── 08_build/
│   ├── groups.asm.fasta      # 最终染色体级组装（核心结果）
│   ├── groups.asm.fasta.fai  # FASTA 索引
│   └── groups.agp            # 挂载清单（contig 到染色体的对应关系）
├── 09_plot/                  # 接触热图
├── 10_jbat_asmkit/
│   ├── groups.assembly       # Juicebox 用的 assembly 文件
│   ├── out.links             # Hi-C 链接文件
│   └── groups.hic            # .hic 接触文件（进 Juicebox 浏览）
└── logs/                     # 运行日志（按时间戳命名）
```

关键文件说明：

- **groups.asm.fasta**：最终染色体级组装，是最重要的结果。
- **groups.agp**：挂载清单，记录每条染色体由哪些 contig、按什么顺序方向拼成，是审查和手工校正的依据。
- **groups.assembly / out.links / groups.hic**：JBAT 三件套，一起载入 Juicebox Assembly Tools 做可视化手工微调。

## 结果解读 | Interpreting Results

**通俗理解|In plain words:** 看挂载好坏，核心是「染色体数对不对、每条染色体内部是否连续、有没有 contig 大量丢失」。

- **groups.asm.fasta**：序列条数应接近 -k 指定的染色体数（可能略多，因为有些短片段没挂上去）。序列数远多于染色体数，说明有 contig 没挂上去。
- **groups.agp**：每条染色体下面排的 contig 越多、越碎，说明组装连续性越差；反之 contig 又长又少说明挂载干净。
- **09_plot 里的热图**：理想情况是沿对角线出现若干整齐的方块（每条染色体一个方块），方块外很「干净」；方块外噪点多说明有错配或嵌合。
- **07_optimize/groupN.tour**：每个组的排序方案，若某组 tour 文件为空或缺失，说明该组优化失败（日志里会报「组 N 优化失败」）。

## 参数选择建议 | Parameter Guidance

- **-k 染色体数**：按物种真实染色体数填，填错会导致分区错误，务必确认。
- **-e 酶切位点**：默认 GATC 对应 MboI/DpnII 建库；若建库用 HindIII（AAGCTT）等，务必改，否则接触矩阵算错。
- **-t 线程数**：BWA 比对与排序是耗时大头，机器核多可调大（如 24~48）。
- **--bin-size**：大基因组（>1 Gb）可保持 500k；小基因组可适当调小让图更细，但不必一开始就调。
- **--skip-***：首次运行保持全开（都不 skip）；某步失败排查后，可用 skip 跳过已成功的上游步骤继续。

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--reference, -r` | 必填 | Path | 参考基因组文件｜Reference genome file |
| `--read1, -1` | 必填 | Path | Hi-C读段1文件｜Hi-C read 1 file |
| `--read2, -2` | 必填 | Path | Hi-C读段2文件｜Hi-C read 2 file |
| `--chr-num, -k` | 必填 | int | 染色体数目｜Number of chromosomes |
| `--enzyme, -e` | `GATC` |  | 酶切位点｜Enzyme motif |
| `--threads, -t` | `12` | int | 线程数｜Number of threads |
| `--workdir, -w` | `./allhic_output` |  | 工作目录｜Working directory |
| `--mapq-step1` | `1` | int | 步骤1比对质量阈值｜Step 1 MapQ threshold |
| `--bin-size` | `500k` |  | Bin大小｜Bin size |
| `--min-bin-size` | `50k` |  | 最小Bin大小｜Minimum bin size |
| `--skip-mapping` | — |  | 跳过步骤1:比对｜Skip Step 1: Mapping |
| `--skip-allele` | — |  | 跳过步骤1.5:等位基因检测｜Skip Step 1.5: Allele Detection |
| `--skip-prune` | — |  | 跳过步骤2:修剪｜Skip Step 2: Pruning |
| `--skip-partition` | — |  | 跳过步骤3:分区｜Skip Step 3: Partition |
| `--skip-extract` | — |  | 跳过步骤3.5:提取矩阵｜Skip Step 3.5: Extract Matrix |
| `--skip-rescue` | — |  | 跳过步骤4:救援｜Skip Step 4: Rescue |
| `--skip-optimize` | — |  | 跳过步骤5:优化｜Skip Step 5: Optimization |
| `--skip-build` | — |  | 跳过步骤6:构建FASTA｜Skip Step 6: Build FASTA |
| `--skip-plot` | — |  | 跳过步骤7:绘制热图｜Skip Step 7: Plot Heatmap |
| `--skip-asmkit` | — |  | 跳过步骤8:JBAT生成｜Skip Step 8: JBAT Generation |
| `--diagnose` | — |  | 诊断模式｜Diagnostic mode |
| `--verbose, -v` | — |  | 详细输出｜Verbose output |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-r, --reference` | 必填 |  | 参考基因组文件｜Reference genome file |
| `-1, --read1` | 必填 |  | Hi-C读段1文件｜Hi-C read 1 file |
| `-2, --read2` | 必填 |  | Hi-C读段2文件｜Hi-C read 2 file |
| `-k, --chr-num` | 必填 | int | 染色体数量｜Number of chromosomes |
| `-e, --enzyme` | `GATC` |  | ALLHiC酶切位点｜ALLHiC enzyme motif |
| `-t, --threads` | `88` | int | CPU线程数｜CPU threads |
| `-w, --workdir` | `./allhic_output` |  | 工作目录｜Working directory |
| `--mapq-step1` | `1` | int | Step 1 MapQ阈值｜Step 1 MapQ threshold |
| `--bin-size` | `500k` |  | Bin大小｜Bin size |
| `--min-bin-size` | `50k` |  | 最小Bin大小｜Minimum bin size |
| `--skip-mapping` | — | store_true | 跳过步骤1: 比对｜Skip Step 1: Mapping |
| `--skip-allele` | — | store_true | 跳过步骤1.5: 等位基因检测｜Skip Step 1.5: Allele Detection |
| `--skip-prune` | — | store_true | 跳过步骤2: 修剪｜Skip Step 2: Pruning |
| `--skip-partition` | — | store_true | 跳过步骤3: 分区｜Skip Step 3: Partition |
| `--skip-extract` | — | store_true | 跳过步骤3.5: 矩阵提取｜Skip Step 3.5: Extract Matrix |
| `--skip-rescue` | — | store_true | 跳过步骤4: 救援｜Skip Step 4: Rescue |
| `--skip-optimize` | — | store_true | 跳过步骤5: 优化｜Skip Step 5: Optimization |
| `--skip-build` | — | store_true | 跳过步骤6: 构建｜Skip Step 6: Build FASTA |
| `--skip-plot` | — | store_true | 跳过步骤7: 绘图｜Skip Step 7: Plot Heatmap |
| `--skip-asmkit` | — | store_true | 跳过步骤8: JBAT生成｜Skip Step 8: JBAT Generation |
| `--diagnose` | — | store_true | 运行诊断模式｜Run diagnostic mode |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- BWA（默认 `~/miniforge3/envs/Population_genetics/bin`，可用环境变量 ALLHIC_BWA_PATH 覆盖）
- ALLHiC 工具集（默认 `~/software/ALLHiC/bin`，含 allhic / ALLHiC_prune / ALLHiC_partition / ALLHiC_rescue / ALLHiC_build / ALLHiC_plot；ALLHIC_SOFTWARE_PATH / ALLHIC_BIN_PATH 覆盖）
- asmkit（默认 `~/software/asmkit/asmkit`，ALLHIC_ASMKIT_PATH 覆盖）
- 3D-DNA 可视化脚本（默认 `~/software/3d-dna/visualize/run-assembly-visualizer.sh`）
- samtools、minimap2（需在 PATH 中）
- conda 环境：运行时自动激活 `biopytools` 环境

## 常见问题 | FAQ

**Q1：断点续传怎么生效？换参数重跑为什么没变化？**
每个步骤都按「关键输出文件是否存在」判断是否跳过（如 sample.clean.bam、alleles.table、groups.asm.fasta 等）。换参数重跑同一个工作目录前，要先删掉对应步骤的旧产物，否则会复用旧结果。

**Q2：提示「缺少必需软件」但软件明明装了？**
依赖检查只做存在性判断，找不到时会打印警告但继续执行（在运行时报真实错误）。确认对应工具路径正确，或通过 ALLHIC_* 环境变量指向实际位置。

**Q3：Step 5 优化报「N 个组优化失败」？**
通常是该组 contig 数量太少或接触矩阵信号不足导致 allhic optimize 无解。检查该组的 groupN.txt 与接触矩阵，确认上游分区/拯救结果正常。

**Q4：groups.hic 没生成？**
Step 8 里 .hic 依赖 3D-DNA 的 run-assembly-visualizer.sh，找不到该脚本时会跳过并提示手动运行；.assembly 和 out.links 仍会生成。

**Q5：运行目录里文件名不是我起的名字？**
ALLHiC 流程内部约定固定文件名（draft.asm.fasta、reads_R1.fastq.gz、sample.clean.bam 等），你的输入会被链接成这些名字，属正常现象。

