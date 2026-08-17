# Hi-C 全基因组热图 | Hi-C Heatmap (HiCPro + PlotHiC)

一句话理解：**把 Hi-C 测序的 FASTQ 先比对到基因组、生成"染色质互作矩阵"，再画成全基因组热图**，一眼看出染色体之间和染色体内部的相互作用强弱。

## 功能概述 | Overview

- 完整流程：HiCPro 比对 + 建矩阵 → PlotHiC 画全基因组热图
- 自动生成染色体大小文件（samtools faidx）和酶切片段文件（按限制酶识别序列扫描）
- 自动检测 / 构建 bowtie2 索引
- 断点续传：染色体大小、酶切片段、HiCPro 矩阵、热图均已存在且非空时跳过（`--force` 重跑）
- 支持 MboI / HindIII / NcoI / EcoRI / BamHI 五种限制酶；多种 bin 分辨率矩阵

## 快速开始 | Quick Start

```bash
biopytools hic-heatmap -i genome.fa -g EcA -1 R1.fq.gz -2 R2.fq.gz -o output
```

最小输入：基因组 FASTA + 基因组 ID（用于命名）+ 一对 Hi-C 测序 FASTQ。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| Hi-C | 一种测序技术，能测出"基因组里哪些区域在空间上靠得近" |
| 互作矩阵(contact matrix) | 每对基因组区间之间的接触次数表 |
| bin | 把基因组切成固定大小的格子（bp），矩阵的最小单元 |
| 分辨率(resolution) | bin 的大小；越小越精细，但矩阵越稀疏 |
| 限制性内切酶 | 在特定序列切 DNA 的酶，Hi-C 实验靠它产生连接位点 |
| 热图(heatmap) | 用颜色深浅表示互作强弱 |

## 输入 | Input

- `-i/--genome`：基因组 FASTA（自动建索引取长度）
- `-g/--genome-id`：基因组 ID，用于输出文件命名（如 hg19、EcA）
- `-1/--fastq-r1`、`-2/--fastq-r2`：Hi-C 双端测序文件（R1/R2）

> 限制酶必须与实验建库时实际使用的酶一致，默认 MboI。

## 参数说明 | Parameters

### 必需参数 | Required

**通俗理解|In plain words:** 四个参数缺一不可。`-g` 只是一个命名用的"代号"（决定输出文件名前缀，也决定 HiC-Pro 参考名），跟基因组实际内容无关，取个短名即可。

### HiCPro 参数 | HiCPro

**通俗理解|In plain words:** `-t` 线程数默认 64、`--max-memory` 内存上限默认 200GB，按你的机器配置调。`--restriction-enzyme` **必须和实验用的酶一致**（默认 MboI）。`--bin-sizes` 是建矩阵时的 bin 大小列表（空格分隔），`--bowtie2-idx` 可指定现成的 bowtie2 索引前缀，不填则自动在输出目录构建。**bin-sizes 一般不用动**，除非需要特定分辨率。

### PlotHiC 参数 | PlotHiC

**通俗理解|In plain words:** `--resolution` 是热图的分辨率（默认 100kb），程序会从 HiCPro 结果里挑最接近的 bin；`--color-map` 配色、`--dpi` 清晰度、`--format` 输出格式、`--bar-max` 颜色条上限（log 变换后）。**默认 100kb 全基因组够用**，想出更细的图再把 resolution 调小（但要求 HiCPro 已算了对应 bin）。

### 工具路径 | Tool paths

**通俗理解|In plain words:** `--hicpro-sif` 是 HiCPro 的 Singularity 镜像路径，**留空则直接调用本机 HiC-Pro**（默认行为）；`--plothic-path` 是画图软件 PlotHiC 的路径。**一般用默认即可**。

### 流程控制 | Control

**通俗理解|In plain words:** `--force` 强制重跑所有步骤（忽略断点续传）；`--verbose`/`--quiet` 控制日志详略。**换参数重跑时用 --force，日常重跑直接跑（自动跳过已完成的）**。

## 分析流程 | Pipeline

```text
基因组FASTA + R1/R2 FASTQ
    │
    ▼
步骤1: 准备 HiCPro 配置
    ├─ samtools faidx → 染色体大小文件
    ├─ 按限制酶扫描 → 酶切片段文件
    └─ 写 hicpro.conf（路径/线程/内存/酶/比对/过滤/bin 等）
    │
    ▼
步骤2: 运行 HiCPro（断点续传：矩阵已有效则跳过）
    ├─ 无 bowtie2 索引则自动构建
    ├─ 复制 FASTQ 到样本目录
    └─ HiC-Pro 比对 + 过滤 + 建矩阵（raw / iced）
    │
    ▼
步骤3: PlotHiC 画热图（断点续传：热图已存在则跳过）
    ├─ 挑最接近 resolution 的 matrix + abs.bed
    └─ plothic 出图并重命名为 {genome_id}_hic_heatmap.{format}
```

## 输出 | Output

```text
output/
├── hicpro.conf                     # 自动生成的 HiCPro 配置文件
├── {genome_id}.chrom.sizes         # 染色体大小表
├── {genome_id}_{酶}_resfrag.bed    # 酶切片段表（BED）
├── hicpro_output/
│   └── hic_results/matrix/{样本}/{raw,iced}/{bin_size}/…   # HiCPro 矩阵
├── plot/
│   └── {genome_id}_hic_heatmap.{format}   # 最终全基因组热图
└── hic_heatmap.log                 # 运行日志
```

## 结果解读 | Interpreting Results

### 1. 全基因组热图（`{genome_id}_hic_heatmap.pdf`）

**通俗理解|In plain words:** 这是核心结果。横纵轴都是基因组位置，颜色越深（红）表示两个位置互作越强。沿对角线的方块 = 每条染色体自身。

- 对角线亮、方块分明：每条染色体内部互作强、染色体间区隔明显（Hi-C 质量好）
- 对角线外（方块之间）也有明显色块：染色体间互作或基因组重排
- 全图发灰、对角线不清晰：数据质量差或分辨率选得过细

### 2. HiCPro 矩阵（`hicpro_output/…`）

**通俗理解|In plain words:** 中间产物，`raw` 是原始矩阵，`iced` 是归一化后的矩阵。PlotHiC 用的就是 iced 矩阵。普通用户无需深究，出现问题时看这里是否有 `.matrix` 和 `.bed` 文件。

### 3. 染色体大小与酶切片段

**通俗理解|In plain words:** 流程自动生成的"底图"文件，染色体大小表用于画图坐标系，酶切片段表记录限制酶的切割位置。核对热图染色体数是否正确时可看 chrom.sizes。

## 参数选择建议 | Parameter Guidance

- `--resolution`：全基因组起步用默认 100kb；局部精细分析再调小（需 bin-sizes 里含对应值）
- `--bin-sizes`：默认 "20000 40000 150000 500000 1000000" 已覆盖常用分辨率，一般不改
- `--restriction-enzyme`：务必与实验一致，否则酶切片段错误、结果不可信
- `-t/--max-memory`：按机器配置调；内存不足时降线程或降 max-memory
- `--force`：只在换参数后重跑时用，日常重跑保持默认（自动续传）

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

- HiC-Pro v3.1.0（conda 环境 `HiC-Pro_v3.1.0`，含 `bowtie2-build`）
- samtools（生成染色体大小）
- BioPython（`Bio.SeqIO`，生成酶切片段）
- PlotHiC（conda 环境 `plothic_v.1.0.0`，默认 `~/miniforge3/envs/plothic_v.1.0.0/bin/plothic`）
- 可选：Singularity（仅 `--hicpro-sif` 指定镜像时使用）

## 常见问题 | FAQ

**Q1：换参数重跑结果没变？**
断点续传按各步产物是否存在判断。改了 `--restriction-enzyme`/`--bin-sizes` 等参数后要加 `--force` 或删掉旧产物，否则复用旧结果。

**Q2：bowtie2 索引会自动构建吗？**
会。找不到 `.1.bt2` 索引时自动在输出目录用 `bowtie2-build` 构建（可能耗时较长）。也可用 `--bowtie2-idx` 指定现成索引前缀。

**Q3：报 HiC-Pro 不存在？**
检查 `~/software/HiC-Pro_v3.1.0/HiC-Pro_3.1.0/bin/HiC-Pro` 或通过环境变量 HICPRO_PATH / 配置文件指定路径；也可用 `--hicpro-sif` 走 Singularity 镜像。

**Q4：为什么热图只画了前几条染色体？**
HiCPro 阶段按基因组全部序列建矩阵，PlotHiC 按 chrom.sizes 全画。若结果异常，先核对 chrom.sizes 里染色体名是否与 FASTA 一致。

**Q5：矩阵已存在但被判无效而重跑？**
程序会检查 matrix 文件是否非空；若上次运行中断留下空文件，会自动重跑该步，属正常保护行为。
