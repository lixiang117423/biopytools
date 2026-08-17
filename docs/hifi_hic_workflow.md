# HiFi+Hi-C 组装挂载工作流 | HiFi+Hi-C Assembly and Scaffolding Workflow

一句话理解：**一条龙把植物基因组从 HiFi reads 拼到染色体级**——HiFi 组装 → Hi-C 挂载 → 染色体重命名 → Hi-C 热图，中间只需给一份参考基因组做「命名对照」。

## 功能概述 | Overview

- 四步串联：HiFi 组装（复用 hifi_hic）→ HapHiC 挂载 → 染色体重命名 → Hi-C 热图
- 参考基因组仅用于「染色体命名」，不参与组装
- 断点续传默认启用，支持 `--skip-*` 跳过单步、`--force` 强制重跑
- 可选 NGS polish
- 染色体数（nchrs）默认从参考基因组自动统计

## 快速开始 | Quick Start

```bash
biopytools hifi-hic-workflow --hifi hifi.fq.gz --hic-r1 R1.fq.gz --hic-r2 R2.fq.gz --ref reference.fa -o workflow_output
```

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| Hi-C 挂载(scaffolding) | 用 Hi-C 的空间邻近信息，把零散的 contig「按正确顺序排成染色体」 |
| scaffold | 用 gap 把多条 contig 串起来的更长的序列 |
| 染色体级基因组 | 每条染色体拼成一条（或接近一条）的最终基因组 |
| 染色体重命名 | 把组装出的序列按参考基因组命名成 chr01/chr02…，方便下游 |
| Hi-C 热图 | 用色块展示序列两两之间的互作强度；对角线清晰、两侧噪音少 = 挂载正确 |
| 参考基因组(仅命名) | 近缘物种已发表的基因组，只用来「对照起名」，不参与拼装 |

## 输入 | Input

### 必需输入

- HiFi reads：`--hifi`
- Hi-C 双端：`--hic-r1` / `--hic-r2`
- 参考基因组 FASTA：`--ref`（仅用于染色体命名）

### 可选输入

- NGS 目录：`--ngs-data`（配合 `--use-ngs-polish`）

## 参数说明 | Parameters

### 必需参数与全局参数 | Required & global

**通俗理解|In plain words:** `--hifi`/`--hic-r1`/`--hic-r2`/`--ref`/`-o` 是五个必需输入输出；`-p` 前缀、`-t` 线程数（默认 64）。

### 流程控制 | Workflow control

**通俗理解|In plain words:** 四个 `--skip-*` 可以只跑/跳过某一步（比如只重跑热图）；`--no-resume` 关断点续传，`--force` 强制重跑。**正常首次运行一个都不用加**。

### HiFi 组装与 NGS polish | HiFi assembly & NGS polish

**通俗理解|In plain words:** 传给底层 hifi_hic 的组装参数。`--genome-size` 报预算（默认 1.45g）、`--n-hap` 倍性（默认 2）、`--purge-level`/`--hom-cov` 一般不用动；`--use-ngs-polish` + `--ngs-data` 开启二代纠错。

### HapHiC 挂载 | HapHiC scaffolding

**通俗理解|In plain words:** `--nchrs` 是染色体条数，不指定会从参考基因组自动统计（一般不用管）。`--haphic-bin`/`--bwa-bin`/`--samtools-bin` 是工具路径，一般不用动。

### 染色体重命名 | Chromosome rename

**通俗理解|In plain words:** 按参考基因组给序列命名。`--naming-min-identity`/`--naming-min-coverage` 是「和参考比对的相似度/覆盖度门槛」，低于门槛的序列不归到对应染色体；**默认 80 一般够用**。`--rename-keep-all` 默认保留所有序列（chr + scaffold）。

### Hi-C 热图 | Hi-C heatmap

**通俗理解|In plain words:** 热图的画法。`--heatmap-resolution` 是分辨率（bp，默认 10 万），越大越粗糙但文件小；`--heatmap-format` 输出格式（默认 pdf）；`--hicpro-enzyme` 是 Hi-C 建库用的限制性内切酶（默认 MboI），**务必和实际建库酶一致**。

## 分析流程 | Pipeline

```text
HiFi reads + Hi-C R1/R2 + 参考基因组(仅命名)
    │
    ▼
步骤1: HiFi 组装(复用 hifi_hic) → 01_hifi_assembly/{prefix}/02_fasta/{prefix}.primary.fa
    │
    ▼
步骤2: HapHiC 挂载 → 02_hic_scaffolding/04_build/{prefix}.fa
    │
    ▼
步骤3: 染色体重命名 → 03_chromosome_rename/{prefix}.renamed.fa
    │
    ▼
步骤4: Hi-C 热图 → 04_hic_heatmap/plot/{prefix}_hic_heatmap.pdf
```

## 输出 | Output

```text
workflow_output/
├── 00_pipeline_info/                   # 各步骤信息文件({step}_info.txt)
├── 01_hifi_assembly/
│   └── {prefix}/02_fasta/{prefix}.primary.fa    # HiFi 组装结果
├── 02_hic_scaffolding/
│   └── 04_build/{prefix}.fa                     # Hi-C 挂载后的 scaffold
├── 03_chromosome_rename/
│   └── {prefix}.renamed.fa                      # 重命名后的染色体级基因组(最终)
├── 04_hic_heatmap/
│   └── plot/{prefix}_hic_heatmap.pdf            # Hi-C 热图
├── logs/
│   └── workflow_时间戳.log                      # 工作流日志
└── workflow_report.txt                          # 流程报告(步骤状态+输出清单)
```

## 结果解读 | Interpreting Results

**通俗理解|In plain words:** 最终基因组是 `03_chromosome_rename/{prefix}.renamed.fa`；热图用来「验收」挂载是否正确。

- **`{prefix}.renamed.fa`**：最终染色体级基因组，序列名应为 chr01/chr02…，可直接下游使用
- **序列数 vs 染色体数**：重命名后染色体条数应接近 `--nchrs`（自动统计值）；若还残留大量 scaffold，说明挂载未完全
- **Hi-C 热图**：理想情况是对角线清晰、对角线两侧（近距离互作）有信号、远距离噪音少；若出现大片「十字」或非对角线块，提示有 mis-join（挂错）或易位
- **`workflow_report.txt`**：汇总每一步完成状态与主要输出路径，先看它快速定位哪步出问题

## 参数选择建议 | Parameter Guidance

- `--nchrs`：一般不用指定，自动从参考基因组统计；若参考不完整可手动给定
- `--genome-size`：报预估大小（宁大勿小）
- `--hicpro-enzyme`：务必与 Hi-C 建库实际使用的酶一致（MboI/DpnII/HindIII 等），否则热图/挂载可能失真
- `--skip-haphic --skip-rename`：只想重新跑热图时，可跳过前面步骤（需已有对应产物）
- `--heatmap-resolution`：默认 100000 bp；想看更细结构可调小（如 50000），文件会变大

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--hifi` | 必填 | Path | HiFi reads文件｜HiFi reads file |
| `--hic-r1` | 必填 | Path | Hi-C R1文件｜Hi-C R1 file |
| `--hic-r2` | 必填 | Path | Hi-C R2文件｜Hi-C R2 file |
| `--ref, --reference` | 必填 | Path | 参考基因组FASTA文件（仅用于命名）｜Reference genome FASTA file (for naming only) |
| `-o, --output` | 必填 | Path | 输出目录｜Output directory |
| `-p, --prefix` | `genome_sample` |  | 样本前缀｜Sample prefix (default: genome_sample) |
| `-t, --threads` | `64` | int | 线程数｜Number of threads (default: 64) |
| `--use-ngs-polish` | — |  | 启用NGS polish｜Enable NGS polish |
| `--ngs-data` | — | Path | NGS二代数据目录｜NGS second-generation data directory |
| `--nchrs` | — | int | 染色体数量（如不指定，从reference统计）｜Number of chromosomes (count from reference if not specified) |
| `--skip-hifi-hic` | — |  | 跳过HiFi组装｜Skip HiFi assembly |
| `--skip-haphic` | — |  | 跳过Hi-C挂载｜Skip Hi-C scaffolding |
| `--skip-rename` | — |  | 跳过重命名｜Skip renaming |
| `--skip-heatmap` | — |  | 跳过热图｜Skip heatmap |
| `--no-resume` | — |  | 禁用断点续传｜Disable resume mode |
| `--force` | — |  | 强制重新运行所有步骤｜Force rerun all steps |
| `-v, --verbose` | — |  | 显示详细日志｜Show verbose logs |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--hifi` | 必填 |  | HiFi reads文件｜HiFi reads file |
| `--hic-r1` | 必填 |  | Hi-C R1文件｜Hi-C R1 file |
| `--hic-r2` | 必填 |  | Hi-C R2文件｜Hi-C R2 file |
| `--ref, --reference` | 必填 |  | 参考基因组FASTA文件（仅用于命名）｜Reference genome FASTA file (for naming only) |
| `-o, --output` | 必填 |  | 输出目录｜Output directory |
| `-p, --prefix` | `genome_sample` |  | 样本前缀｜Sample prefix (default: genome_sample) |
| `-t, --threads` | `64` | int | 线程数｜Number of threads (default: 64) |
| `--skip-hifi-hic` | — | store_true | 跳过HiFi组装｜Skip HiFi assembly |
| `--skip-haphic` | — | store_true | 跳过Hi-C挂载｜Skip Hi-C scaffolding |
| `--skip-rename` | — | store_true | 跳过重命名｜Skip renaming |
| `--skip-heatmap` | — | store_true | 跳过热图｜Skip heatmap |
| `--no-resume` | — | store_true | 禁用断点续传｜Disable resume mode |
| `--force` | — | store_true | 强制重新运行所有步骤｜Force rerun all steps |
| `--genome-size` | `1.45g` |  | 预估基因组大小｜Estimated genome size (default: 1.45g) |
| `--n-hap` | `2` | int | 倍性｜Ploidy (default: 2) |
| `--purge-level` | — | int | Purge level (0=no purging, 1=light, 2/3=aggressive) |
| `--hom-cov` | — | int | Homozygous read coverage (--hom-cov) |
| `--use-ngs-polish` | — | store_true | 启用NGS polish｜Enable NGS polish |
| `--ngs-data` | — |  | NGS二代数据目录｜NGS second-generation data directory |
| `--ngs-high-cov` | `95.0` | float | 高质量contig覆盖度阈值｜High quality contig coverage threshold (default: 95.0) |
| `--ngs-pattern` | `_1.clean.fq.gz` |  | NGS文件匹配模式｜NGS file matching pattern (default: _1.clean.fq.gz) |
| `--nchrs` | — | int | 染色体数量（如不指定，从reference统计）｜Number of chromosomes (count from reference if not specified) |
| `--haphic-bin` | — |  | HapHiC可执行文件路径｜HapHiC executable path |
| `--bwa-bin` | — |  | BWA可执行文件路径｜BWA executable path |
| `--samtools-bin` | — |  | Samtools可执行文件路径｜Samtools executable path |
| `--rename-keep-all` | `True` | store_true | 保留所有序列（chrNN + scaffolds）｜Keep all sequences (default: True) |
| `--naming-min-identity` | `80.0` | float | 最小序列一致性｜Min sequence identity %% (default: 80.0) |
| `--naming-min-coverage` | `80.0` | float | 最小覆盖度｜Min coverage %% (default: 80.0) |
| `--naming-minimap2-preset` | `asm5` | asm5/asm10/asm20 | minimap2预设｜minimap2 preset (default: asm5) |
| `--hicpro-enzyme` | `MboI` |  | HiCPro限制性内切酶｜HiCPro restriction enzyme (default: MboI) |
| `--heatmap-resolution` | `100000` | int | 热图分辨率｜Heatmap resolution in bp (default: 100000) |
| `--heatmap-colormap` | `YlOrRd` |  | 热图颜色方案｜Heatmap color scheme (default: YlOrRd) |
| `--heatmap-format` | `pdf` |  | 热图输出格式｜Heatmap output format (default: pdf) |
| `-v, --verbose` | — | store_true | 显示详细日志｜Show verbose logs |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

本工作流复用以下模块与软件：

- hifi_hic 模块（依赖 hifiasm `asm`、seqkit `misc`、samtools `align`）
- HapHiC（`--haphic-bin`，默认 `haphic`）、bwa、samtools（Hi-C 比对与挂载）
- minimap2（染色体重命名时与参考比对）
- HiC-Pro（默认路径 `~/software/HiC-Pro_v3.1.0/.../HiC-Pro`）、PlotHiC（conda 环境 `plothic_v.1.0.0`）

## 常见问题 | FAQ

**Q1：支持断点续传吗？**
支持，默认启用。每步完成后写 `00_pipeline_info/{step}_info.txt` 并按产物存在性跳过；`--no-resume` 关闭，`--force` 强制全部重跑。只重跑某一步可用 `--skip-*` 组合。

**Q2：参考基因组是拿来干嘛的？**
只用于「染色体命名」（对照近缘参考给 chr01/chr02 起名）和自动统计染色体条数，不参与组装，也不会把你的组装「替换」成参考序列。

**Q3：染色体数怎么来的？**
默认从 `--ref` 的 FASTA 里统计 `>` 开头的序列条数；若参考是未去冗余的（含 scaffold），可手动用 `--nchrs` 指定真实染色体数。

**Q4：热图看起来乱怎么办？**
先确认 `--hicpro-enzyme` 与建库酶一致；再看 `--nchrs` 是否合理；分辨率过低也可能掩盖结构。可尝试 `--skip-haphic --skip-rename --force` 单独重跑热图。

**Q5：只想做前两步（组装+挂载），不要热图？**
加 `--skip-heatmap` 即可；同理 `--skip-rename` 跳过重命名。

