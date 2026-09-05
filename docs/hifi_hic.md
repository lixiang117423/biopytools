# HiFi+Hi-C 基因组组装 | HiFi+Hi-C Genome Assembly (hifi-hic)

一句话理解：**用 HiFi 数据（可加 Hi-C）组装基因组，再自动做去冗余**。输入一份 HiFi reads（可选 Hi-C、可选 NGS），输出主组装、两套单倍型，以及去冗余后的最终基因组。

## 功能概述 | Overview

- 基于 hifiasm 组装，支持「仅 HiFi」或「HiFi + Hi-C」两种模式
- 可选 NGS polish：用二代数据校正组装
- 默认启用 Purge_Dups 去冗余，去掉因杂合/重复导致的冗余序列
- 断点续传默认启用，按各步骤产物存在性自动跳过
- 按倍性（n-hap）动态输出 primary / hap1 / hap2 等 FASTA

## 快速开始 | Quick Start

```bash
biopytools hifi-hic -i hifi.fq -p sample1
```

最小输入：一份 HiFi reads（FASTQ/FASTA）；加 `-1`/`-2` 可附带 Hi-C，加 `--ngs` 可做 NGS polish。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| HiFi reads | PacBio 高保真长读长 |
| Hi-C | 一种能反映 DNA 在空间上「谁挨着谁」的测序，可辅助把 contig 挂到染色体级 |
| 单倍型(haplotype) | 二倍体来自父/母的两套基因组，hifiasm 可分开输出 hap1/hap2 |
| primary / alternate | 主组装（主要那套）与备选组装（另一套） |
| 倍性(n-hap) | 一套基因组的染色体组数；二倍体=2 |
| 去冗余(Purge_Dups) | 把因杂合或重复而被「多拼一份」的序列去掉，让结果更接近真实一套基因组 |
| NGS polish | 用精确的二代短 reads 校正长读组装的碱基错误 |
| 覆盖度(coverage) | 某条 contig 被 reads 覆盖的次数，用于区分「真序列」和「冗余/污染」 |

## 输入 | Input

### HiFi reads

HiFi 测序数据，FASTQ/FASTA：

```text
@read1
ACGTACGT...
+
IIIIIIII...
```

### 可选输入

- Hi-C：`--hic-r1` / `--hic-r2`（两端成对）
- NGS：`--ngs` 指定二代数据目录，`--ngs-pattern` 指定 R1 文件匹配模式（默认 `_1_clean.fq.gz`）

## 参数说明 | Parameters

### 必需参数与基本参数 | Required & basic

**通俗理解|In plain words:** `-i` 是 HiFi 数据，`-p` 是前缀。`-g` 基因组大小、`--n-hap` 倍性、`-t` 线程数（默认 88，偏大，按机器调整）。

### 组装参数 | Assembly parameters

**通俗理解|In plain words:** `--purge-level` 控制 hifiasm 内部对重复/冗余的清理力度（不指定时 hifiasm 自动选），`--hom-cov` 指定纯合覆盖度（不指定时自动）。**这两个一般不用动**，只有对高杂合/特殊倍性物种才需要手调。

### 去冗余参数 | Purge_Dups

**通俗理解|In plain words:** 默认会自动做 Purge_Dups 去冗余。`--high-cov`/`--medium-cov-min` 决定「多少覆盖算真、多少算冗余」，**一般不用动**；如果只想直接要 hifiasm 原始结果、不做去冗余，用 `--no-purge-dups` 关掉。

### NGS polish 与执行控制 | NGS polish & execution

**通俗理解|In plain words:** `--ngs` 开启 NGS 纠错；`--no-resume` 关闭断点续传（强制全部重跑）。

## 分析流程 | Pipeline

```text
HiFi reads (+ Hi-C 可选)
    │
    ▼
步骤1: hifiasm 组装(01_raw_output, GFA)
    │
    ▼
步骤2: GFA → FASTA 转换(02_fasta, primary/hap1/hap2/alternate)
    │
    ▼
步骤3: 生成 contig-reads 映射(02_fasta/*_contig_reads.tsv)
    │
    ▼
步骤4: NGS polish(可选, 03_ngs_polish)
    │
    ▼
步骤5: Purge_Dups 去冗余(默认, 04_purge_dups → *_purged.purge.fa)
```

## 输出 | Output

```text
assembly_output/
└── {prefix}/
    ├── 00_pipeline_info/
    │   └── software_versions.yml        # 软件版本记录
    ├── 01_raw_output/                   # hifiasm 原始 GFA
    │   ├── {prefix}.hic.p_ctg.gfa       # (Hi-C 模式) primary contigs
    │   ├── {prefix}.hic.hap1.p_ctg.gfa  # 单倍型1
    │   ├── {prefix}.hic.hap2.p_ctg.gfa  # 单倍型2
    │   └── {prefix}.hic.a_ctg.gfa       # alternate
    ├── 02_fasta/                        # 转换后的 FASTA
    │   ├── {prefix}_primary.fa          # 主组装(最常用)
    │   ├── {prefix}_hap1.fa / {prefix}_hap2.fa / {prefix}_alternate.fa
    │   └── {prefix}_p_ctg_contig_reads.tsv  # contig→reads 映射
    ├── 03_ngs_polish/                   # (仅给 --ngs 时)
    │   └── {prefix}_polished.fa         # NGS 校正后的基因组
    ├── 04_purge_dups/                   # 去冗余(默认)
    │   └── sequences/{prefix}_purged.purge.fa  # 去冗余最终结果
    └── 99_logs/
        └── hifiasm_assembly.log         # 组装日志
```

> 注：仅 HiFi 模式（不给 Hi-C）时，GFA/FASTA 文件名为 `{prefix}.p_ctg.gfa`、`{prefix}.hap1.p_ctg.gfa` 等（无 `.hic.` 段）。hifiasm 使用 purge 时文件名可能带 `.bp.` 前缀，模块会自动识别两种命名。

## 结果解读 | Interpreting Results

**通俗理解|In plain words:** 先看 `04_purge_dups/sequences/{prefix}_purged.purge.fa`（去冗余后的最终基因组），再看 `02_fasta/{prefix}_primary.fa`（主组装）。

- **去冗余后的 purged.fa**：最终推荐使用的基因组，冗余少、更接近真实一套
- **primary.fa vs hap1/hap2.fa**：primary 是主组装；hap1/hap2 是拆分出的两套单倍型，供研究等位差异用
- **序列数/长度**：去冗余后序列数应明显下降（冗余被合并/删除）；若几乎没变化，说明原始组装本身杂合/冗余低
- `software_versions.yml` 记录 hifiasm/seqkit/samtools 的版本，写论文 Methods 时直接抄

## 参数选择建议 | Parameter Guidance

- `--n-hap`：二倍体=2（默认）；单倍体=1；多倍体按染色体组数
- `--genome-size`：**不传则 hifiasm 自动估计**（推荐，`--hg-size auto`）；手动指定时报预估大小（宁大勿小），单位 `g`/`m`——早期版本写死默认 `1.45g`，对小于/大于 1.45Gb 的基因组会跑偏覆盖度推断，已移除
- `--no-purge-dups`：想要 hifiasm 原始结果、或后续自己控制去冗余时关掉
- `--threads`：默认 88 偏大，按机器核数调整；去冗余默认复用组装线程数（`--purge-dups-threads` 可单独指定）
- `--high-cov`/`--medium-cov-min`：**一般不用动**，只有覆盖分布异常（如极高覆盖污染）时才调

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--hifi, -i` | 必填 | Path | HiFi数据文件路径｜Path to HiFi data file |
| `--hic-r1, -1` | — | Path | Hi-C R1文件路径（可选）｜Path to Hi-C R1 file (optional) |
| `--hic-r2, -2` | — | Path | Hi-C R2文件路径（可选）｜Path to Hi-C R2 file (optional) |
| `--prefix, -p` | `genome_sample` | str | 样本前缀｜Sample prefix |
| `--threads, -t` | `88` | int | 线程数｜Number of threads |
| `--genome-size, -g` | — | str | 预估基因组大小,不传则hifiasm自动估计｜Estimated genome size (e.g., 1.2g, 250m); omit for hifiasm auto |
| `--n-hap` | `2` | int | 倍性｜Ploidy (haploid count) |
| `--purge-level, -l` | — | int | Purge level (0=no purging, 1=light, 2/3=aggressive) [default: 3 for unzip] |
| `--hom-cov` | — | int | Homozygous read coverage (--hom-cov) [default: auto] |
| `--output, -o` | `./assembly_output` | Path | 输出目录｜Output directory |
| `--ngs` | — | Path | NGS二代数据目录（可选）｜NGS second-generation data directory (optional) |
| `--ngs-pattern` | `_1_clean.fq.gz` | str | NGS文件匹配模式｜NGS file matching pattern (default: _1_clean.fq.gz) |
| `--high-cov` | `95.0` | float | 高质量contig覆盖度阈值｜High quality contig coverage threshold (default: 95.0) |
| `--medium-cov-min` | `30.0` | float | 中等质量contig最小覆盖度｜Medium quality contig minimum coverage (default: 30.0) |
| `--no-purge-dups` | `False` |  | 禁用Purge_Dups去冗余｜Disable Purge_Dups deduplication (enabled by default) |
| `--purge-dups-path` | `~/miniforge3/envs/purge_dups_v.1.2.6` | str | Purge_Dups软件路径｜Purge_Dups software path (default: ~/miniforge3/envs/purge_dups_v.1.2.6) |
| `--purge-dups-threads` | — | int | 去冗余线程数｜Deduplication threads (default: same as assembly threads) |
| `--purge-dups-read-type` | `hifi` | pacbio/hifi/illumina | 去冗余reads类型｜Deduplication reads type (default: hifi) |
| `--no-resume` | `False` |  | 禁用断点续传（强制重新运行所有步骤）｜Disable resume mode (force rerun all steps) |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--hifi, -i` | 必填 |  | HiFi数据文件路径｜Path to HiFi data file |
| `--hic-r1, -1` | — |  | Hi-C R1文件路径｜Path to Hi-C R1 file |
| `--hic-r2, -2` | — |  | Hi-C R2文件路径｜Path to Hi-C R2 file |
| `--prefix, -p` | `genome_sample` |  | 样本前缀｜Sample prefix |
| `--threads, -t` | `88` | int | 线程数｜Number of threads |
| `--genome-size, -g` | — |  | 预估基因组大小,不传则hifiasm自动估计｜Estimated genome size (e.g., 1.2g, 250m); omit for hifiasm auto |
| `--n-hap` | `2` | int | 倍性｜Ploidy (haploid count) |
| `--purge-level, -l` | — | int | Purge level (0=no purging, 1=light, 2/3=aggressive) [default: 3 for unzip] |
| `--hom-cov` | — | int | Homozygous read coverage (--hom-cov) [default: auto] |
| `--output, -o` | `./assembly_output` |  | 输出目录｜Output directory |
| `--ngs` | — |  | NGS二代数据目录｜NGS second-generation data directory (optional) |
| `--ngs-pattern` | `_1_clean.fq.gz` |  | NGS文件匹配模式｜NGS file matching pattern (default: _1_clean.fq.gz) |
| `--high-cov` | `95.0` | float | 高质量contig覆盖度阈值｜High quality contig coverage threshold (default: 95.0) |
| `--medium-cov-min` | `30.0` | float | 中等质量contig最小覆盖度｜Medium quality contig minimum coverage (default: 30.0) |
| `--no-purge-dups` | — | store_true | 禁用Purge_Dups去冗余｜Disable Purge_Dups deduplication (enabled by default) |
| `--purge-dups-path` | `~/miniforge3/envs/purge_dups_v.1.2.6` |  | Purge_Dups软件路径｜Purge_Dups software path (default: ~/miniforge3/envs/purge_dups_v.1.2.6) |
| `--purge-dups-threads` | — | int | 去冗余线程数｜Deduplication threads (default: same as assembly threads) |
| `--purge-dups-read-type` | `hifi` | pacbio/hifi/illumina | 去冗余reads类型｜Deduplication reads type (default: hifi) |
| `--no-resume` | — | store_true | 禁用断点续传（强制重新运行所有步骤）｜Disable resume mode (force rerun all steps) |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- hifiasm（conda 环境 `asm`，默认路径 `~/miniforge3/envs/asm/bin/hifiasm`）
- seqkit（conda 环境 `misc`，用于 GFA→FASTA 序列格式化）
- samtools（conda 环境 `align`）
- Purge_Dups（默认目录 `~/miniforge3/envs/purge_dups_v.1.2.6`）
- bwa（仅 NGS polish 时，conda 环境 `align`）

## 常见问题 | FAQ

**Q1：支持断点续传吗？**
支持，默认启用。按 `primary.fa`、contig-reads 映射、NGS polish、去冗余各步骤产物是否存在来跳过；用 `--no-resume` 可强制全部重跑。换参数重跑前建议删除旧产物。

**Q2：为什么文件名里有时带 `.bp.`？**
hifiasm 使用 purge 时会在前缀后加 `.bp.`（如 `{prefix}.bp.hic.p_ctg.gfa`）。模块会自动识别并兼容两种命名，无需手动改。

**Q3：去冗余后结果和 primary.fa 有什么不同？**
Purge_Dups 会基于覆盖度把「多拼的一份」（杂合冗余）删掉，purged.fa 通常更短、序列更少、更接近真实单套基因组，推荐作为最终结果。

**Q4：只想组装、不想去冗余怎么办？**
加 `--no-purge-dups` 关闭去冗余，最终结果用 `02_fasta/{prefix}_primary.fa`。

**Q5：NGS polish 需要什么命名？**
`--ngs` 目录里的二代数据需按 `--ngs-pattern`（默认 `_1_clean.fq.gz`）匹配 R1，程序自动推导 R2。

