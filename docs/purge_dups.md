# Purge_Dups 基因组去冗余 | Purge_Dups Genome Deduplication

一句话理解：**用测序深度当「尺子」，把基因组组装里重复出现的序列和单倍型冗余(haplotig)识别出来并去掉**，得到一套更干净的单倍体主组装。
输入一条组装 FASTA 和一套测序 reads，输出「去冗余后的主序列」和「被拆下来的冗余序列」两个 FASTA。

## 功能概述 | Overview

- 基于测序深度识别并去除组装中的单倍型(haplotig)与重叠(overlap)冗余，得到单倍体主组装
- 完整五步流程：计算深度 → 计算阈值 → 分割自比对 → 去冗余 → 取序列
- 支持 PacBio / HiFi 两种长读长 reads（`illumina` 类型当前未实现，见 FAQ）
- 可通过 `-s/--step 1~5` 单独运行某一步，方便中途失败后只补跑后半段
- 默认只去除 contig 末端的冗余(`--ends-only`)，保留中间区域

## 快速开始 | Quick Start

```bash
biopytools purge_dups -i assembly.fa -r hifi_reads.fq.gz -o purged_output
```

最小输入：一条组装 FASTA + 一套 HiFi/PacBio reads，输出目录里得到去冗余后的 `seqs/*_purged.purge.fa`。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| 单倍体组装 | 基因组只保留「一套」，像只留一副牌而不是两副 |
| 单倍型冗余(haplotig) | 同一段区域的两套版本都被装出来了，工具要丢掉多余的一套 |
| 重叠(overlap) | 两条 contig 首尾是同一段序列，重复算了两次 |
| 测序深度 | 每个位置被测到多少次；真实单拷贝区域和重复区域深度会不一样 |
| 阈值(cutoffs) | 一条分界线：深度明显偏高/偏低都可能是冗余，用来判定留不留 |
| 自比对(self-alignment) | 把组装序列跟自己比，找出哪里和哪里是「同一段」 |
| PAF | 比对结果的一种文本格式，记录两段序列哪里对上了 |

## 输入 | Input

### 组装 FASTA

待去冗余的基因组组装，标准 FASTA 格式。序列 ID 会用来给输出文件命名（去掉 `.fa`/`.fasta` 及 `.primary`/`.asm` 后缀后的名字作为 `{genome}` 前缀）。

### 测序 reads

- **PacBio**（`--read-type pacbio`）：连续长读长，内部用 minimap2 `-xmap-pb` 比对
- **HiFi**（`--read-type hifi`，默认）：高精度长读长，内部用 minimap2 `-xmap-hifi` 比对
- **Illumina**（`--read-type illumina`）：当前**未实现**，代码里会直接报错提示需先自行 bwa 比对（见 FAQ）

## 参数说明 | Parameters

### 必需参数与数据类型 | Required & read type

**通俗理解|In plain words:** `--read-type` 告诉程序你的 reads 是哪种，好选对子比对参数。HiFi 数据（现在的 PacBio 高精度模式）保持默认 `hifi` 即可；老一代 CLR 连续读长才用 `pacbio`。选错会导致深度统计和去冗余不准。

### purge_dups 核心阈值 | Core thresholds

**通俗理解|In plain words:** `--min-fraction` 是「一段序列要被判定为冗余，覆盖它的深度占比至少要多少」，默认 0.8，**一般不用动**；`--two-round-chaining` 是两轮链式匹配，能更细地找出短冗余，默认开启，**一般不用动**。这些都是 purge_dups 原生的精细参数，无生信基础的用户保持默认即可。

### get_seqs 输出选项 | Output options

**通俗理解|In plain words:** `--ends-only` 默认开启，表示「只把 contig 两端的冗余裁掉」，这是最常用、最安全的用法；`--no-ends-only` 会连 contig 中间的冗余一起处理，更激进，可能多删。`--min-primary-length` 是「短于这个长度的主序列不值得保留」，默认 10000bp，**一般不用动**。`--manual-cutoffs` 供想手工指定深度阈值文件的高级用户使用。

### 步骤控制 | Step control

**通俗理解|In plain words:** 正常跑全部五步即可。若某一步失败想只补跑后半段，用 `-s/--step 1~5` 单独指定；`--split-by-n` 让分割时按 N 缺口切，仅在组装里 N 缺口很多导致分割过长时才需要开。

## 分析流程 | Pipeline

```text
组装 FASTA + reads
    │
    ▼
步骤1: 计算深度(minimap2 比对 → pbcstat/ngscstat 统计)  → coverage/*_base.cov, *.stat
    │
    ▼
步骤2: 计算阈值(calcuts)                                  → coverage/cutoffs
    │
    ▼
步骤3: 分割基因组并自比对(split_fa + minimap2 -xasm5)      → split_aln/*_split_self_paf.gz
    │
    ▼
步骤4: 去冗余(purge_dups, 输出冗余区间)                     → purge_dups/dups.bed
    │
    ▼
步骤5: 取序列(get_seqs, 拆出主序列和冗余序列)               → seqs/*_purged.purge.fa + *.haplotigs.fa
```

## 输出 | Output

```text
purged_output/
├── coverage/
│   ├── {genome}_base.cov          # 每个碱基的测序深度
│   ├── {genome}.stat              # 深度统计(用于算阈值)
│   └── cutoffs                    # 深度分界阈值
├── split_aln/
│   ├── {genome}.split             # 分割后的序列
│   └── {genome}_split_self_paf.gz # 自比对结果
├── purge_dups/
│   ├── dups.bed                   # 被判定为冗余的区间清单
│   └── purge_dups.log             # purge_dups 子命令日志
├── seqs/
│   ├── {genome}_purged.purge.fa      # 去冗余后的主组装(最终产物)
│   └── {genome}_purged.haplotigs.fa  # 被拆下来的冗余/单倍型序列
└── purge_dups.log                 # 流程主日志
```

- `seqs/{genome}_purged.purge.fa`：最终主组装，拿去下游分析的就是它
- `seqs/{genome}_purged.haplotigs.fa`：被判定为冗余的单倍型序列，别直接删，先看两眼确认没误删
- `purge_dups/dups.bed`：冗余区间坐标，想复核「哪些位置被判为冗余」就查这个

## 结果解读 | Interpreting Results

- **主序列变短是正常现象**：去冗余必然让总长下降（去掉了多装的一份）。关注的是「主序列是否还保留了预期大小」，而不是和原始组装一样大
- **haplotigs.fa 里不该有核心单拷贝基因**：如果发现关键基因被丢进了冗余文件，说明阈值偏严，考虑调 `--manual-cutoffs` 或检查 `--read-type` 是否选对
- **dups.bed 区间应集中在 contig 末端**（默认 `--ends-only`）：若中间大量被判冗余，多为数据或参数问题

## 参数选择建议 | Parameter Guidance

- **HiFi 数据**：默认参数直接跑，最省心
- **PacBio CLR 老数据**：`--read-type pacbio` 必改，否则比对和深度统计不对
- **只想裁末端、保守去冗余**：保持默认 `--ends-only`
- **中间也有冗余、想更彻底**：`--no-ends-only`
- **中途失败只补后半段**：`-s 3` 从分割比对重跑，或 `-s 5` 只重跑取序列
- **阈值想自己定**：先用默认跑一遍拿到 `coverage/*.stat`，再据此写 `--manual-cutoffs` 文件

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input, -i` | 必填 |  | 基因组组装文件(FASTA格式)｜Genome assembly file in FASTA format |
| `--reads, -r` | 必填 |  | 测序文件(PacBio/HiFi/illumina)｜Sequencing reads file |
| `--output-dir, -o` | `./purge_dups_output` | Path | 输出目录｜Output directory |
| `--threads, -t` | `12` | int | 线程数｜Number of threads |
| `--read-type` | `hifi` | pacbio/hifi/illumina | 测序数据类型｜Sequencing data type |
| `--min-fraction` | `0.8` | float | 最小比例阈值｜Minimum fraction threshold |
| `--two-round-chaining` | `True` |  | 启用两轮链式匹配｜Enable two-round chaining |
| `--no-two-round-chaining` | — |  | 禁用两轮链式匹配｜Disable two-round chaining |
| `--ends-only` | `True` |  | 只去除contig末端的冗余｜Only remove duplications at contig ends |
| `--no-ends-only` | — |  | 也去除contig中间的冗余｜Also remove duplications in contig middle |
| `--min-primary-length` | `10000` | int | 最小主contig长度｜Minimum primary contig length |
| `--step, -s` | — | 1/2/3/4/5 | 运行指定步骤｜Run only specified step: 1: 计算深度｜Calculate coverage 2: 计算阈值｜Calculate cutoffs 3: 分割比对｜Split and align 4: 去冗余｜Purge duplications 5: 获取序列｜Get sequences |
| `--split-by-n` | — |  | split_fa按N分割｜split_fa split by N |
| `--manual-cutoffs` | — | Path | 手动指定阈值文件｜Manual cutoffs file |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 基因组组装文件(FASTA格式)｜Genome assembly file in FASTA format |
| `-r, --reads` | 必填 |  | 测序文件(PacBio/HiFi/illumina)｜Sequencing reads file |
| `-o, --output-dir` | `./purge_dups_output` |  | 输出目录｜Output directory |
| `-t, --threads` | `12` | int | 线程数｜Number of threads |
| `--read-type` | `hifi` | pacbio/hifi/illumina | 测序数据类型｜Sequencing data type |
| `--min-fraction` | `0.8` | float | 最小比例阈值｜Minimum fraction threshold |
| `--two-round-chaining` | `True` | store_true | 启用两轮链式匹配｜Enable two-round chaining |
| `--no-two-round-chaining` | — | store_false | 禁用两轮链式匹配｜Disable two-round chaining |
| `--ends-only` | `True` | store_true | 只去除contig末端的冗余｜Only remove duplications at contig ends |
| `--no-ends-only` | — | store_false | 也去除contig中间的冗余｜Also remove duplications in contig middle |
| `--min-primary-length` | `10000` | int | 最小主contig长度｜Minimum primary contig length |
| `-s, --step` | — | 1/2/3/4/5 | 只运行指定步骤｜Run only specified step (1: 计算深度｜coverage, 2: 计算阈值｜cutoffs, 3: 分割比对｜split&align, 4: 去冗余｜purge, 5: 获取序列｜get seqs) |
| `--split-by-n` | — | store_true | split_fa按N分割｜split_fa split by N |
| `--manual-cutoffs` | — |  | 手动指定阈值文件｜Manual cutoffs file |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- conda 环境 `purge_dups_v.1.2.6`（默认路径 `~/miniforge3/envs/purge_dups_v.1.2.6`），内含 purge_dups 全套工具：`purge_dups`、`split_fa`、`get_seqs`、`calcuts`、`pbcstat`、`ngscstat`、`minimap2`
- 无需额外单独安装其它软件

## 常见问题 | FAQ

**Q1：`--read-type illumina` 能用吗？**
不能。代码里 Illumina 分支未实现，会报错提示「需要先 bwa 比对生成 BAM」。Illumina 数据请改用其它去冗余工具，或自行走 purge_dups 的 ngscstat 流程。

**Q2：重跑会不会跳过已完成的步骤？**
不会。本模块**没有断点续传**——完整重跑会从头再做一遍并覆盖旧产物。想只补跑某步用 `-s/--step` 单独指定。

**Q3：换参数重跑，结果怎么没变？**
因为没有跳过逻辑、旧产物位置固定，换参数重跑容易造成混淆。建议换参数时**换一个新的输出目录**，保证结果干净、可对比。

**Q4：主序列比原来短了很多，正常吗？**
正常。去冗余就是去掉「多装的那一份」，总长必然下降。重点看关键基因是否还在主序列里，而不是追求和原来一样大。

**Q5：输出文件名里的 `{genome}` 是怎么来的？**
取输入 FASTA 的文件名（去掉 `.fa`/`.fasta`，再剥掉 `.primary`/`.asm` 后缀）。例如 `sample.primary.fa` 会得到 `sample` 前缀。
