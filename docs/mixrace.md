# mixrace — WGS 混合小种检测 | WGS Mixed-race Detection

一句话理解：**输入一堆重测序数据和一个参考基因组，自动判断每个样品是「只有一个菌株」还是「多个菌株混在一起」，并给出一份能直接发给人看的判读报告**。

## 功能概述 | Overview

- 面向**单倍体病原**(如根肿菌 *Plasmodiophora brassicae*，静息孢子 n=20)——它们无法纯培养，一个田间编号常是多个基因型(小种)的混合
- 全流程：bwa-mem2 比对 → samtools markdup → **freebayes 单倍体 calling**(-p 1，保留 2% 低频等位)→ 过滤 → 等位频率谱(AFS)分析 → smudgescope k-mer 谱 → 判读 + 报告 + 系统发育树
- 核心产出：`summary/mixrace_report.html`(**单文件可分享**，图全部内嵌)，含每样品判读(疑似纯/疑似混合/不确定)、依据、指标通俗解释、样品聚类树

## 快速开始 | Quick Start

```bash
biopytools mixrace -i fastq_dir -g ref.fa -o out_dir/
```

(已做过 QC 的 clean 数据用 `--clean-fastq-dir` 替代 `-i`，可跳过 fastp。)

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗解释 |
|---|---|
| 单倍体(haploid) | 每个位置只有一份遗传信息。像一个班级里每个学生只举**一只手**——同时出现两只手(两种碱基)就是「混进了别的班的人」 |
| 杂合率(het_rate) | 「两只手同时举」的位置占整个基因组的比例。越低越纯；<0.01% 基本纯，>1% 明显不纯 |
| 等位频率谱(AFS) | 把所有「两只手」位置的少数派比例画成直方图。像看一个聚会的人数构成：只有一类人(纯)、两伙人各半(0.5 峰)、一主一从、大杂烩 |
| 优势株占比 | 若为混合，最大那伙人约占总人数的百分比 |
| 测序深度(depth) | 每个位置平均被读了**几遍**。像数人数时每个人被点了几次名；≥50x 才可靠 |
| 系统发育树 | 按基因组相似度给样品「排家谱」：分支越近越像 |
| VAF | 某个位置少数碱基的读数占比 |

## 输入 | Input

- `-i` 原始 fastq 目录(自动配对 `_1/_2`、`_R1/_R2`)或 `--clean-fastq-dir` 已清洗 fastq(二选一)。
- `-g` 参考基因组 FASTA(如 e3)。
- 可选 `--repeat-bed`：重复/低复杂度区域 BED，给出则过滤时排除，不给则跳过(不影响流程)。
- 可选 `--pure-samples`：已知纯样品，逗号分隔(如 `--pure-samples P1,P2`)，用于自动校准杂合率判读阈值。

## 参数说明 | Parameters

**通俗理解|In plain words:** 绝大多数参数**一般不用动**——默认值是针对本类分析调好的。只有判读阈值(想严一点/松一点)和 `--pure-samples`(手上有已知纯样品)才需要考虑。

| 参数 | 默认 | 说明 |
|---|---|---|
| `-i/--input` | — | 原始 fastq 目录(与 `--clean-fastq-dir` 二选一) |
| `--clean-fastq-dir` | — | 已清洗 fastq(跳过 QC) |
| `-g/--genome` | 必填 | 参考基因组 |
| `-o/--output-dir` | `mixrace_out` | 输出目录 |
| `--repeat-bed` | 无 | 重复区 BED(可选) |
| `-t/--threads` | 12 | 线程数 |
| `-k/--kmer-size` | 21 | k-mer 大小(smudgescope) |
| `-l/--read-length` | 150 | 读长 |
| `--step` | 全跑 | 只跑 1-7 某一步(断点续传友好) |
| `--min-qual` | 30 | 变异位点质量下限 |
| `--min-dp` | 15 | 位点深度下限 |
| `--min-alt-reads` | 3 | ALT 最少支持读条数 |
| `--min-coverage` | 30 | freebayes 最小覆盖 |
| `--min-alt-fraction` | 0.02 | freebayes 保留 2% 低频等位 |
| `--pure-samples` | 无 | 已知纯样品(逗号分隔)，自动按 mean+2SD 校准 het 阈值 |
| `--skip-tree` | 关 | 跳过系统发育树 |
| `--no-checkpoint` / `--dry-run` | 关 | 禁用断点 / 只打印命令 |

## 分析流程 | Pipeline

```text
01 索引+QC        02 比对+markdup       03 freebayes -p 1     04 过滤
bwa-mem2 index    bwa-mem2 mem→        (单倍体,保低频,        QUAL/DP/去repeat
fastp(可跳)       sort-n→fixmate→      AF 用 AO/RO)          保留多等位
                  sort→markdup→index
                                                                    ↓
07 判读+报告 ←── 06 smudgescope ←────────────── 05 AFS 分析
(het率+AFS形态   (k-mer 谱)                (杂合率+形态+优势株占比
 +优势株+深度                                    +VAF 直方图)
 +系统发育树)
```

## 输出 | Output

```text
out_dir/
├── 00_pipeline_info/   software_versions.yml, checkpoints/, index/
├── 01_qc/              clean fastq(raw 输入时)
├── 02_alignment/       {sample}.sorted.markdup.bam
├── 03_variants/        {sample}.raw.vcf.gz(freebayes)
├── 04_filtered/        {sample}.filtered.vcf.gz, filter_stats.tsv
├── 05_vaf/             {sample}.vaf.tsv, {sample}.vaf_histogram.png
├── 06_kmer/{sample}/   smudgescope(02_genomescope/, 03_smudgeplot/)
├── 07_report/          {sample}.report.md
├── 08_tree/            tree.png(ggtree 静态树,叶名=样品[判读]杂合率),
│                       tree.R/tree.ann.tsv(画树脚本与注解), merged.vcf.gz,
│                       vcf2tree/(IQ-TREE2 输出,含 .nwk)
├── alignment_qc/       {sample}.stats.txt(平均深度)
├── summary/            mixrace_report.html(单文件报告), verdict_summary.tsv/.html
└── 99_logs/            mixrace.log
```

## 结果解读 | Interpreting Results

- **看 `summary/mixrace_report.html`**：第一节汇总表一眼看全部样品判读；点开每个样品看依据和三张图。
- **杂合率**是最主要的判据：<0.01% 纯；0.01–0.1% 需排查(repeat 误比对/旁系同源/轻微污染)；0.1–1% 可疑；>1% 不纯。
- **AFS 直方图**：峰贴在两头=纯；中间有峰=有混合；从 0 到 1 糊成一片=多基因型大杂烩。
- **系统发育树**：分支越近基因组越像；叶名格式 `编号[判读]杂合率`(如 `Pb5 [纯] 0.0005%`)。注意混合样品可能表现为**长枝**——因为它同时带多个基因型的信号。
- **判读是参考不是定论**：阈值默认值基于经验，建议用 `--pure-samples` 提供已知纯样品校准，或在汇总表里读纯样品的背景基线自己定界。

## 参数选择建议 | Parameter Guidance

| 场景 | 建议 |
|---|---|
| 常规分析 | 全部默认 |
| 手上有已知纯样品 | `--pure-samples P1,P2` 自动校准 het 阈值 |
| 怀疑稀有基因型被漏 | `--min-alt-fraction 0.02`(默认已保留 2%)可再降 |
| 样品数 <4 或不要树 | `--skip-tree` |
| 中途失败重跑 | 直接重跑原命令(断点续传自动跳过已完成步骤；树注解变化会自动重画) |

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input, -i` | — |  | 原始FASTQ目录(与--clean-fastq-dir二选一)｜Raw FASTQ dir (or --clean-fastq-dir) |
| `--clean-fastq-dir` | — |  | 已清洗FASTQ目录(给则跳过QC)｜Clean FASTQ dir (skip QC) |
| `--genome, -g` | 必填 |  | 参考基因组FASTA｜Reference genome FASTA |
| `--output-dir, -o` | `mixrace_out` |  | 输出目录｜Output directory |
| `--repeat-bed` | — |  | 重复/低复杂度区域BED(可选)｜Repeat/low-complexity BED (optional) |
| `--threads, -t` | `12` | int | 线程数｜Threads |
| `--kmer-size, -k` | `21` | int | K-mer大小｜K-mer size |
| `--read-length, -l` | `150` | int | 测序读长｜Read length |
| `--step` | — | int | 只跑指定步骤1-7(默认全跑)｜Run single step 1-7 (default all) |
| `--no-checkpoint` | `False` |  | 禁用断点续传｜Disable checkpoint |
| `--dry-run` | `False` |  | 只打印命令不执行｜Print commands only |
| `--min-qual` | `30` | int | 变异QUAL下限｜Min QUAL |
| `--min-dp` | `15` | int | 位点深度下限｜Min DP |
| `--min-alt-reads` | `3` | int | ALT支持reads下限｜Min ALT reads |
| `--min-coverage` | `30` | int | freebayes --min-coverage(默认30)｜freebayes min-coverage (default 30) |
| `--min-alt-fraction` | `0.02` | float | freebayes --min-alternate-fraction(默认0.02)｜freebayes min-alternate-fraction (default 0.02) |
| `--pure-samples` | — |  | 已知纯样品(逗号分隔,校准het阈值)｜Known-pure samples (comma-sep, calibrate) |
| `--skip-tree` | `False` |  | 跳过系统发育树｜Skip phylogenetic tree |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | — |  |  |
| `--clean-fastq-dir` | — |  |  |
| `-g, --genome` | 必填 |  |  |
| `-o, --output-dir` | `mixrace_out` |  |  |
| `--repeat-bed` | — |  |  |
| `-t, --threads` | `12` | int |  |
| `-k, --kmer-size` | `21` | int |  |
| `-l, --read-length` | `150` | int |  |
| `--step` | — | int |  |
| `--no-checkpoint` | — | store_false |  |
| `--dry-run` | — | store_true |  |
| `--min-qual` | `30` | int |  |
| `--min-dp` | `15` | int |  |
| `--min-alt-reads` | `3` | int |  |
| `--min-coverage` | `30` | int | freebayes --min-coverage(默认30) |
| `--min-alt-fraction` | `0.02` | float | freebayes --min-alternate-fraction(默认0.02,保低频等位) |
| `--pure-samples` | — |  | 已知纯样品(逗号分隔,校准het阈值)｜known-pure samples (calibrate) |
| `--skip-tree` | — | store_true | 跳过系统发育树｜skip phylogenetic tree |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- bwa-mem2(`cphasing`)、samtools/bcftools/bedtools/freebayes(`align`)、R+ggplot2(`WGCNA_v.1.73`)、**R+ggtree(`r`，画树)**、smudgescope 自带 envs

## 常见问题 | FAQ

- **Q：为什么所有样品都判「混合」？** 检查 mean_depth——深度 <50x 时低频信号不可靠；或参考基因组本身来自混合孢子(群体多样性)会抬高背景杂合率。
- **Q：树上的样品名和表格对不上？** 树叶名是 `编号[判读]杂合率`，与汇总表同源；若刚改过阈值重跑，树会自动重画(注解比对守卫)。
- **Q：为什么没有树？** 样品 <4、filtered VCF 缺失、`--skip-tree`、或 IQ-TREE/vcf2tree 失败(见日志)。
- **Q：报告能发给别人看吗？** `mixrace_report.html` 是单文件，所有图内嵌，直接发送即可。
