# mixrace — WGS 混合小种检测(三分支判读)

一句话:对一群**单倍体病原菌样品**(如油菜根肿菌)的重测序数据做"纯不纯"体检——谁是可以直接保存的纯菌、谁是还能用的优势菌株、谁是必须重新分离纯化的混杂菌株,并给出证据和实验建议。

## 功能概述 | Overview

- **GTX 联合 calling 后端**:fastp QC →(可选寄主剔除)→ `fastq2vcf-gtx` 一步比对+联合分型,全群体共用一张 VCF
- **四层杂合评估**(导师 v4 方法学):L1 AD/DP 排测序错误 → L2 shared/private ALT + **混合伴侣分析** → L3 100kb 窗口分布 → L4 共享热点排除
- **三分支判读 + 实验建议**:纯菌(可保存)/ 优势菌株·参考差异型(可保存,可强制纯合化)/ 混杂菌株(需再分离纯化,附成分推断)
- **reads 账本**:寄主占比、病原 mapping 率、**污染 reads**(两步都没比对上)全程可追溯
- **k-mer 分析只用比对上的 reads**:genomescope 基因组估计 + smudgeplot 不被污染 reads 干扰
- **群体结构 + 全套图**:SNP 距离矩阵、PCA、NJ 树、杂合热图、Manhattan、三面板评估等 9 图

## 快速开始 | Quick Start

示例|Examples: biopytools mixrace -i fastq_dir -g ref.fa --host-genome host.fa -o out_dir/

## 零基础概念速览 | Concepts in plain words

| 术语<br>Term | 通俗解释<br>In plain words |
|---|---|
| 杂合率 (het_rate) | 根肿菌是单倍体,每个位置只该有一种碱基;出现两种(0/1)就像"一个人的体检报告里出现了两个人的指标"——杂合率就是混入"另一份基因组的指纹"的比例 |
| 稳健杂合率 (robust) | 只数"两种碱基都有相当数量reads支持"(altAD≥5 且 altfrac≥0.2)的杂合位点,把低深度的测序假信号剔除后的杂合率 |
| shared/private ALT | 你的杂合位点上的变异碱基,别的样品里也能看到(shared)还是只有你有(private)。private 的多半是错误或无关污染,不算数 |
| 混合伴侣分析 | 在群体里找"与你共享变异的另一份基因组":我的杂合位点上,某个样品大量携带 ALT 甚至纯合(1/1),说明我 ≈ 它 + 参考型的混合。像 Pb9 ≈ 88% Pb22 型 + 12% 参考型 |
| 伴侣 1/1 占比 | 区分两种情况的关键:伴侣在我的杂合位点是**纯合**(我是"它+参考"的混合,即 Pb9-Pb22 模式)还是**杂合**(它自己也是同群混合) |
| 热点区域 | 所有混杂样品一致偏高的窗口——是远缘菌株与参考基因组的系统性差异区,不是样品特异的混杂信号,评估时可剔除 |
| DP 检验 (dp_ratio) | 杂合位点的深度 ÷ 纯合位点深度。>1.5 说明杂合位点反而测得更深——真混合,不是低覆盖错误 |
| 污染 reads | 既没比对上寄主基因组、也没比对上病原基因组的 reads——可能来自其他微生物或低质量 |
| 强制纯合化 | 对判"可保存但有杂合"的样品,取每个位点占多数的碱基当纯合用,适合高精度下游 |

## 输入 | Input

| 输入<br>Input | 要求<br>Requirement |
|---|---|
| `-i` 原始 fastq 目录 | 配对 R1/R2,支持 `_1/_2`、`_R1/_R2`、`.clean.` 等命名;跑 fastp QC |
| `--clean-fastq-dir` | 已清洗 fastq(与 `-i` 二选一,跳过 QC) |
| `-g` 病原参考基因组 | FASTA;GTX 比对与联合 calling 的参考 |
| `--host-genome` 寄主基因组(可选) | 给则先剔除寄主 reads 再进主流程 |
| 样品数 | 建议 ≥4(伴侣分析与群体结构依赖多样本;单样品仍可跑但只有 L1) |

## 参数说明 | Parameters

**寄主剔除与 reads 口径**

**通俗理解|In plain words:** 寄主剔除管"把混进来的植物 DNA 扔掉再分析";`--min-mapq` 是"什么算真的比对上"的统一门槛,一般不用动。

| 参数<br>Parameter | 默认<br>Default | 说明<br>Description |
|---|---|---|
| `--host-genome` | 无(关) | 寄主基因组 FASTA;给则整对剔除寄主 reads 并报告寄主占比 |
| `--min-mapq` | 20 | 比对质量阈值:寄主判定 + mapped reads 提取 + 统计口径;0=不过滤 |

**三分支判读阈值**

**通俗理解|In plain words:** 这三个数决定"多杂算杂"。导师口径:杂合率 0.1% 分纯/不纯;伴侣同时满足"携带 ALT≥80%"和"其中纯合 1/1≥50%"才判混杂(Pb9-Pb22 模式)。没把握就不动,真实数据校验后可调。

| 参数<br>Parameter | 默认<br>Default | 说明<br>Description |
|---|---|---|
| `--pure-het-threshold` | 0.001 | 总杂合率低于此值(0.1%)判纯菌 |
| `--partner-alt-rate` | 0.8 | 混合伴侣:我的杂合位点上伴侣携带 ALT 的比例阈值 |
| `--partner-hom-rate` | 0.5 | 混合伴侣:伴侣携带中纯合 1/1 的占比阈值 |
| `--min-sites` | 1000 | 有 GT 位点低于此数判 uncertain(数据不足) |

**热点识别(L4)**

**通俗理解|In plain words:** 自动找"所有混杂样品一起偏高"的基因组窗口并剔除,让杂合率不被参考差异区虚抬。窗口大小、几倍算偏高都有默认,一般不用动。

| 参数<br>Parameter | 默认<br>Default | 说明<br>Description |
|---|---|---|
| `--window-size` | 100000 | 窗口大小 bp(100kb) |
| `--hotspot-fold` | 2.0 | 窗口杂合率 > 该倍数×自身全基因组率算偏高 |
| `--hotspot-min-median` | 0.10 | 窗口在候选样品中的中位杂合率下限 |
| `--repeat-bed` | 无 | 手动追加排除区域 BED(与自动热点取并集) |

**执行**

| 参数<br>Parameter | 默认<br>Default | 说明<br>Description |
|---|---|---|
| `-t/--threads` | 12 | 线程数(samtools/bwa/GTX/smudgescope 全部生效) |
| `--sample-parallel` | 1 | 样本级并行数:寄主剔除/mapped提取/reads统计同时跑 N 个样本,每 worker 线程=threads/N。25 样本推荐 4-8 |
| `--step` | 全跑 | 1=QC+寄主剔除 2=GTX 3=评估判读 4=k-mer 5=图+报告 |
| `--no-checkpoint` / `--dry-run` | 关 | 禁用断点 / 只打印命令 |
| `-k/--kmer-size`、`-l/--read-length` | 21 / 150 | smudgescope 参数 |

## 分析流程 | Pipeline

```
raw fastq ─fastp→ clean ─[寄主剔除]→ nohost fastq
   │                                      │
   │              ①寄主占比(host_rate)    ▼
   │                          ②GTX 比对+联合calling → joint VCF + per-sample BAM
   │                                      │
   │                     ③污染reads = 总-寄主-mapped
   │                                      ▼
   │   bcftools query GT/AD/DP 长表 → L1排错 → L2伴侣 → L3窗口 → L4热点
   │                                      │
   │                        群体: 距离矩阵+PCA+NJ → 三分支判读+实验建议
   │                                      ▼
   └────────────────→ ④mapped reads → smudgescope(genomescope+smudgeplot)
                                          ▼
                          ⑤全套图 + 判读汇总 + 单样品报告 + HTML
```

## 输出 | Output

```
out/                          # by-step:所有样本共享编号步骤目录,文件名用 {sample}_ 前缀区分
├── 00_pipeline_info/   index_host/ checkpoints/ software_versions.yml(含全部阈值)
├── 01_qc/              fastp 产物
├── 02_host_filter/     (仅给 --host-genome 时存在){sample}_1/2.nohost.fq.gz · host_filter.tsv · host_stats.tsv(reads账本)
├── 03_gtx/             03_mapping/bam/{sample}.bam · 04_joint_calling/gtx_joint_raw.vcf.gz
├── 04_het_eval/        gt_ad_dp.tsv · L1_杂合统计 · L2_shared_private · L2_shared_only评估
│                       混合伴侣矩阵 · 混合伴侣top · L3_窗口杂合率 · L4_共享热点窗口
│                       hotspots.bed · L4_排除热点前后对比 · 距离矩阵 · PCA坐标 · nj_tree.nwk
│                       verdict_table.tsv(判读+证据链+建议) ·
│                       alignment_qc/{sample}.stats.txt(深度缓存)
├── 05_kmer/            mapped_fastq/{sample}_1/2.mapped.fq.gz + smudgescope 输出
├── 06_figures/         9 张图(热图/Manhattan/距离/PCA/NJ/altfrac/三面板等)
├── 07_report/          {sample}.report.md(证据链)
├── summary/            verdict_summary.tsv · mixrace_report.html(自包含)
├── tmp/                临时文件(运行中,结束清理)
└── 99_logs/
```

## 结果解读 | Interpreting Results

**看 `summary/verdict_summary.tsv`(或 HTML)的 verdict + advice 两列:**

| 判读<br>Verdict | 含义<br>Meaning | 建议(给做实验的同学)<br>Advice | 关键证据<br>Key evidence |
|---|---|---|---|
| **pure 纯菌** | 杂合率 <0.1% | **可保存** | het_rate、robust_rate 都极低 |
| **divergent 优势菌株/参考差异型** | 有杂合但无强伴侣 | **可保存**;高精度下游可强制纯合化 | shared% 高、伴侣多为 0/1、排除热点后仍稳;`轻度`=robust<0.1% |
| **contaminated 混杂菌株** | 强混合伴侣(Pb9-Pb22 模式) | **需再分离纯化** | top_partner 两率超阈值;mix_proportion 给成分(如 88% B型+12% 参考型) |
| uncertain | 位点不足 | 补数据 | n_sites < min_sites |

**辅助判据:** `dp_ratio`>1.5(杂合位点反而更深=真混合);`排除热点后杂合率`仍高=证据坚实;`host_rate/污染reads占比`高=样品制备问题。

## 参数选择建议 | Parameter Guidance

- **默认阈值即导师 v4 口径**,先用默认跑;判读边界样品(如恰在 0.1% 附近)再结合证据链人工复核
- 样品只有 2-3 个:伴侣分析统计力弱,判读以 L1/热点为主(报告会显示证据不足)
- 寄主污染重的样品务必给 `--host-genome`——寄主 reads 会同时虚抬"污染reads"并干扰 k-mer 基因组估计
- 想看"不剔除热点"的原始杂合率:对比 `l4_hotspot_excluded_compare.tsv` 两列即可,无需重跑

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
| `--repeat-bed` | — |  | 额外排除区域BED(与自动热点并集)｜Extra exclude BED (merged with hotspots) |
| `--host-genome` | — |  | 寄主基因组FASTA(给则比对寄主并整对剔除寄主reads,报告寄主占比)｜Host genome FASTA (deplete host reads, report host rate) |
| `--min-mapq` | `20` | int | 比对质量阈值:mapped reads提取+统计口径(0=不过滤)｜Min MAPQ: mapped-read extraction + stats (0=off) |
| `--threads, -t` | `12` | int | 线程数｜Threads |
| `--sample-parallel` | `1` | int | 样本级并行数(每worker线程=threads/N)｜Per-sample parallelism |
| `--kmer-size, -k` | `21` | int | K-mer大小｜K-mer size |
| `--read-length, -l` | `150` | int | 测序读长｜Read length |
| `--step` | — | int | 只跑指定步骤1-5(1=QC+寄主剔除 2=GTX 3=评估判读 4=k-mer 5=图+报告)｜Run single step 1-5 (default all) |
| `--no-checkpoint` | `False` |  | 禁用断点续传｜Disable checkpoint |
| `--dry-run` | `False` |  | 只打印命令不执行｜Print commands only |
| `--pure-het-threshold` | `0.001` | float | 总杂合率低于此值判纯菌(0.001=0.1%)｜Pure threshold |
| `--partner-alt-rate` | `0.8` | float | 混合伴侣:ALT携带率阈值｜Partner ALT-carrier threshold |
| `--partner-hom-rate` | `0.5` | float | 混合伴侣:伴侣纯合1/1占比阈值｜Partner homozygous threshold |
| `--min-sites` | `1000` | int | 最低有GT位点数,低于判uncertain｜Min called sites |
| `--window-size` | `100000` | int | 热点窗口大小bp｜Hotspot window size |
| `--hotspot-fold` | `2.0` | float | 热点:窗口杂合率>该倍数x自身全基因组率｜Hotspot fold |
| `--hotspot-min-median` | `0.1` | float | 热点:窗口候选中位杂合率下限｜Hotspot min median rate |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | — |  |  |
| `--clean-fastq-dir` | — |  |  |
| `-g, --genome` | 必填 |  |  |
| `-o, --output-dir` | `mixrace_out` |  |  |
| `--repeat-bed` | — |  | 额外排除区域BED(与自动热点并集)｜extra exclude BED (merged with hotspots) |
| `--host-genome` | — |  | 寄主基因组FASTA(给则比对寄主并整对剔除寄主reads,报告寄主占比)｜host genome FASTA (deplete host reads, report host rate) |
| `--min-mapq` | `20` | int | 比对质量阈值:mapped reads提取+统计口径(0=不过滤)｜min MAPQ for mapped-read extraction + stats (0=off) |
| `-t, --threads` | `12` | int |  |
| `--sample-parallel` | `1` | int | 样本级并行数(寄主剔除/mapped提取/reads统计同时跑N个样本,每worker线程=threads/N;默认1串行)｜per-sample parallelism |
| `-k, --kmer-size` | `21` | int |  |
| `-l, --read-length` | `150` | int |  |
| `--step` | — | int |  |
| `--no-checkpoint` | — | store_false |  |
| `--dry-run` | — | store_true |  |
| `--pure-het-threshold` | `0.001` | float | 总杂合率低于此值判纯菌(默认0.001=0.1%%)｜pure threshold |
| `--partner-alt-rate` | `0.8` | float | 混合伴侣:ALT携带率阈值(默认0.8)｜partner ALT-carrier threshold |
| `--partner-hom-rate` | `0.5` | float | 混合伴侣:伴侣纯合1/1占比阈值(默认0.5)｜partner homozygous threshold |
| `--min-sites` | `1000` | int | 最低有GT位点数,低于判uncertain(默认1000)｜min called sites |
| `--window-size` | `100000` | int | 热点窗口大小bp(默认100kb)｜hotspot window size |
| `--hotspot-fold` | `2.0` | float | 热点:窗口杂合率>该倍数×自身全基因组率(默认2)｜hotspot fold |
| `--hotspot-min-median` | `0.1` | float | 热点:窗口在候选中的中位杂合率下限(默认0.1)｜hotspot min median rate |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

| 工具/库<br>Tool/Lib | 用途<br>Purpose | 环境<br>Env |
|---|---|---|
| fastq2vcf-gtx(biopytools 模块) | 比对+gVCF+联合 calling | `~/software/gtx`(第三方二进制,直调) |
| bwa-mem2 / samtools / bcftools | 寄主索引 / 提取计数 / query | cphasing / align 域环境 |
| smudgescope(biopytools 模块) | genomescope + smudgeplot | smudgescope 环境 |
| numpy / matplotlib / biopython | 评估引擎 / 全套图 / NJ 树 | biopytools env |

## 常见问题 | FAQ

- **Pb9 那种样品会怎么判?** contaminated,partner=Pb22,报告写明"≈88% Pb22型+12% 参考型",建议再分离纯化
- **群2/群3 那种 4-5% 杂合的近缘样品呢?** divergent(伴侣互为 0/1,纯合占比不过线)→ 可保存,需要高精度时强制纯合化
- **旧 v0.2 输出目录能接着跑吗?** 不能,目录结构与后端都变了(02_alignment/03_variants 已不存在),换新目录重跑
- **`--step 3` 要重跑但 GTX 太慢?** GTX 有断点;VCF 已在 `03_gtx/` 就直接 `--step 3`,秒级起评估
- **图里中文变方框?** 系统无中文字体时模块自动退英文标签,不影响数据
