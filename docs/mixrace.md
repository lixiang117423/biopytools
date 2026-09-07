# mixrace — WGS 混合小种检测(三分支判读)

一句话:对一群**单倍体病原菌样品**(如油菜根肿菌)的重测序数据做"纯不纯"体检——谁是可以直接保存的纯菌、谁是还能用的优势菌株、谁是必须重新分离纯化的混杂菌株,并给出证据和实验建议。

## 功能概述 | Overview

- **GTX 联合 calling 后端**:fastp QC →(可选寄主剔除)→ `fastq2vcf-gtx` 一步比对+联合分型,全群体共用一张 VCF
- **四层杂合评估**(导师 v4 方法学):L1 AD/DP 排测序错误 → L2 shared/private ALT + **混合伴侣分析** → L3 100kb 窗口分布 → L4 共享热点排除
- **三分支判读 + 实验建议**:纯菌(可保存)/ 优势菌株·参考差异型(可保存,可强制纯合化)/ 混杂菌株(需再分离纯化,附成分推断)
- **reads 账本**:寄主占比、病原 mapping 率、**污染 reads**(两步都没比对上)全程可追溯
- **kraken2+bracken 污染评估**(step 6):clean reads 逐样本做物种分类,回答"样品里除了病原和寄主,还混了哪些微生物、各占多少"
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
| 分类 (classified/unclassified) | kraken2 把每条 read 贴标签:"像哪个物种"。贴不上的叫 unclassified——目标菌根肿菌不在公共库里,它的 reads 基本都落在这里,属正常 |
| 物种丰度 (bracken) | 光数 reads 会把"长得像的近缘种"数错;bracken 按各物种基因组特征重新分摊,给出更接近真实的占比 |
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

**AF 直接判定杂合度(论文式对照,默认开)**

**通俗理解|In plain words:** 这是一把"备用尺子"。默认的总杂合率要靠软件先给每个位点判基因型(GT),而 GATK 这类软件的先验倾向"变异碱基占一半才算杂合"——如果两个菌株按 7:3 混合,变异只占 30%,软件会误判成纯合,杂合率被低估。论文式做法不看软件判的基因型,直接数每个变异位点上"参考碱基几条 read、变异碱基几条 read":变异占比在 5%-95% 之间且变异 reads ≥3 条、深度 ≥10,就算杂合位点。**默认随全流程自动跑**(`--no-af-het-eval` 关闭),verdict_table 和汇总报告会多一列"AF口径杂合率"与总杂合率**并排对比**(判读本身不变,仍用 GT 口径)。注意两者**分母不同**:AF 口径 = 杂合/(杂合+纯合变异),不含纯合参考位点,数值天然偏高,只用于对照"GT 先验有没有压低非对称混合的信号",不能和总杂合率直接比大小。三个阈值一般不用动。

| 参数<br>Parameter | 默认<br>Default | 说明<br>Description |
|---|---|---|
| `--no-af-het-eval` | 开(默认跑) | 关闭论文式 AF 直接判定(Cao et al. 2026);关闭后不再产出 het_rate_af 列 |
| `--af-het-min-frac` | 0.05 | 杂合判定的 alt 比例下限(杂合区间 [0.05, 0.95] 闭区间;须在 0-0.5 之间) |
| `--af-het-min-depth` | 10 | 位点参与判定的最低深度 |
| `--af-het-min-alt-ad` | 3 | 杂合判定的最低变异 reads 数(压测序错误) |

**污染评估(kraken2+bracken,step 6)**

**通俗理解|In plain words:** 管的是"给每条 read 查户口——你到底是谁家的"。结果回答:样品里除了目标菌,还混了哪些微生物、各占百分之几。**数据库默认用超算上的 `~/database/kraken2`(238GB 大库),默认模式要把整个库装进内存——整个作业必须提交到 ≥250GB 内存的节点**;内存不够就加 `--kraken-memory-mapping`(慢一些但省内存)。不想要这一步就 `--skip-kraken2`,其余参数一般不用动。

| 参数<br>Parameter | 默认<br>Default | 说明<br>Description |
|---|---|---|
| `--skip-kraken2` | 不跳过 | 关闭污染评估(默认全流程自动跑) |
| `--kraken2-db` | `~/database/kraken2` | kraken2/bracken 数据库 |
| `--kraken-memory-mapping` | 关 | 省内存模式(适合内存不足节点,速度变慢) |
| `--bracken-level` | S | 丰度统计层级:S=种 / G=属 / F=科(粗排查用 G 更稳) |

**执行**

| 参数<br>Parameter | 默认<br>Default | 说明<br>Description |
|---|---|---|
| `-t/--threads` | 12 | 线程数(samtools/bwa/GTX/smudgescope/kraken2 全部生效) |
| `--sample-parallel` | 1 | 样本级并行数:寄主剔除/mapped提取/reads统计同时跑 N 个样本,每 worker 线程=threads/N。25 样本推荐 4-8。注意 kraken2 因内存按样本串行,不受此参数影响 |
| `--step` | 全跑 | 1=QC+寄主剔除 2=GTX 3=评估判读 4=k-mer 5=图+报告 6=kraken2污染评估 |
| `--no-checkpoint` / `--dry-run` | 关 | 禁用断点 / 只打印命令 |
| `-k/--kmer-size`、`-l/--read-length` | 21 / 150 | smudgescope 参数;`-l` 同时是 bracken 读长(自动吸附到库内可用档位 50/75/100/150/200/250/300) |

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

01_qc clean ──→ ⑥kraken2 逐样本分类(寄主剔除前,寄主占比也可见)
                      └→ bracken 物种丰度 → contamination_summary/detail
```

## 输出 | Output

```
out/                          # by-step:所有样本共享编号步骤目录,文件名用 {sample}_ 前缀区分
├── 00_pipeline_info/   index_host/ checkpoints/ software_versions.yml(含全部阈值)
├── 01_qc/              fastp 产物
├── 02_host_filter/     (仅给 --host-genome 时存在){sample}_1/2_nohost.fq.gz · host_filter.tsv · host_stats.tsv(reads账本)
├── 03_gtx/             03_mapping/bam/{sample}.bam · 04_joint_calling/gtx_joint_raw.vcf.gz
├── 04_het_eval/        gt_ad_dp.tsv · L1_杂合统计 · L2_shared_private · L2_shared_only评估
│                       混合伴侣矩阵 · 混合伴侣top · L3_窗口杂合率 · L4_共享热点窗口
│                       hotspots.bed · L4_排除热点前后对比 · 距离矩阵 · PCA坐标 · nj_tree.nwk
│                       verdict_table.tsv(判读+证据链+建议;默认多一列 het_rate_af) ·
│                       af_ad_dp.tsv · l1_het_af_based.tsv(论文式 AF 直接判定,不依赖GT;--no-af-het-eval 时不产出) ·
│                       alignment_qc/{sample}_stats.txt(深度缓存)
├── 05_kmer/            mapped_fastq/{sample}_1/2_mapped.fq.gz + smudgescope 输出
├── 06_figures/         9 张图(热图/Manhattan/距离/PCA/NJ/altfrac/三面板等)
├── 07_report/          {sample}_report.md(证据链)
├── 08_contamination/   kraken2/{sample}_k2_report.txt · bracken/{sample}_bracken_S.txt ·
│                       contamination_summary.tsv(每样本一行) · contamination_detail.tsv(样本×物种长表)
├── summary/            verdict_summary.tsv · verdict_summary.xlsx(中文/英文双sheet) · mixrace_report.html(自包含交互报告) · verdict_summary.html(汇总表独立页)
├── tmp/                临时文件(运行中,结束清理)
└── 99_logs/
```

### mixrace_report.html 报告结构 | HTML report layout

> **通俗理解|In plain words:** 双击打开就能看,不需要联网;像一份体检报告——先给结论总览,再给全景图,最后是每个样品的完整证据。

- **页头总览**:判读构成条 + 四张计数卡(纯菌/优势菌株/混杂/不确定各多少份),一眼看清整批判读构成
- **01 判读汇总表**:点击表头排序(数值列按原始数值排序,不受显示格式影响);横向滚动时表头和样品列固定不动;建议列超长省略号显示,鼠标悬浮看全文
- **02 图版**:全部 9 张图按阅读顺序排列(全景→定位→关系→机理),每张带中文标题和一句话图注
- **03 逐样品证据链**:默认折叠;顶部搜索框按样品名过滤(汇总表同步过滤);「全部展开」后可直接打印(打印时自动展开全部,搜索框和按钮自动隐藏)
- 颜色语义与图一致:绿=纯菌、橙=优势菌株/参考差异型、红=混杂、灰=不确定;所有颜色都伴随文字标签,不单独依赖颜色区分

## 结果解读 | Interpreting Results

**看 `summary/verdict_summary.tsv`(或 HTML)的 verdict + advice 两列:**

| 判读<br>Verdict | 含义<br>Meaning | 建议(给做实验的同学)<br>Advice | 关键证据<br>Key evidence |
|---|---|---|---|
| **pure 纯菌** | 杂合率 <0.1% | **可保存** | het_rate、robust_rate 都极低 |
| **divergent 优势菌株/参考差异型** | 有杂合但无强伴侣 | **可保存**;高精度下游可强制纯合化 | shared% 高、伴侣多为 0/1、排除热点后仍稳;`轻度`=robust<0.1% |
| **contaminated 混杂菌株** | 强混合伴侣(Pb9-Pb22 模式) | **需再分离纯化** | top_partner 两率超阈值;mix_proportion 给成分(如 88% B型+12% 参考型) |
| uncertain | 位点不足 | 补数据 | n_sites < min_sites |

**辅助判据:** `dp_ratio`>1.5(杂合位点反而更深=真混合);`排除热点后杂合率`仍高=证据坚实;`host_rate/污染reads占比`高=样品制备问题。

### 两把杂合率尺子怎么对照 | GT-based vs AF-based

> **通俗理解|In plain words:** 一个看"软件判的基因型",一个看"reads 里数出来的碱基比例"。AF 口径列默认就有(`--no-af-het-eval` 才会关)。

- **两数接近** → 该样品的混合接近对称(变异占比≈50%),GT 口径没吃亏
- **AF口径 >> 总杂合率** → 存在大量 5%-45% 的非对称混合信号被 GATK 先验压成"纯合"了——非对称多菌株混合(如 7:3)的典型表现;这类样品优先按"需再分离纯化"处理
- **AF口径也低** → 杂合确实少,GT 口径低估的担忧可以排除
- **AF口径杂合率显示"—"(空)** → 小分母守卫:该样品 杂合+纯合变异 位点总数低于 `--min-sites`(默认1000),论文分母在这些样品上退化为噪声(与参考几乎一致的纯菌只有几十个变异位点,几十个位点的比率没有统计意义)。看 `l1_het_af_based.tsv` 的 n_het_af/n_hom_alt_af 计数即可,这些样品本身就是"与参考几乎一致"的结论
- 两个口径**分母不同**(AF 口径分母不含纯合参考位点),数值本身不可直接比大小,看的是**量级差距**;完整逐样本数字见 `04_het_eval/l1_het_af_based.tsv`(n_sites_eval/n_het_af/n_hom_alt_af 全字段)

### 污染评估怎么看 | Reading contamination results

**看 `08_contamination/contamination_summary.tsv`(每样本一行):**

| 列<br>Column | 通俗解释<br>In plain words |
|---|---|
| classified_pct / unclassified_pct | 成功查到户口的 read 占比 / 查不到的占比。**根肿菌不在公共库里,所以 unclassified 的大头就是目标菌,unclassified 高是好事** |
| top_species / top_species_pct | 占比第一的物种。样品没除寄主时它多半是**油菜/拟南芥(寄主)**;除净寄主后它就是**最大污染源** |
| other_classified_pct | 除第一物种外的已分类占比——这个数大,说明污染是"多种微生物混入"而非单一来源 |
| n_species_ge_1pct | 占比 ≥1% 的物种个数。干净样品通常 ≤2(寄主+一个杂菌),连续两位数就是脏了 |

**判读口径:** 目标菌(根肿菌)和寄主之外的任何具体物种占比 ≥1% 就值得追查;≥5% 基本可确认污染,对照 `contamination_detail.tsv`(样本×物种长表,含 0.1% 以上全部物种)定位具体是谁。注意 top_species_pct 的分母是**全部 reads(含 unclassified)**,与 bracken 原生输出(分母只有已分类 reads)口径不同。

## 参数选择建议 | Parameter Guidance

- **默认阈值即导师 v4 口径**,先用默认跑;判读边界样品(如恰在 0.1% 附近)再结合证据链人工复核
- 样品只有 2-3 个:伴侣分析统计力弱,判读以 L1/热点为主(报告会显示证据不足)
- 寄主污染重的样品务必给 `--host-genome`——寄主 reads 会同时虚抬"污染reads"并干扰 k-mer 基因组估计
- 想看"不剔除热点"的原始杂合率:对比 `l4_hotspot_excluded_compare.tsv` 两列即可,无需重跑
- 旧结果想补 AF 口径对照:直接 `--step 3` 重跑(AF 分支默认开)——GT 侧四层评估有断点会跳过,只跑 AF 分支并把 `het_rate_af` 列整合进已有 `verdict_table.tsv`;想刷新报告再 `--step 5`

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
| `--step` | — | int | 只跑指定步骤1-6(1=QC+寄主剔除 2=GTX 3=评估判读 4=k-mer 5=图+报告 6=kraken2污染评估)｜Run single step 1-6 (default all) |
| `--skip-kraken2` | `False` |  | 跳过kraken2+bracken污染评估｜Skip contamination assessment |
| `--kraken2-db` | `~/database/kraken2` |  | kraken2/bracken数据库(默认模式内存需约DB大小,~240GB)｜kraken2 db (RAM ~ DB size in default mode) |
| `--kraken-memory-mapping` | `False` |  | kraken2省内存模式(慢)｜kraken2 --memory-mapping (slower, less RAM) |
| `--bracken-level` | `S` | D/P/C/O/F/G/S | bracken丰度层级(S=种)｜bracken abundance rank (S=species) |
| `--no-checkpoint` | `False` |  | 禁用断点续传｜Disable checkpoint |
| `--dry-run` | `False` |  | 只打印命令不执行｜Print commands only |
| `--pure-het-threshold` | `0.001` | float | 总杂合率低于此值判纯菌(0.001=0.1%)｜Pure threshold |
| `--partner-alt-rate` | `0.8` | float | 混合伴侣:ALT携带率阈值｜Partner ALT-carrier threshold |
| `--partner-hom-rate` | `0.5` | float | 混合伴侣:伴侣纯合1/1占比阈值｜Partner homozygous threshold |
| `--min-sites` | `1000` | int | 最低有GT位点数,低于判uncertain｜Min called sites |
| `--window-size` | `100000` | int | 热点窗口大小bp｜Hotspot window size |
| `--hotspot-fold` | `2.0` | float | 热点:窗口杂合率>该倍数x自身全基因组率｜Hotspot fold |
| `--hotspot-min-median` | `0.1` | float | 热点:窗口候选中位杂合率下限｜Hotspot min median rate |
| `--no-af-het-eval` | `False` |  | 关闭论文式AF直接判定杂合度(默认开启,verdict_table追加het_rate_af列与GT口径并列对比)｜Disable paper-style AF-based het eval (on by default) |
| `--af-het-min-frac` | `0.05` | float | AF杂合:alt比例下限,杂合区间[min,1-min]闭｜Min alt fraction |
| `--af-het-min-depth` | `10` | int | AF杂合:参与判定的最低深度｜Min depth to evaluate |
| `--af-het-min-alt-ad` | `3` | int | AF杂合:杂合判定最低alt reads数｜Min alt AD |

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
| `--no-af-het-eval` | — | store_false | 关闭论文式AF直接判定杂合度(默认开启;Cao et al. 2026,不依赖GT,从联合VCF取AD/DP,verdict_table 追加 het_rate_af 列与GT口径并列对比)｜disable paper-style AF-based het (on by default) |
| `--af-het-min-frac` | `0.05` | float | AF杂合:alt比例下限,杂合区间[min,1-min]闭(默认0.05)｜min alt fraction |
| `--af-het-min-depth` | `10` | int | AF杂合:参与判定的最低深度(默认10)｜min depth to evaluate |
| `--af-het-min-alt-ad` | `3` | int | AF杂合:杂合判定最低alt reads数(默认3)｜min alt AD |
| `--skip-kraken2` | — | store_false | 跳过 kraken2+bracken 污染评估(默认跑)｜skip contamination assessment |
| `--kraken2-db` | `~/database/kraken2` |  | kraken2/bracken 数据库(默认~/database/kraken2,内存需约DB大小)｜kraken2 db (RAM ~ DB size) |
| `--kraken-memory-mapping` | — | store_true | kraken2 省内存模式(慢,适合内存不足节点)｜kraken2 --memory-mapping (slower) |
| `--bracken-level` | `S` | D/P/C/O/F/G/S | bracken 丰度层级(默认S=种)｜bracken rank (S=species) |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

| 工具/库<br>Tool/Lib | 用途<br>Purpose | 环境<br>Env |
|---|---|---|
| fastq2vcf-gtx(biopytools 模块) | 比对+gVCF+联合 calling | `~/software/gtx`(第三方二进制,直调) |
| bwa-mem2 / samtools / bcftools | 寄主索引 / 提取计数 / query | cphasing / align 域环境 |
| kraken2 2.17 / bracken 3.0 | 物种分类 / 丰度重估 | kraken_v.2.17 |
| kraken2 数据库 | PlusPF 大库(238GB,含真菌/原生/植物) | `~/database/kraken2` |
| smudgescope(biopytools 模块) | genomescope + smudgeplot | smudgescope 环境 |
| numpy / matplotlib / biopython | 评估引擎 / 全套图 / NJ 树 | biopytools env |

## 常见问题 | FAQ

- **Pb9 那种样品会怎么判?** contaminated,partner=Pb22,报告写明"≈88% Pb22型+12% 参考型",建议再分离纯化
- **群2/群3 那种 4-5% 杂合的近缘样品呢?** divergent(伴侣互为 0/1,纯合占比不过线)→ 可保存,需要高精度时强制纯合化
- **AF口径杂合率比总杂合率高很多,是bug吗?** 不是。分母口径不同(AF 口径分母=杂合+纯合变异位点,不含大量纯合参考位点,天然偏高),且 AF 口径把 5%-45% 的非对称混合信号也计为杂合而 GT 口径会漏掉。两列看的是量级差距,不是精确比值;详见「两把杂合率尺子怎么对照」
- **纯菌样品的 AF口径杂合率显示"—"?** 小分母守卫:纯菌与参考几乎一致,杂合+纯合变异位点常常只有几十个,论文比率在这点分母上是纯噪声(实测 25 个位点里 24 个杂合会显示 96%)。分母不足 `--min-sites`(默认1000)时该列置空、计数照常保留,防误导
- **旧 v0.2 输出目录能接着跑吗?** 不能,目录结构与后端都变了(02_alignment/03_variants 已不存在),换新目录重跑
- **`--step 3` 要重跑但 GTX 太慢?** GTX 有断点;VCF 已在 `03_gtx/` 就直接 `--step 3`,秒级起评估
- **kraken2 报内存不足/OOM 被杀?** 默认模式要把 238GB 数据库整装进内存,作业须提交到 ≥250GB 内存节点;上不了大内存节点就加 `--kraken-memory-mapping`(省内存但明显变慢)
- **unclassified 占比很高正常吗?** 正常。目标菌根肿菌**不在公共数据库里**,它的 reads 全落 unclassified;反而要警惕 classified 高的样品——说明混了大量库内微生物
- **top_species 是油菜?** step 6 吃的是寄主剔除**前**的 clean reads,寄主占比本来就该看见;想看剔除后的组成,拿 `02_host_filter/` 的 nohost fastq 自行跑 `--step 6 --clean-fastq-dir`
- **读长 120 会被怎么处理?** bracken 只认库内档位(50/75/100/150/200/250/300),模块自动吸附到最近的 100 并打 WARNING,无需手动换算
- **图里中文变方框?** 系统无中文字体时模块自动退英文标签,不影响数据
