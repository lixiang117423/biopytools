# annorefine - 基因组注释精修(端到端) | Genome annotation refinement (end-to-end)

一句话理解：**给一条命令，自动完成「基因结构预测(BRAKER) + 用近缘蛋白查漏补缺 + 拆开错误合并的基因」，最后给你一个尽量完整的基因注释文件(GFF3)**。
它解决基因组注释里最常见的两类遗漏——基因没被预测出来(漏检)、以及几个基因被错误合并成一个(折叠)——并用「蛋白相似性 + 是否真的能翻译 + 是否有转录表达」三重证据把假基因挡在门外。

## 功能概述 | Overview { #overview }

- **端到端**：一个命令跑完 BRAKER 注释和同源查漏补缺两个阶段，直接输出整合后的 GFF3
- **同源查漏补缺**：用近缘物种蛋白(miniprot 比对)找出 BRAKER 没预测到的新基因
- **合并拆分**：识别并拆开 BRAKER 把多个基因错误折叠成的一个基因
- **三重质控**：完整 ORF 检查(ATG 开头 + 终止密码子 + 长度 3 倍数) + 蛋白相似度/覆盖度 + RNA-seq 表达证据(唯一比对 reads 的深度与覆盖广度)
- **小蛋白回收通道**：可选，找回被最小 CDS 长度阈值整类丢掉的小蛋白(如效应子)
- **断点续传**：BRAKER 和证据扫描都按产物存在性跳过，中断后重跑不重复计算

## 快速开始 | Quick Start { #quick-start }

```bash
biopytools annorefine -g genome.fa -s psojae -p prot.fa -o out/
```

最小输入：一个未屏蔽(mask)的原始基因组 FASTA + 一个物种名 + 一个近缘蛋白文件(或目录)。有 RNA-seq 时强烈建议加上 `--rnaseq-dirs`，表达证据会让查漏补缺准很多。

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语 | 通俗理解 |
|------|----------|
| 基因注释 | 在基因组这条「长字符串」上标出哪里是基因、哪里是外显子、能翻译成什么蛋白 |
| BRAKER | 一款主流的基因结构预测软件，靠「蛋白相似 + 转录本比对」自动学出基因模型 |
| 屏蔽(mask) | 把基因组里高度重复的区域(如转座子)盖住，避免这些「噪音」干扰基因预测 |
| 重复序列/TE | 基因组里会自我复制、搬家的序列，常被误当成基因 |
| 查漏补缺 | 预测软件总会漏掉一些基因，用近缘物种的已知蛋白反查一遍把漏的补上 |
| miniprot | 把蛋白序列「回贴」到基因组上的工具，用来定位蛋白证据落在哪里 |
| 折叠基因 | 两个相邻基因被预测软件错当成一个，需要拆开 |
| 效应子 | 病原菌(如疫霉)分泌出来攻击寄主的小蛋白，常藏在重复区、表达量低、基因短 |
| GFF3 | 基因组注释的标准文件格式，记录每个基因/外显子的坐标和结构 |
| 断点续传 | 中断后重跑，已经完成的步骤自动跳过，不重复算 |

## 输入 | Input { #input }

### 基因组(-g)

**必须未屏蔽(mask)**：BRAKER 阶段内部会自己做屏蔽，查漏补缺阶段要用原始序列(屏蔽会把真基因也盖掉)。FASTA 格式。

### 近缘蛋白(-p)

一个 FASTA 文件或一个目录。给目录时程序自动识别并合并目录里的蛋白文件；给文件时程序会自动清理序列里非标准氨基酸字符(如 `.`、数字)。越近缘、越全，查漏补缺越准。

### RNA-seq(可选，强烈建议)

`--rnaseq-dirs` 逗号分隔的目录列表，目录里是成对的二代测序数据，默认按 `_1.clean.fq.gz` / `_2.clean.fq.gz` 识别双端文件。

### 三代转录本(可选)

`--isoseq` 三代全长转录本，文件或目录。

## 参数说明 | Parameters { #parameters }

### 必需参数 | Required

**通俗理解|In plain words:** 这四个是启动的最小集合：基因组、物种名、近缘蛋白、输出目录。物种名只用于 BRAKER 输出文件的命名，不参与生物学判断；输出目录下会按步骤编号生成一系列子目录。

### BRAKER 证据与通用参数 | BRAKER evidence & general

**通俗理解|In plain words:** 这里控制「用哪些证据帮 BRAKER 找基因」。蛋白证据(`-p`)必填；`--rnaseq-dirs` 和 `--isoseq` 可选但强烈建议至少给一个转录组证据。`--fungus` 是真菌/卵菌模式，疫霉这类物种保持默认开(即不要用 `--no-fungus`)。`--threads` 只影响速度不影响结果，按机器核数给。`--singularity-image` / `--no-singularity` 是 BRAKER3 的容器运行方式，**一般不用动**。

### 重复序列屏蔽 | Repeat masking

**通俗理解|In plain words:** 决定要不要先屏蔽重复区再做预测。默认**做屏蔽**(对多数物种更稳)；`--skip-repeat` 跳过屏蔽。`--skip-repeat-filter` 跳过「重复库过滤」(默认过滤)，`--skip-rescue` 控制「证据还原」(默认不还原，即 rescue 关闭)。**一般保持默认即可**；只在重复区特别多、或明确知道要保住重复区的基因时才调。

### 查漏补缺判据 | Gap-filling criteria

**通俗理解|In plain words:** 控制「什么样的蛋白命中才算一个新基因、要不要拆合并基因」。`--split-min-copy-coverage` 是拆合并基因的判据——一个命中要覆盖蛋白多大比例才算「完整拷贝」，默认 80，**一般不用动**。`--no-split` 关掉合并拆分、`--no-gap-fill` 关掉纯漏检填补，只在想只看某一条路径时用。`--repeat-out` 是 RepeatMasker 的 .out 文件，用于排除「真转座子区」的假基因；默认自动找 BRAKER 产物，**一般不用手动给**。`--exclude-te-gap` 开启后会把落在 TE 区的候选基因也过滤掉(默认不排，因为真基因也可能在 TE 区)。

### 生物学质控 | Biological QC

**通俗理解|In plain words:** 这是查漏补缺的「打假」三道关，默认全开，**几乎不用动**。① 真实完整 ORF(`--no-real-orf` 关闭)：候选基因必须能从头到尾翻译成完整蛋白，堵住片段；② 坐标零重叠(`--no-coord-zero-overlap` 关闭)：跟已有 BRAKER 基因坐标有交集的，不算新基因；③ 唯一比对表达(`--no-unique-reads` 关闭)：只统计「唯一比对」的 reads 算表达，多比对的 reads 多是重复区假象。`--min-expression-depth`(默认 1.0)与 `--min-coverage-breadth`(默认 50)是表达的下限，调低会更宽松(召回更多但假阳性变多)。

### 小蛋白回收通道 | Small-protein recovery lane

**通俗理解|In plain words:** 默认关闭。开启 `--recover-small-proteins` 后，会放宽长度限制，专门找回被最小 CDS 长度(默认 300bp)整类丢掉的小蛋白(如疫霉效应子)。这一组参数只在开启该通道后才有意义：`--small-max-cds-len`(450bp=150 个氨基酸)是「多短算小蛋白」的上限；identity/coverage 阈值在有表达证据时放宽到 50/50；`--small-strong-homology-identity`(95)是「强同源直通」——相似度极高(近乎自比对)的命中直接绕过 TE/表达过滤。**普通物种不用开，专门做效应子/小蛋白时才需要。**

## 分析流程 | Pipeline { #pipeline }

```text
阶段1: BRAKER 注释(内部自动完成, 断点续传)
  1) 重复序列屏蔽(RepeatModeler 建库 + RepeatMasker 屏蔽)
  2) 转录组比对(三代 minimap2 / 二代 HISAT2)
  3) BRAKER3 基因结构预测 -> braker.gtf + braker.gff3
        |
        v
阶段2: 同源查漏补缺(annorefine 核心)
  1) miniprot 证据扫描: 近缘蛋白回贴到未屏蔽基因组
  2) 漏检/合并分析:
     - 纯漏检填补: 找 BRAKER 完全没预测到的新基因
     - 合并拆分: 找 BRAKER 折叠的基因并拆开
  3) 三重质控: 完整ORF + 同源相似度/覆盖度 + 表达(深度/广度) [+ TE排除]
  4) 建基因模型 -> gap_filled.gff3
  5) 与 braker.gff3 合并 -> merged.gff3(最终结果)
```

## 输出 | Output { #output }

```text
out/
├── 01_repeat_masking/               # BRAKER 阶段: 屏蔽后的基因组
│   └── genome.fa.masked             # 屏蔽后基因组(软屏蔽)
├── 02_long_reads/                   # 三代转录本比对产物(给了 isoseq 才有)
│   └── isoseq.sorted.bam
├── 03_short_reads/                  # 二代 RNA-seq 比对产物(给了 rnaseq-dirs 才有)
│   └── rnaseq.sorted.bam            # 也用作查漏补缺的表达证据
├── 04_braker_annotation/            # BRAKER 预测结果
│   ├── braker.gtf                   # 基因结构(GTF 格式)
│   ├── braker.gff3                  # 基因结构(GFF3 格式)
│   ├── braker.aa                    # 预测蛋白序列
│   ├── braker.codingseq             # 预测 CDS 序列
│   └── hintsfile.gff                # 证据提示文件
├── logs/                            # BRAKER 阶段日志
└── 05_gap_filling/                  # 查漏补缺结果(本模块核心产物)
    ├── 00_pipeline_info/
    │   └── software_versions.yml    # 软件版本与关键参数
    ├── 01_evidence_scan/
    │   └── <prefix>.miniprot.gff3   # 蛋白证据扫描结果
    ├── 02_gap_analysis/
    │   └── <prefix>.gap_report.tsv  # 每个补入基因的证据验证报告
    ├── 03_gap_filled/
    │   └── <prefix>.gap_filled.gff3 # 只含查漏补缺新增的基因
    ├── 04_merged/
    │   └── <prefix>.merged.gff3     # 最终整合结果(BRAKER + 新增)
    └── 99_logs/
        └── annorefine.log           # 查漏补缺日志
```

## 结果解读 | Interpreting Results { #results }

### 1. 最终结果(`04_merged/<prefix>.merged.gff3`)

**通俗理解|In plain words:** 这是要交给下游分析的最终答案——BRAKER 的预测加上查漏补缺新增的基因，合在一个文件里。基因 ID 前缀是 `<prefix>_gap_N`(新增常规基因)或 `<prefix>_small_gap_N`(新增小蛋白)。

- 用 `<prefix>.gap_filled.gff3` 单独看「这次补了哪些基因」
- 用 `<prefix>.merged.gff3` 做下游分析(功能注释、比较基因组等)

### 2. 证据验证报告(`02_gap_analysis/<prefix>.gap_report.tsv`)

**通俗理解|In plain words:** 每个补入基因的「简历」——它靠什么证据被认定是真基因。一列一个证据，方便人工复核可疑的条目。

关键列：`gap_id`(基因 ID)、坐标、`cds_bp`(CDS 长度)、`prot_identity`(蛋白相似度 %)、`prot_coverage`(覆盖度 %)、`rnaseq_mean_depth`(RNA-seq 平均深度)、`coverage_breadth`(覆盖广度 %)、`fpkm`/`tpm`(表达量)、`te_overlap_pct`(TE 重叠 %)、`te_family`(TE 家族)。

- 有表达(深度 > 0 且广度达标)的条目可信度高
- `te_overlap_pct` 高且没有表达、没有强同源的条目要警惕是转座子假象

### 3. 漏检/合并统计

**通俗理解|In plain words:** 日志(以及报告行数)告诉你「补了多少、拆了多少」——补得多说明 BRAKER 漏得多(常见于注释困难物种)，拆得多说明折叠严重。数量异常大或异常小都值得回头检查参数。

## 参数选择建议 | Parameter Guidance { #guidance }

- **一般物种**：只给最小四参数 + `--rnaseq-dirs`，其余全默认
- **疫霉/卵菌**：直接用 `braker4phyto`(默认不屏蔽重复)；用 annorefine 做卵菌时考虑 `--skip-repeat`
- **专门找效应子/小蛋白**：加 `--recover-small-proteins`，必要时 `--no-small-exclude-te`(效应子常在 TE 区)、`--small-min-expression-depth 0.1`(低表达)
- **只想拆合并、不想补漏**：加 `--no-gap-fill`
- **换参数重跑**：先删旧的 `04_braker_annotation/braker.gtf` 或 `05_gap_filling` 下对应产物，否则断点续传会复用旧结果(见 FAQ)

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-g, --genome` | 必填 |  | 未mask原始基因组｜Unmasked genome |
| `-s, --species` | 必填 |  | 物种名｜Species name |
| `-p, --prot-seq` | 必填 |  | 近缘蛋白(文件或目录)｜Protein file/dir |
| `-o, --output-dir` | 必填 |  | 输出目录｜Output dir |
| `--rnaseq-dirs` | — |  | 二代RNA-seq目录(逗号分隔)｜RNA-seq dirs |
| `--isoseq` | — |  | 三代转录本(文件或目录)｜Iso-seq file/dir |
| `-t, --threads` | `12` | int | 线程数｜Threads |
| `--fungus/--no-fungus` | `True` |  | 真菌模式(疫霉适用)｜Fungus mode |
| `--no-singularity` | — |  | 不用Singularity｜No singularity |
| `--skip-repeat` | — |  | 跳过repeat屏蔽｜Skip repeat masking |
| `--skip-repeat-filter` | — |  | 跳过repeat库过滤｜Skip repeat filter |
| `--skip-rescue/--no-skip-rescue` | `True` |  | 证据还原(默认关)｜Rescue (default off) |
| `--split-min-copy-coverage` | `80` | float | 保守合并判据:完整拷贝覆盖率%｜Split copy coverage |
| `--no-split` | — |  | 关闭合并拆分｜Disable split |
| `--repeat-out` | — |  | RepeatMasker .out(filling真TE排除)｜RepeatMasker out |
| `--exclude-te-gap` | — |  | 质控排除TE区gap(默认不排)｜exclude TE-overlap gaps |
| `--no-real-orf` | — |  | 关闭真实完整ORF检查(默认开)｜disable real-ORF check (default on) |
| `--no-coord-zero-overlap` | — |  | 关闭gap坐标零重叠(默认开)｜disable coord-zero-overlap (default on) |
| `--no-unique-reads` | — |  | 关闭唯一比对过滤(默认开)｜disable unique-read filter (default on) |
| `--min-unique-mapq` | `20` | int | 唯一比对MAPQ兜底阈值｜unique MAPQ fallback |
| `--min-expression-depth` | `1.0` | float | 唯一reads平均深度下限｜min unique-read depth |
| `--min-coverage-breadth` | `50.0` | float | CDS覆盖广度%下限｜min coverage breadth |
| `--no-gap-fill` | — |  | 关闭纯漏检填补(只保留合并拆分)｜disable pure gap-fill (split only) |
| `--recover-small-proteins` | — |  | 开启小蛋白回收通道(默认关, 通用)｜enable small-protein lane (default off) |
| `--small-max-cds-len` | `450` | int | 小蛋白CDS上限bp｜small max CDS len |
| `--small-min-identity` | `50.0` | float | 小蛋白放宽identity%(有表达时)｜small min identity |
| `--small-min-coverage` | `50.0` | float | 小蛋白放宽coverage%(有表达时)｜small min coverage |
| `--small-min-expression-depth` | `1.0` | float | 小蛋白表达深度下限(effector低表达可调低如0.1)｜small min expression depth |
| `--small-min-coverage-breadth` | `60.0` | float | 小蛋白CDS覆盖广度%下限｜small min coverage breadth |
| `--no-small-exclude-te` | — |  | 关闭小蛋白TE区排除(effector常在TE区可关)｜disable small-protein TE exclusion |
| `--small-strong-homology-identity` | `95.0` | float | 强同源直通identity%阈值(≥此值绕过TE/表达过滤)｜strong-homology bypass identity |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-g, --genome` | 必填 |  | 未mask原始基因组(braker 内部 mask, filling 用未mask)｜Unmasked genome |
| `-s, --species` | 必填 |  | 物种名(braker 输出命名)｜Species name |
| `-p, --prot-seq` | 必填 |  | 近缘蛋白(文件或目录, braker+filling 共用)｜Protein file/dir |
| `-o, --output-dir` | 必填 |  | 输出目录｜Output dir |
| `--rnaseq-dirs` | — |  | 二代RNA-seq目录(逗号分隔)｜RNA-seq dirs |
| `--isoseq` | — |  | 三代转录本(文件或目录)｜Iso-seq file/dir |
| `-t, --threads` | `12` | int | 线程数｜Threads (default 12) |
| `--fungus` | `True` |  | 真菌模式(默认开, --no-fungus 关)｜Fungus mode (default on) |
| `--singularity-image` | `~/software/singularity/braker3_devel.sif` |  | Singularity镜像｜Singularity image |
| `--no-singularity` | — | store_true | 不用Singularity｜No singularity |
| `--skip-repeat` | — |  |  |
| `--skip-repeat-filter` | — | store_true | 跳过repeat库过滤(默认开)｜Skip repeat filter |
| `--skip-rescue` | `True` |  | 跳过证据还原(默认关, --no-skip-rescue 开)｜Skip rescue (default on) |
| `--split-min-copy-coverage` | `80` | float | 保守合并判据:完整拷贝覆盖率%%｜Split copy coverage (default 80) |
| `--no-split` | — | store_true | 关闭合并拆分｜Disable merged-gene split |
| `--repeat-out` | — |  | RepeatMasker .out(默认自动找braker产物)｜RepeatMasker out |
| `--exclude-te-gap` | — | store_true | 质控排除TE区gap(默认不排)｜exclude TE-overlap gaps |
| `--gap-min-identity` | `70` | float | filling identity%%(default 70) |
| `--gap-min-coverage` | `80` | float | filling coverage%%(default 80) |
| `--no-real-orf` | — | store_true | 关闭真实完整ORF检查(ATG+stop+3倍数,默认开)｜disable real-ORF check (default on) |
| `--no-coord-zero-overlap` | — | store_true | 关闭gap坐标零重叠(默认开:与BRAKER基因坐标相交不算新基因)｜disable coord-zero-overlap (default on) |
| `--no-unique-reads` | — | store_true | 关闭唯一比对过滤(默认开:多比对reads不算表达)｜disable unique-read filter (default on) |
| `--min-unique-mapq` | `20` | int | 唯一比对MAPQ兜底阈值(samtools无-e时)｜unique MAPQ fallback (default 20) |
| `--min-expression-depth` | `1.0` | float | 唯一reads平均深度下限(>0)｜min unique-read depth (default 1.0) |
| `--min-coverage-breadth` | `50.0` | float | CDS被唯一reads覆盖广度%%下限｜min coverage breadth (default 50) |
| `--no-gap-fill` | — | store_true | 关闭纯漏检填补(只保留合并拆分)｜disable pure gap-fill (split only) |
| `--recover-small-proteins` | — | store_true | 开启小蛋白回收通道(默认关, 放宽长度找回短蛋白)｜enable small-protein lane (default off) |
| `--small-max-cds-len` | `450` | int | 小蛋白CDS上限bp(默认450=150aa)｜small max CDS len (default 450) |
| `--small-min-identity` | `50.0` | float | 小蛋白放宽identity%%(默认50, 有表达时)｜small min identity (default 50, with expr) |
| `--small-min-coverage` | `50.0` | float | 小蛋白放宽coverage%%(默认50, 有表达时)｜small min coverage (default 50, with expr) |
| `--small-min-expression-depth` | `1.0` | float | 小蛋白表达深度下限(默认1.0; effector低表达可调低如0.1)｜small min expression depth (default 1.0) |
| `--small-min-coverage-breadth` | `60.0` | float | 小蛋白CDS覆盖广度%%下限(默认60)｜small min coverage breadth (default 60) |
| `--no-small-exclude-te` | — | store_true | 关闭小蛋白TE区排除(默认排; effector常在TE区可关)｜disable small-protein TE exclusion (default on) |
| `--small-strong-homology-identity` | `95.0` | float | 强同源直通identity%%阈值(默认95, ≥此值绕过TE/表达过滤)｜strong-homology bypass identity (default 95) |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies { #dependencies }

- **BRAKER3**(含 GeneMark-ETP、AUGUSTUS、TSEBRA)：通过 Singularity 镜像 `~/software/singularity/braker3_devel.sif` 运行
- **RepeatModeler / RepeatMasker / BuildDatabase**：conda 环境 `repeat`
- **minimap2**：conda 环境 `align`(三代转录本比对)
- **HISAT2**：conda 环境 `rna`(二代转录本比对)
- **miniprot / StringTie**：conda 环境 `annot`
- **samtools**：默认 `~/miniforge3/envs/annot/bin/samtools` 或系统自动检测
- **Singularity**：`~/miniforge3/envs/singularity_v.3.8.7/bin/singularity`

## 常见问题 | FAQ { #faq }

**Q1：中断后重跑会从头再来吗？**
不会。BRAKER 阶段以 `braker.gtf` 存在性判断完成(存在则整段跳过)，证据扫描以 `miniprot.gff3` 存在性判断。中断重跑会自动跳过已完成的步骤。

**Q2：换参数重跑，为什么结果没变？**
断点续传按产物存在性判断。换了查漏补缺阈值(如 identity/coverage 或 split 判据)重跑前，先删除 `05_gap_filling` 下对应的旧产物(至少删 `03_gap_filled` 和 `04_merged`)，否则会复用旧结果。

**Q3：`--gap-min-identity` / `--gap-min-coverage` 在命令行里为什么没有？**
这两个阈值(默认 70 / 80)只在直接调用 `python -m biopytools.annorefine` 时暴露，click 包装器用默认值。需要改它们时走模块直调。

**Q4：没有 RNA-seq，查漏补缺还准吗？**
能跑，但会退化为「只靠 ORF + 坐标 + 同源相似度」判断，表达过滤不生效，假阳性会多一些。强烈建议提供 RNA-seq。

**Q5：为什么要求未屏蔽的基因组？**
查漏补缺要在原始序列上找回被 BRAKER 屏蔽阶段盖掉的基因；屏蔽后的基因组会把重复区的真基因一起盖掉。

**Q6：miniprot 命中了但被质控删掉，正常吗？**
正常。质控三层(完整 ORF、同源相似度/覆盖度、表达)就是用来删假阳性命中的。想放宽可在报告里看具体是哪一关没过，再针对性调参。
