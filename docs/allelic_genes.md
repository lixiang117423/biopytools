# allelic_genes - 等位基因批量鉴定 | Allelic Gene Batch Identification(已废弃)

一句话理解：**这个命令已经废弃**——它原来的活儿（批量找出两个基因组之间"同源等位"的基因对）已整体搬到 `jcvi allelic`，请直接改用 `biopytools jcvi allelic`，本命令不再维护。

## 功能概述 | Overview

- 本模块在命令注册表中被标注为 `[已废弃,请使用jcvi allelic]|[DEPRECATED, use jcvi allelic]`，**已无独立源码**（`biopytools/allelic_genes/` 与 CLI 包装文件均已移除），不要再用 `biopytools allelic-genes`
- 原有功能：读入多个样本的基因组序列(FASTA)与注释(GFF)，两两做共线性比对，自动提取并汇总**等位基因对**
- 替代命令 `biopytools jcvi allelic` 功能完全等价，且并入 JCVI 工具集统一维护（`jcvi` 还提供 mcscan/macro/micro 等共线性工具）
- 典型用途：比较两个近缘物种/品系的基因组，找出"同一个基因"在两个基因组里各自的坐标，用于等位基因鉴定、基因家族收缩扩张等分析

## 快速开始 | Quick Start

~~~bash
biopytools jcvi allelic -i data -o output
~~~

最小输入：一个目录 `data/`，里面每个样本放一对 `样本前缀.fa` 和 `样本前缀.gff`（至少 2 个样本）。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| 等位基因对 | 两个基因组里"其实是同一个基因"的那一对，像同一本书的中文版和英文版 |
| 共线性 | 一段染色体上基因的"排队顺序"在另一个基因组里也基本保持，像两段旋律音高一致 |
| 锚点(anchor) | 两基因组间被判定为同源的一对基因，是拼共线性拼图的"榫卯点" |
| C-score | 对每个锚点"周围邻居也同源"的置信度打分，越高越可信，低分会被过滤掉 |
| lifted 锚点 | 借道第三方基因组"间接"确认的同源对；不带 L 的 direct 是直接比对得到的同源对 |
| LAST / diamond | 两种序列比对引擎，本工具用它们给基因"找对应"，默认 LAST |

## 输入 | Input

输入是一个目录，每个样本必须同时提供同前缀的 FASTA 和 GFF：

~~~text
data/
├── spA.fa          # 物种A基因组序列
├── spA.gff         # 物种A基因注释
├── spB.fa
└── spB.gff
~~~

- 要求至少 **2 个样本**；默认每个样本与其余所有样本两两比较，也可用 `--pairs` 指定只比某几对
- 默认按蛋白模式处理（从 GFF 的 CDS 提取蛋白序列再比对）；若注释是核酸 CDS 可用 `--dbtype nucl`

## 参数说明 | Parameters

### 基础参数 | Basic

**通俗理解|In plain words:** `-i` 是输入目录、`-o` 是输出目录，这两个必填。`-t` 是线程数（默认 24），机器核多可以调大加速，核少就保持默认即可。

### 序列类型与比对引擎 | dbtype & align engine

**通俗理解|In plain words:** `--dbtype` 决定比的是蛋白还是核酸：`prot`（默认）用蛋白序列比对，容错性好、适合远缘物种；`nucl` 用 CDS 核酸比对，只适合关系很近的样本。`--align-soft` 选比对引擎：`last`（默认）快而省内存，`diamond_blastp` 更敏感但只支持蛋白模式。**绝大多数情况这两个都不用动。**

### 过滤与配对 | Filter & pairs

**通俗理解|In plain words:** `--cscore` 是共线性置信度门槛（默认 0.7），调高更严格（结果更少但更可信）、调低更宽松（召回更多但可能混入假阳性），一般不用动。`--pairs` 用来只算指定的样本对（如 `--pairs "A,B A,C"`），样本很多、只想比其中几对时才需要，不指定就是全两两。

## 分析流程 | Pipeline

**通俗理解|In plain words:** 先给每个样本"翻译"出蛋白/CDS 序列并整理成床文件，再两两比对、过滤、扫共线性，最后从共线性锚点里挑出等位基因对并汇总。

~~~text
输入目录(.fa + .gff, 每个样本一对)
    │
    ▼
步骤1: 发现样本(前缀配对, 至少2个)
步骤2: gffread 提取蛋白/CDS序列 → 01_pep/
步骤3: GFF → BED + 去重 → 02_bed/
步骤4: 两两 LAST/diamond 比对 → .last
步骤5: blastfilter 过滤(按 C-score) → .last.filtered
步骤6: synteny scan 扫共线性 → .anchors + .lifted.anchors
    │
    ▼
步骤7: 提取等位基因对 → {A}_{B}.allelic_pairs.txt
步骤8: 汇总全部样本对 → all_allelic_pairs.txt
~~~

## 输出 | Output

~~~text
output/
├── 01_pep/                              # 每个样本提取的蛋白(.pep)或CDS(.cds)
│   ├── spA.pep
│   └── spB.pep
├── 02_bed/                              # 基因坐标 BED(去重后 .uniq.bed)
├── 03_pairwise/
│   └── spA_vs_spB/
│       ├── spA.spB.last                 # LAST 原始比对
│       ├── spA.spB.last.filtered        # C-score 过滤后
│       ├── spA.spB.anchors              # 共线性锚点
│       ├── spA.spB.lifted.anchors       # 提升后的锚点
│       └── spA_spB.allelic_pairs.txt    # 该样本对的等位基因对
├── all_allelic_pairs.txt                # 全部样本对汇总(核心结果)
└── 99_logs/
    └── allelic_genes.log                # 运行日志
~~~

- `all_allelic_pairs.txt`：最终汇总表，前面多两列 `sample_a`、`sample_b` 标明来源，后面是各样本对的等位基因对明细
- 各样本对的 `*.allelic_pairs.txt`：该对的等位基因对清单

## 结果解读 | Interpreting Results

**通俗理解|In plain words:** 核心文件 `all_allelic_pairs.txt` 每行是一对"同一个基因"，直接看 `gene_a` 和 `gene_b` 两列就是两个基因组里互相对应的基因 ID。

各列含义：

| 列 | 含义 |
|----|------|
| sample_a / sample_b | 两个样本名 |
| gene_a / gene_b | 两个基因组里互为等位的基因 ID |
| chrom_a / start_a / end_a / strand_a | 基因 a 的染色体、起止坐标、方向 |
| chrom_b / start_b / end_b / strand_b | 基因 b 的染色体、起止坐标、方向 |
| score | 共线性置信度打分，越高越可信 |
| pair_type | `direct`=直接比对得到；`lifted`=借道提升得到，可信度略低 |

- 判断好坏：`pair_type=direct` 且 score 较高的行最可信；大量 `lifted` 行说明两基因组亲缘关系较远或注释质量一般
- 若某样本对"未找到等位基因"，会记录 warning 且该对不出现在汇总里（正常，可能确实无共线性）

## 参数选择建议 | Parameter Guidance

- **近缘物种/品系**：默认 `prot` + `last` + `cscore 0.7` 即可，几乎不用改
- **关系很远的物种**：保持 `prot`（蛋白比对更稳），必要时把 `--align-soft` 换成 `diamond_blastp` 提高敏感性
- **只关心某几对样本**：用 `--pairs "A,B A,C"` 跳过无关组合，省时间
- **注释是核酸 CDS 而非蛋白**：加 `--dbtype nucl`
- **换 conda 环境**：在 `biopytools jcvi` 上用 `--conda-env` 指定（默认 `JCVI_v.1.5.6`）

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

_未找到 CLI 参数定义|No CLI definitions found_

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- JCVI 工具集（conda 环境默认 `JCVI_v.1.5.6`，可用 `biopytools jcvi --conda-env` 覆盖）
- gffread（从 GFF 提取蛋白/CDS）
- LAST（默认比对引擎）；可选 diamond（`--align-soft diamond_blastp`）

## 常见问题 | FAQ

**Q1：为什么运行 `biopytools allelic-genes` 报错/无反应？**
本模块已废弃，源码已移除。请改用 `biopytools jcvi allelic`，参数与用法见本文档。

**Q2：会不会断点续传？**
会。JCVI 管道按"输出文件已存在且非空"跳过已完成步骤（如已生成的 .pep、.bed、.last.filtered、.anchors 等）。中途失败重跑会自动续上，不会从头重算。

**Q3：提示"至少需要2个有效样本"？**
说明输入目录里配对的 `.fa` + `.gff` 不足 2 对。检查文件名前缀是否一致（必须是 `同前缀.fa` 和 `同前缀.gff`）。

**Q4：想自定义 conda 环境怎么做？**
`biopytools jcvi --conda-env 你的环境名 allelic -i data -o output`。

**Q5：`lifted` 和 `direct` 有什么区别？**
`direct` 是两基因组直接比对得到；`lifted` 是 synteny scan 借道提升得到的锚点，通常是补充找回的同源对，可信度稍低，解读时优先看 `direct`。
