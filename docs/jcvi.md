# JCVI 共线性分析工具集 | JCVI Synteny Analysis Toolkit (MCscan)

一句话理解：**一套"基因排列顺序"比对与可视化工具，帮你找出多个物种基因组里"顺序相似的同源基因块"，并画成宏观（染色体核型）或微观（局部基因）共线性图**。

## 功能概述 | Overview

- 四个子命令：`mcscan`（两两共线性分析）、`allelic`（等位基因批量鉴定）、`macro`（宏观核型可视化）、`micro`（微观局部可视化）
- 共享管道：发现样本 → gffread 提取序列 → GFF 转 BED + 去重 → LAST/blast 比对 → blastfilter 过滤 → synteny scan 出 anchors
- 支持 prot（蛋白，默认）/ nucl（CDS）两种 dbtype；last / diamond_blastp 两种比对软件
- 断点续传：每个步骤输出文件存在且非空即自动跳过
- macro 支持多物种（3+），自动生成 seqids/layout，`--replot` 只重绘图

## 快速开始 | Quick Start

```bash
biopytools jcvi mcscan -i data -o output
```

最小输入：一个目录，里面成对放着 `*.fa` 和同名 `*.gff`（至少 2 个样本）。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| 共线性(synteny) | 两个基因组里"同源基因排列顺序相似"的保守块 |
| anchors | 两个基因组间"同源基因对"的连线锚点文件 |
| C-score | 过滤比对的质量阈值，越高越严格 |
| LAST / diamond | 序列比对软件（LAST 默认，diamond 更快） |
| 核型(karyotype) | 染色体级别的宏观排列图 |
| GFF | 基因注释文件，记录基因/转录本的坐标 |
| seqids / layout | JCVI 绘图用的染色体顺序与排版配置 |
| 等位基因(allelic) | 同一物种不同基因组间对应位置上的同源基因 |

## 输入 | Input

- `-i/--input`：输入目录，内含成对的 `*.fa`（基因组/蛋白）与同名 `*.gff`（如 `A.fa` + `A.gff`）
- macro 额外：`--species A,B,C`（有序物种列表）
- micro 额外：`--pairs A,B` + `--region-a chr:start-end` + `--region-b chr:start-end`

> 样本"前缀" = 去掉 `.gff` 的路径，程序自动要求同前缀存在对应 `.fa`。

## 参数说明 | Parameters

### 通用参数 | Common

**通俗理解|In plain words:** 四个子命令共用的核心开关。`-i` 输入目录、`-o` 输出目录、`-t` 线程数（默认 24）。`--dbtype` 选序列类型：`prot`（蛋白，默认，最常用）或 `nucl`（CDS 核酸）；`--cscore` 是比对过滤阈值（默认 0.7，越高越严、block 越少）；`--align-soft` 选比对软件（`last` 默认，`diamond_blastp` 更快但仅 prot）；`--pairs` 指定只比哪些配对（如 `A,B A,C`），不填则全两两比较。**一般用默认即可**。

### macro 专属 | Macro-specific

**通俗理解|In plain words:** `--species` 指定要画哪些物种及顺序（逗号分隔，2 个或 3+ 个）。`--minspan`（默认 30）是 screen 步骤的最小跨度，调大会只保留更大的共线性块、图更简洁；`--min-chr-genes`（默认 20）过滤掉基因数过少的 scaffold，避免图里一堆碎片；`--gff-key` 是 GFF 里基因 ID 的属性键（默认 ID）。`--figsize`/`--shadestyle`/`--chrstyle` 是绘图样式。`--replot` 用已有 seqids/layout 只重绘图（改样式时省去重新比对）。**物种顺序错了会影响排版，务必按想要的顺序写**。

### micro 专属 | Micro-specific

**通俗理解|In plain words:** `--pairs` 是物种对，`--region-a`/`--region-b` 是两个基因组上要放大看的区间（`chr:start-end`）。`--genes-a`/`--genes-b` 可进一步只显示指定基因（文件，一行一个基因 ID）；`--glyph-style`（arrow/box）和 `--shadestyle` 是绘图样式。**区域跨度别太小**，否则区间内基因太少画不出共线性。

## 分析流程 | Pipeline

```text
输入目录(成对 .fa + .gff)
    │
    ▼
步骤1: 发现样本（找 .gff，校验对应 .fa）
    │
    ▼
步骤2: gffread 提取序列（prot→.pep / nucl→.cds）
    │
    ▼
步骤3: GFF 转 BED + uniq 去重
    │
    ▼
步骤4-6: 两两比较（断点续传，逐文件跳过）
    ├─ LAST/blast 比对 → .last
    ├─ blastfilter 过滤 → .last.filtered
    └─ synteny scan → .anchors / .lifted.anchors
    │
    ▼
各子命令专属步骤：
    ├─ mcscan: 汇总 blocks/gene pairs → all_collinearity_summary.txt
    ├─ allelic: 提取等位基因对 → all_allelic_pairs.txt
    ├─ macro: screen → seqids → layout → karyotype.pdf
    └─ micro: 取区间 block → 过滤 → blocks.layout → synteny 图(.pdf) + collinear_pairs.xlsx
```

## 输出 | Output

```text
output/
├── 01_pep/{样本}.pep（或 .cds）          # 提取的蛋白/CDS 序列
├── 02_bed/{样本}.bed、{stem}.uniq.bed    # GFF 转 BED 及去重结果
├── 03_pairwise/{A}_vs_{B}/
│   ├── {stemA}_{stemB}.last              # 原始比对
│   ├── {stemA}_{stemB}.last.filtered     # C-score 过滤后
│   ├── {stemA}_{stemB}.anchors           # 共线性锚点（直接比对）
│   ├── {stemA}_{stemB}.lifted.anchors    # liftover 推断后的锚点
│   └── {A}_{B}.allelic_pairs.txt         # 等位基因对（allelic 子命令）
├── all_collinearity_summary.txt          # mcscan 汇总
├── all_allelic_pairs.txt                 # allelic 汇总
├── seqids / layout + karyotype.pdf       # macro 可视化
├── micro/{…}.pdf + collinear_pairs.xlsx  # micro 可视化
└── 99_logs/*.log                         # 各子命令日志
```

## 结果解读 | Interpreting Results

### 1. anchors 文件（共线性锚点）

**通俗理解|In plain words:** 每个 block 以 `###` 分隔，block 内每行是一个同源基因对。block 越多 = 两基因组共线性块越多；block 越大 = 连续保守的基因越多。`all_collinearity_summary.txt` 已把 blocks / gene pairs / 平均与最大 block 大小汇总成表。

### 2. 汇总表（`all_collinearity_summary.txt` / `all_allelic_pairs.txt`）

**通俗理解|In plain words:** mcscan 的输出是"两两比较的共线性统计表"；allelic 的输出是"等位基因对清单"（含坐标、链、score、pair_type）。gene pairs 越多、blocks 越大，两基因组亲缘越近/保守性越强。

### 3. macro 核型图（`karyotype.pdf`）

**通俗理解|In plain words:** 多条染色体并排，彩色连线表示物种间的共线性区块。连线平行成束 = 大段共线；交叉 = 倒位/易位。

### 4. micro 共线性图（`micro/{…}.pdf`）

**通俗理解|In plain words:** 放大到基因水平，箭头表示基因及方向，连线表示同源基因。适合核对某个具体区间的基因排列是否保守。`collinear_pairs.xlsx` 是配套的配对明细表。

## 参数选择建议 | Parameter Guidance

- `--dbtype`：有蛋白注释用默认 `prot`；只有 CDS 用 `nucl`
- `--align-soft`：样本多/想快用 `diamond_blastp`（需 prot）；默认 `last` 更稳
- `--cscore`：默认 0.7；block 太少就降到 0.5-0.6，噪声太多就升到 0.8
- `--pairs`：样本多时只指定关心的配对，避免全两两耗时爆炸
- macro `--min-chr-genes`：scaffold 碎片多时调大（如 50）让图更干净
- micro 区域跨度：建议几十 kb 到几百 kb，太小会因基因太少而画不出

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

_未找到 CLI 参数定义|No CLI definitions found_

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- JCVI（conda 环境 `JCVI_v.1.5.6`，默认；可用 `--conda-env` 改）
- gffread（提取蛋白/CDS 序列）
- LAST（默认比对软件）
- 可选：diamond（`--align-soft diamond_blastp` 时）

## 常见问题 | FAQ

**Q1：换参数重跑，结果没变？**
断点续传按各步骤产物是否存在判断。改了 `--cscore`/`--dbtype`/`--align-soft` 后，要删掉对应旧产物（如 `03_pairwise/…/.last.filtered` 和 `.anchors`）再重跑，否则复用旧结果。

**Q2：报"至少需要 2 个有效样本"？**
输入目录里要成对出现 `*.fa` 和同名 `*.gff`。缺 `.fa` 的 `.gff` 会被跳过并警告。

**Q3：macro 的多物种顺序怎么定？**
`--species A,B,C` 的顺序决定图上染色体的排列，相邻物种之间画共线性。按你想要的进化顺序或展示顺序写。

**Q4：micro 报 region 格式错误？**
区间必须是 `chr:start-end`（含冒号和连字符），且染色体名要和 BED/GFF 里一致。

**Q5：--replot 报找不到 seqids/layout？**
`--replot` 依赖之前完整运行生成的 seqids/layout 文件，必须先不带 `--replot` 跑过一次。

**Q6：diamond_blastp 报错？**
diamond 只支持 prot 模式（`--dbtype prot`）；用 nucl 时只能选 `last`。
