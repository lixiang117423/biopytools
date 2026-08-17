# 泛基因组变异分析 | Pan-genome Variant Analysis

一句话理解：把「多个基因组构成的泛基因组」里的变异（SNP、小片段插入缺失、大片段结构变异）全部列出来并分类统计，再画一条「基因组数越多、泛基因组总大小怎么涨」的曲线，算出 gamma 值，判断这个物种的泛基因组是「还在扩张」还是「已经饱和」。

## 功能概述 | Overview { #overview }

- 输入一个泛基因组图（GBZ）或变异表（VCF），自动识别类型：`.gbz` 先做 `vg deconstruct` 提取变异，`.vcf` 直接进入统计
- 把每条变异按长度分成三类：SNP、InDel（插入缺失，长度差 50bp 以内）、SV（结构变异，长度差 >50bp）
- 输出一张逐变异的分类统计表，含每个样本实际取到的等位基因
- 用 R/ggplot2 画泛基因组增长曲线并拟合 Heaps' Law，给出 gamma 值（判断泛基因组开放/闭合）
- 断点续传：仅 deconstruct 步骤（VCF 已存在则跳过），统计与绘图每次重跑

## 快速开始 | Quick Start { #quick-start }

```bash
biopytools panvar -i graph.gbz -P T2T -o output/
```

最小输入：一个 GBZ 图（配 `-P` 参考路径前缀）或一个 VCF 文件。

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语 | 通俗理解 |
|------|----------|
| 泛基因组 | 一个物种「所有个体」基因的总和；有的基因人人都有，有的只有部分个体有 |
| 变异图 / GBZ | 把多条基因组「揉成一张图」，相同的地方共用一条路，不同的地方开叉；GBZ 是这种图的压缩存储格式 |
| VCF | 记录「哪里和参考不一样」的标准表格，一行一个变异 |
| SNP | 单个碱基的差异，像一句话里改了一个字 |
| InDel | 插入或缺失一小段（长度差 50bp 以内），像一句话里多/少了几个字 |
| SV | 结构变异，大片段（长度差 >50bp）的插入/缺失等 |
| deconstruct | 把「图」反向解包成「表格」的过程，即从变异图导出 VCF |
| gamma 值 | Heaps' Law 的指数：<1 泛基因组趋于「闭合」（加新基因组不再涨），>1 趋于「开放」（持续涨） |
| 置换（permutation） | 把样本顺序随机打乱很多遍重算，得到曲线的波动范围（阴影带） |

## 输入 | Input { #input }

支持两种输入，程序按文件后缀自动判断（`.gbz` 结尾算图，其余算 VCF）：

### GBZ 图（`.gbz`）

- 需同时给 `-P/--ref-path` 指定参考路径前缀（如 `T2T`），`vg deconstruct` 用它定位参考路径
- 程序先执行 `vg deconstruct` 导出 VCF，再做统计

### VCF 文件（`.vcf` / `.vcf.gz`）

- 直接进入统计，跳过 deconstruct
- VCF 里必须有样本（header 含 `#CHROM ... FORMAT <样本列>`），否则打印警告并跳过统计

## 参数说明 | Parameters { #parameters }

### 必需参数 | Required

**通俗理解|In plain words:** 告诉程序「读哪个文件、结果放哪」。GBZ 输入时还缺一个「参考路径前缀」，就像解包一个打包文件时要先知道从哪条主线开始还原，缺了它 deconstruct 无法定位参考路径，程序会直接报错。

### 运行资源 | Runtime

**通俗理解|In plain words:** 线程数管「用几个 CPU 并行」，调大通常更快，但 vg 与 R 的并行收益有限，默认 12 一般够用。`--vg-env` 与 `--r-path` 是软件位置设置，只在你的 vg / Rscript 不叫默认名字时才需要改，一般不用动。

### 统计与绘图 | Statistics & plotting

**通俗理解|In plain words:** `--ref-size` 是参考基因组大小（Mb），默认 0 让程序从 VCF 坐标自动推算；只有当自动推算明显偏小（比如只统计了部分染色体）时才手动指定。`--permutations` 决定增长曲线随机置换多少次，次数越多曲线阴影带越稳、越慢，默认 100 是快慢折中，一般不用动。

## 分析流程 | Pipeline { #pipeline }

```text
输入 GBZ / VCF
  -> [GBZ] vg deconstruct 导出 VCF（已存在则跳过）
  -> pysam 逐条分类统计（SNP / InDel / SV）
  -> 生成 pangenome_variants_summary.tsv
  -> R/ggplot2 随机置换 + 拟合 Heaps' Law
  -> 生成增长曲线图 + gamma 值
```

## 输出 | Output { #output }

```text
output_dir/
├── 01_deconstruct/
│   └── <输入名>.vcf                    # 仅 GBZ 输入时生成，deconstruct 导出的变异
├── 02_summarize/
│   ├── pangenome_variants_summary.tsv  # 逐变异分类统计表（核心结果）
│   ├── pangenome_growth_curve.png      # 泛基因组增长曲线图
│   └── plot_growth.R                   # 绘图脚本，可复现/改参重画
└── 99_logs/
    └── panvar.log                      # 运行日志
```

## 结果解读 | Interpreting Results { #interpreting-results }

### 1. 变异统计表（pangenome_variants_summary.tsv）

每一行是一条变异，关键列：

- `TYPE`：SNP / InDel / SV（分类依据见概念表）
- `LENGTH`：变异长度差（`len(alt) - len(ref)`），SV 时绝对值较大
- `ALT_COUNT` / `ALT_ALLELES`：有几个不同的等位基因、分别是什么
- `FREQ_DETAIL`：每个 ALT 等位基因在样本里的出现次数与占比（如 `ALT1:3(60.0%)`）
- 后续各列是每个样本实际取到的等位基因，`.` 表示缺失

程序日志里会汇总打印总变异数，以及 SNP / InDel / SV 各有多少。

### 2. 增长曲线图（pangenome_growth_curve.png）

- 横轴 = 基因组数量，纵轴 = 泛基因组总大小（Mb），点与线是均值，蓝色阴影带是随机置换的波动范围
- 副标题里有关键数字：参考基因组大小、`gamma = ...`、`R-squared = ...`
- **gamma < 1**：曲线趋于平缓，泛基因组「闭合」，核心序列已基本饱和
- **gamma > 1**：曲线仍明显上扬，泛基因组「开放」，新样本还会持续贡献新序列
- R-squared 越接近 1，拟合越可信

## 参数选择建议 | Parameter Guidance { #parameter-guidance }

- GBZ 输入一定记得 `-P`；VCF 输入不需要，给了也无妨
- 快速预览可 `--permutations 20`，正式结果用默认 100，追求平滑可加到 500
- 自动推断的参考大小明显不对（如只测了少量 contig）时，用 `--ref-size` 手动指定（Mb）
- `--threads`、`--vg-env`、`--r-path` 绝大多数场景保持默认即可

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | [FILE] GBZ图文件或VCF文件｜GBZ graph or VCF file |
| `-o, --output` | 必填 |  | [DIR] 输出目录｜Output directory |
| `-P, --ref-path` | — |  | [STR] 参考路径前缀，GBZ输入时必需 (如T2T)｜Reference path prefix (required for GBZ) |
| `-t, --threads` | `12` | int | [INT] 线程数 (default: 12) |
| `--ref-size` | `0.0` | float | [FLOAT] 参考基因组大小Mb (0=自动推断, default: 0) |
| `--permutations` | `100` | int | [INT] 增长曲线随机置换次数 (default: 100) |
| `--vg-env` | — |  | [STR] vg conda环境名 (default: vg_v.1.7.0) |
| `--r-path` | — |  | [FILE] Rscript路径｜Rscript binary path |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies { #dependencies }

- vg（conda 环境 `vg_v.1.7.0`，通过 `--vg-env` 可改）
- Rscript + R 包 ggplot2（通过 `--r-path` 指定 Rscript）
- Python 包 pysam（读 VCF）

## 常见问题 | FAQ { #faq }

**Q1：GBZ 输入报「需要参考路径前缀」？**
GBZ 输入必须给 `-P/--ref-path`，否则 `vg deconstruct` 不知道以哪条参考路径为基准，程序会在校验阶段直接报错。

**Q2：断点续传怎么生效？**
只有 deconstruct 步骤会跳过（VCF 已存在且非空时）；统计表和增长曲线每次都会重跑。所以 VCF 输入时实际上每次都是全新统计。

**Q3：R 报「there is no package called ggplot2」？**
Rscript 所在环境没装 ggplot2，需在对应 R 环境里安装（`install.packages('ggplot2')`），或换一个装了 ggplot2 的 Rscript 路径。

**Q4：VCF 里没有样本会怎样？**
程序打印警告并跳过统计，不产出 summary.tsv。确认 VCF header 里有样本列。

**Q5：参考大小自动推断为什么有时是 87.0 Mb？**
当从坐标推算出的参考大小 < 1.0 Mb（如只统计了很短区间）时，程序用 87.0 Mb 兜底，避免除出离谱结果；需要精确时请用 `--ref-size` 手动指定。
