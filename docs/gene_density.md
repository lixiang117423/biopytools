# 基因密度计算 | Gene Density Calculation

一句话理解：**把每条染色体切成等长的小窗口，数每个窗口里有多少基因、折算出「每 Mb 几个基因」**——用一张表加一张图，直观看出基因组哪些区域基因密集、哪些区域是「荒漠」。

## 功能概述 | Overview { #overview }

- 按固定大小窗口沿染色体滑动，统计每窗口基因数与基因密度（genes/Mb）
- 每个基因按其「中点」归入一个窗口，不重复计数
- 输出 TSV 密度表 + 每染色体一张子图的 PNG 密度图（matplotlib）
- 染色体长度可选 .fai 或 FASTA 提供，缺省回退到 GFF 最大基因坐标
- 纯 Python 实现，无外部比对软件

## 快速开始 | Quick Start { #quick-start }

```bash
biopytools gene-density -i genome.gff3 -g genome.fa.fai -w 100000
```

最小输入：一个 GFF3 注释文件（`-i`）；`-g` 和 `-w` 都可省。

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语<br>Term | 通俗理解 |
|---|---|
| 窗口 | 把染色体切成的一段等长区间，像把一条路按固定米数分段 |
| 基因密度 | 单位长度里的基因个数（这里用「每 Mb 几个基因」），衡量基因挤不挤 |
| 中点 | 基因起止坐标的中间位置，一个基因只算在它中点所在的那个窗口里 |
| .fai 文件 | 记录每条染色体长度的索引文件（samtools faidx 生成），用来算染色体真实长度 |
| feature 类型 | GFF 第 3 列的类型标签，默认统计 `gene` 行 |

## 输入 | Input { #input }

### GFF3 注释（必需，`-i`）

标准 GFF3，统计第 3 列为 `--feature-type`（默认 `gene`）的行。第 1、4、5 列分别取染色体名、起始、终止。

### 染色体长度（可选，`-g`）

`.fai` 文件（samtools faidx 生成）或 FASTA，用于确定每条染色体总长。不给的话，染色体长度回退为「该染色体上最后一个基因的终止坐标」，末尾窗口可能不完整。

## 参数说明 | Parameters { #parameters }

### 窗口与统计对象 | Window & feature

**通俗理解|In plain words:** `-w` 是窗口大小（bp），调小图更细但更碎、调大更平滑但会糊掉局部细节，默认 10 万 bp 是常用尺度；`--feature-type` 决定统计 GFF 里哪类行，默认 `gene`，要统计别的（如 exon）改这里。

相关参数：`-w/--window-size`（默认 100000）、`--feature-type`（默认 gene）。

### 染色体长度与绘图 | Length & plot

**通俗理解|In plain words:** `-g` 给 .fai/FASTA 能让末尾窗口长度准确（否则最后一截按最大基因坐标估算）；`--no-plot` 关掉画图、只要数据表；`--prefix` 自定义输出文件前缀，默认取 GFF 文件名。

相关参数：`-g/--genome`、`--no-plot`、`--prefix`。

## 分析流程 | Pipeline { #pipeline }

```text
GFF3 + 可选染色体长度(.fai/FASTA)
    │
    ▼
解析目标 feature 坐标 → 按中点分箱
    │
    ▼
按窗口算基因数 / 密度(genes/Mb)
    │
    ▼
写 TSV 表 + 画密度图 + 写版本信息
```

## 输出 | Output { #output }

```text
gene_density_output/
├── 00_pipeline_info/
│   └── software_versions.yml        # 版本与参数记录
├── 01_density/
│   └── {prefix}_gene_density.tsv    # 密度表(主结果)
├── 02_plot/
│   └── {prefix}_gene_density.png    # 每条染色体一张子图
└── 99_logs/
    └── gene_density.log             # 运行日志
```

密度表列：`chrom`、`start`、`end`、`gene_count`、`genes_per_Mb`。

## 结果解读 | Interpreting Results { #interpreting-results }

- `gene_count`：该窗口内基因数；`genes_per_Mb`：标准化密度（基因数 / 窗口实际长度 × 1e6），跨染色体可比
- 图里每条染色体一张子图，纵轴 genes/Mb，横轴位置；峰值区=基因密集区，谷底=基因稀疏/荒漠区
- 基因密集区常见于端粒附近、GC 富集区；着丝粒附近往往密度骤降，属正常现象
- 末尾窗口若明显异常（如密度异常低），多半是没给 `-g`、染色体长度被低估所致

## 参数选择建议 | Parameter Guidance { #parameter-guidance }

- 宏观趋势：默认 `-w 100000` 即可；精细结构用 10–50 kb，超长范围用 500 kb–1 Mb
- 需要准确的末端窗口：务必给 `-g genome.fa.fai`（`samtools faidx genome.fa` 生成）
- 只要数据不要图：`--no-plot`（跳过 matplotlib，加快且少一个依赖）
- 多个 GFF 对比：用 `--prefix` 区分输出文件名，避免互相覆盖

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --gff` | 必填 | Path | GFF3注释文件｜GFF3 annotation file |
| `-o, --output-dir` | `./gene_density_output` |  | 输出目录｜Output directory |
| `-w, --window-size` | `100000` | int | 窗口大小(bp)｜Window size (bp) |
| `--feature-type` | `gene` |  | 统计的GFF feature类型｜GFF feature type to count |
| `-g, --genome` | — |  | 染色体长度来源(.fai或FASTA)｜Chromosome length source (.fai or FASTA) |
| `--prefix` | — |  | 输出文件前缀(默认GFF文件名)｜Output file prefix (default: GFF stem) |
| `--no-plot` | — |  | 不绘制密度图｜Skip density plot |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --gff` | 必填 |  | GFF3注释文件｜GFF3 annotation file |
| `-o, --output-dir` | `./gene_density_output` |  | 输出目录｜Output directory |
| `-w, --window-size` | `100000` | int | 窗口大小(bp)｜Window size (bp) |
| `--feature-type` | `gene` |  | 统计的GFF feature类型(GFF第三列)｜GFF feature type to count (column 3) |
| `-g, --genome` | — |  | 染色体长度来源(.fai或FASTA,可选,提升末尾窗口精度)｜Chromosome length source (.fai or FASTA, optional, improves last-window accuracy) |
| `--prefix` | — |  | 输出文件前缀(默认GFF文件名stem)｜Output file prefix (default: GFF stem) |
| `--no-plot` | — | store_true | 不绘制密度图｜Skip density plot |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies { #dependencies }

- Python 3
- matplotlib（仅绘图需要，缺失时自动跳过绘图并告警）
- PyYAML（写版本信息）
- 无 conda 环境、无外部生信软件依赖

## 常见问题 | FAQ { #faq }

**Q1：会断点续传吗？**
不会。单步计算，重跑直接覆盖输出。

**Q2：为什么末尾窗口密度偏低？**
没给 `-g` 时染色体长度回退到最大基因坐标，末尾窗口可能被截断或缺失。给 `-g genome.fa.fai` 即可修复。

**Q3：密度怎么算的？**
每个基因按中点归入唯一窗口（不重复计数），密度 = 窗口内基因数 /（窗口实际长度 / 1e6）。末尾窗口若不足一个窗口长度，按实际长度折算。

**Q4：`--feature-type` 找不到记录？**
确认 GFF 第 3 列确实是目标类型（区分大小写）。日志会打印统计到的 feature 总数，为 0 时会输出空表并警告。
