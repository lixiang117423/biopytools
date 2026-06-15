# JCVI共线性分析工具集 | JCVI Synteny Analysis Toolkit

**MCscan 共线性 + 等位基因鉴定 + 宏观/微观共线性可视化 | MCscan collinearity, allelic identification, macro/micro synteny visualization.**

## 功能概述 | Overview

`jcvi` 模块封装了 JCVI (Python-based comparative genomics toolkit) 的核心工作流，提供四个子命令：

| 子命令 | 说明 |
|------|------|
| `mcscan`  | MCscan 批量两两共线性分析（输入 CDS/蛋白 + GFF，输出共线性块） |
| `allelic` | 基于共线性的等位基因批量鉴定 |
| `macro`   | 宏观共线性（核型）可视化，自动调用 mcscan，支持多物种 |
| `micro`   | 微观共线性（局部区块）可视化，自动调用 mcscan |

所有子命令都通过 conda 环境调用 JCVI，默认 `JCVI_v.1.5.6`，可用 `--conda-env` 覆盖。底层比对软件可选 `last` 或 `diamond_blastp`。

## 快速开始 | Quick Start

输入目录约定：每个物种提供一个 `*.fa`（蛋白或 CDS）和一个 `*.gff` 文件，文件名（去掉扩展名）作为物种名。

```bash
# 1. 批量两两 MCscan 共线性分析
biopytools jcvi mcscan -i data/ -o output/ -t 24

# 2. 等位基因鉴定
biopytools jcvi allelic -i data/ -o output/

# 3. 两物种宏观共线性核型图
biopytools jcvi macro -i data/ -o output/ --species A,B

# 4. 三物种宏观共线性
biopytools jcvi macro -i data/ -o output/ --species A,B,C --shadestyle gradient

# 5. 微观共线性（局部区块）
biopytools jcvi micro -i data/ -o output/ \
    --pairs A,B \
    --region-a chr1:10000-50000 \
    --region-b chr3:20000-80000
```

## 参数说明 | Parameters

### 公共参数（所有子命令）| Common Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-i, --input` | 必需 | 输入目录（含 `*.fa` 和 `*.gff`） |
| `-o, --output` | 必需 | 输出目录 |
| `-t, --threads` | `24` (macro/mcscan/allelic) 或 `12` (micro) | 线程数 |
| `--conda-env` | `JCVI_v.1.5.6` | JCVI 所在 conda 环境（顶层选项） |
| `--dbtype` | `prot` | 序列类型：`prot` / `nucl` |
| `--cscore` | `0.7` | C-score 阈值 |
| `--align-soft` | `last` | 比对软件：`last` 或 `diamond_blastp` |
| `--pairs` | - | 指定配对，空格分隔，如 `A,B A,C`（mcscan/allelic）；`A,B`（micro，必需） |

### macro 子命令特有 | macro-specific

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--species` | 必需 | 有序物种列表，逗号分隔，如 `A,B,C` |
| `--gff-key` | `ID` | GFF 属性键 |
| `--minspan` | `30` | 最小共线性块跨度 |
| `--min-chr-genes` | `20` | 染色体最小基因数，少于此值的 scaffold 被过滤 |
| `--figsize` | 自动 | 画布大小，如 `14x10` |
| `--shadestyle` | `line` | 阴影样式：`line`/`gradient`/`solid` |
| `--chrstyle` | `rect` | 染色体样式：`rect`/`roundrect` |
| `--replot` | off | 仅重新绘图，复用已有 seqids/layout |

### micro 子命令特有 | micro-specific

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--region-a` / `--region-b` | 必需 | 物种 A / B 区间，格式 `chr:start-end` |
| `--genes-a` / `--genes-b` | - | 高亮基因 ID 列表文件（一行一个） |
| `--glyph-style` | `arrow` | 基因 glyph 样式：`arrow`/`box` |
| `--shadestyle` | `line` | 阴影样式：`line`/`gradient`/`solid` |

（运行 `biopytools jcvi <subcommand> -h` 查看完整参数列表）

## 输出 | Output

- `mcscan`/`allelic`: 共线性块文件、LAST/DIAMOND 比对结果、filtered 同源对
- `macro`: PDF/PNG 核型图、`seqids`/`layout` 文件、`blocks`/`anchors`
- `micro`: 局部共线性 PDF、`bed`/`layout` 中间文件

所有运行日志在 `99_logs/` 下。

## 依赖 | Dependencies

- conda 环境 `JCVI_v.1.5.6`（或自定义），内含 `jcvi`、`lastal`、`diamond`
- `matplotlib`, `biopython`（JCVI 依赖）

## 引用 | Citation

- Tang H. et al. Synteny and collinearity in plant genomes. Science. 2008. (MCscan)
- Tang H. et al. JCVI: A comparative genomics analysis toolkit. (JCVI)
- Kiełbasa SM et al. Adaptive seeds tame genomic sequence comparison. Genome Res. 2011. (LAST)
- Buchfink B. et al. Fast and sensitive protein alignment using DIAMOND. Nature Methods. 2021.

## 相关链接 | References

- [JCVI GitHub](https://github.com/tanghaibao/jcvi)
- [项目主页](https://github.com/lixiang117423/biopytools)
