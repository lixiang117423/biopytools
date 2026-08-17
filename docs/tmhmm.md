# TMHMM - 跨膜螺旋预测 | TMHMM Transmembrane Helix Prediction

一句话理解：**给每条蛋白数一数「有几段横穿细胞膜」，并估计它是不是膜蛋白**。输入一个蛋白质 FASTA，输出一张汇总表，告诉你跨膜螺旋个数、位置及膜蛋白可能性。

## 功能概述 | Overview

- 基于经典 TMHMM 隐马尔可夫模型预测跨膜螺旋（TM helix）
- 输出每个蛋白的跨膜螺旋个数、跨膜区氨基酸期望数、N 端跨膜概率、信号序列提示
- 把 TMHMM 原始文本解析成干净 TSV；可选 `--plot` 生成图形
- 直接调用 tmhmm（无需 conda 包装），采用原子写入避免断点续传时读到半截文件
- 断点续传：raw 与 clean 两个输出都存在则跳过

## 快速开始 | Quick Start

```bash
biopytools tmhmm -i proteins.fa -o output_dir
```

最小输入：一个蛋白质 FASTA 文件 + 输出目录。输出前缀默认取输入文件名。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| 跨膜螺旋(TMH) | 蛋白质「横穿」细胞膜的一段螺旋，像一根钉子钉穿木板 |
| 隐马尔可夫模型(HMM) | 一种「按序列逐位打标签」的统计模型，判断每个氨基酸属于膜内/跨膜/膜外 |
| N-in 概率 | 蛋白 N 端在细胞膜「内侧」的概率，用于判断膜的朝向 |
| 信号序列提示 | TMHMM 对 N 端是否像信号肽/信号锚的提示，非严格信号肽预测 |

## 输入 | Input

蛋白质 FASTA 文件：

```text
>protein_01
MKKLLIAAMMAAALAACSQEAKTEVFSKSADEGGAPK...
```

- 序列必须是氨基酸（蛋白质）
- 输入文件名自动作为默认输出前缀

## 参数说明 | Parameters

### 必需参数 | Required

**通俗理解|In plain words:** 两个必填：输入蛋白文件、输出目录。没有默认值。

### 绘图与前缀 | Plot & prefix

**通俗理解|In plain words:** `--plot` 生成跨膜结构示意图（默认关闭，`-noplot`）；`--prefix` 决定输出文件名，默认取输入文件名，**一般不用动**。

## 分析流程 | Pipeline

**通俗理解|In plain words:** 先查 raw 和 clean 两个输出是否都在（在就收工），否则跑 tmhmm 生成原始输出，再解析成汇总表。

```text
蛋白质 FASTA
    │
    ▼
(断点续传) raw + clean 都存在? → 是则跳过
    │
    ▼
运行 tmhmm (默认 -noplot) → 原子写入 {prefix}_tmhmm_raw.txt
    │
    ▼
解析原始输出 → {prefix}_tmhmm.tsv
```

## 输出 | Output

```text
output_dir/
├── 00_pipeline_info/
│   └── software_versions.yml     # 软件版本与参数记录
├── {prefix}_tmhmm_raw.txt        # TMHMM 原始输出
├── {prefix}_tmhmm.tsv            # 主结果:解析后的汇总表
└── 99_logs/
    └── tmhmm.log                 # 运行日志
```

- `{prefix}` 默认是输入文件名
- 主结果 `{prefix}_tmhmm.tsv` 列：`ID / Length / Pred_TMHs / Exp_AAs_in_TMHs / Exp_AAs_first60 / Prob_N_in / Signal_Seq`

## 结果解读 | Interpreting Results

### 1. 汇总表（`{prefix}_tmhmm.tsv`）

**通俗理解|In plain words:** 每行一个蛋白，看 `Pred_TMHs` 一列就能定性是不是膜蛋白。

- `Pred_TMHs`：预测的跨膜螺旋个数。0 个=非跨膜蛋白；≥1 个=膜蛋白（跨膜蛋白）
- `Exp_AAs_in_TMHs`：预计落在跨膜区内的氨基酸数，越大说明跨膜部分越多
- `Exp_AAs_first60`：前 60 个氨基酸里预计在跨膜区内的个数，辅助判断 N 端是否跨膜
- `Prob_N_in`：N 端在膜内侧的概率，用于判断朝向
- `Signal_Seq`：`yes` 表示 N 端可能有信号序列（信号肽/信号锚提示），否则为空

### 2. 好坏判据 | Judgment

- **跨膜螺旋 1 个且 N 端有信号序列提示**：单次跨膜受体/信号锚定蛋白的常见模式
- **跨膜螺旋 ≥2 个**：多次跨膜蛋白（转运蛋白、通道、G 蛋白偶联受体等）
- **跨膜螺旋 0 个**：可溶性蛋白或外周膜蛋白（需注意 TMHMM 对信号肽与 N 端跨膜区偶有混淆）

## 参数选择建议 | Parameter Guidance

- **`--plot`**：需要可视化跨膜结构时开启，批量统计时关闭（默认）以省时省空间
- **`--prefix`**：同输入复用想区分结果时设置，否则用默认
- 本工具参数极少，一条最简命令即可覆盖绝大多数场景

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 输入蛋白质FASTA文件｜Input protein FASTA file |
| `-o, --output-dir` | 必填 | Path | 输出目录｜Output directory |
| `--plot` | `False` |  | 生成图形(默认-noplot)｜Generate plots (default -noplot) |
| `--prefix` | — |  | 输出文件前缀｜Output file prefix |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | [FILE] 输入蛋白质FASTA文件｜Input protein FASTA file |
| `-o, --output-dir` | 必填 |  | [DIR] 输出目录｜Output directory |
| `--plot` | `False` | store_true | [FLAG] 生成图形(默认不生成)｜Generate plots (off by default) |
| `--prefix` | — |  | [STR] 输出文件前缀(默认使用输入文件名)｜Output file prefix |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- TMHMM（`tmhmm`，默认 `~/miniforge3/envs/protein/bin/tmhmm`，直接调用，无需 conda 包装）
- Python 3 + PyYAML（记录版本）

## 常见问题 | FAQ

**Q1：TMHMM 为什么不走 conda 包装？**
TMHMM 的 Perl 驱动仅依赖核心 Getopt::Long 模块、内部 decodeanhmm 是静态链接二进制，不依赖 conda 环境的动态库，因此直接调用即可。

**Q2：换参数重跑，为什么结果没变？**
断点续传要求 `{prefix}_tmhmm_raw.txt` 和 `{prefix}_tmhmm.tsv` 都存在才跳过。想强制重跑，删掉这两个文件（或换输出目录）。

**Q3：报「tmhmm 不存在」？**
默认路径 `~/miniforge3/envs/protein/bin/tmhmm` 里没有该程序。可设环境变量 `TMHMM_PATH` 或在 `~/.config/biopytools/config.yml` 里配置 tmhmm 路径。

**Q4：TMHMM 和 Phobius 有什么区别？**
都预测跨膜螺旋。Phobius 把信号肽和跨膜区联合建模，区分二者更准；TMHMM 是经典跨膜预测工具、速度快。需要同时精确判断信号肽和跨膜时优先 Phobius；只关心跨膜螺旋个数时 TMHMM 足够。