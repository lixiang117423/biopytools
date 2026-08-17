# DeepTMHMM - 跨膜螺旋与信号肽预测 | DeepTMHMM Transmembrane Helix & Signal Peptide Prediction

一句话理解：**用深度学习给每条蛋白判断「有没有跨膜区、有没有信号肽」**，并标出这些结构在序列上的具体位置。输入一个蛋白质 FASTA，输出一张汇总表，告诉你每条蛋白是膜蛋白还是分泌蛋白。

## 功能概述 | Overview

- 基于 DeepTMHMM 1.0 深度学习模型，同时预测跨膜螺旋（TM helix）与信号肽
- 输出每个蛋白的跨膜螺旋个数与坐标区间、信号肽有无及位置
- 除了原始输出（3line 拓扑、GFF3 注释、Markdown 报告），还整理出一张干净的 TSV 汇总表
- 通过 conda 环境运行官方 predict.py，自动处理其「输出目录不能预先存在」的限制
- 断点续传：三个关键输出齐全则直接跳过，不重复计算

## 快速开始 | Quick Start

```bash
biopytools deeptmhmm -i proteins.fa -o output_dir
```

最小输入：一个蛋白质 FASTA 文件。输出前缀默认取输入文件名。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| 跨膜螺旋 | 蛋白质「横穿」细胞膜的一段螺旋结构，像一根钉子钉穿木板；有几个这样的「钉子」就叫有几个跨膜螺旋 |
| 信号肽 | 蛋白质 N 端的一小段「地址标签」，引导蛋白被分泌或定位，通常随后被切掉 |
| 拓扑(topology) | 蛋白质相对细胞膜「怎么穿」的示意图：哪段在膜内、哪段在膜外 |
| GFF3 | 一种标准注释格式，用「哪条序列的第几到第几位是什么特征」的方式记录结构位置 |

## 输入 | Input

蛋白质 FASTA 文件，支持 .fa / .faa / .fasta（含 .gz 压缩）：

```text
>protein_A
MKKLLIAAMMAAALAACSQEAKTEVFSKSADEGGAPK...
```

- 序列必须是氨基酸（蛋白质）；一条序列可包含多个跨膜螺旋
- 输入文件名会被用作默认输出前缀（自动剥离 .gz 和扩展名）

## 参数说明 | Parameters

### 必需参数 | Required

**通俗理解|In plain words:** 两个必填：输入蛋白文件、输出目录。没有默认值。

### 输出前缀 | Output prefix

**通俗理解|In plain words:** 决定输出文件叫什么名。不填就默认用输入文件名，**一般不用动**；只有想给不同批次的结果区分命名时才指定。

### 运行环境 | Runtime environment

**通俗理解|In plain words:** 这两个指向 DeepTMHMM 装在哪：conda 环境名、安装目录。部署时已配好默认值，**普通用户一般不用动**；只有你的环境名或安装位置不同时才需要指定。

## 分析流程 | Pipeline

**通俗理解|In plain words:** 先查三个关键输出是否都在（在就直接收工），否则在输出目录下建临时目录跑 predict.py，再把结果搬回并按规范改名、解析成汇总表。

```text
蛋白质 FASTA
    │
    ▼
(断点续传) summary.tsv + topologies.3line + tmr.gff3 是否齐全?
    ├─ 是 → 直接结束
    └─ 否 → 建临时目录 → conda run predict.py --fasta 输入 --output-dir 临时目录
              │
              ▼
         搬运并改名: predicted_topologies.3line / TMRs.gff3 / deeptmhmm_results.md
    │
    ▼
解析 3line + GFF3 → 生成 {prefix}_deeptmhmm_summary.tsv
```

## 输出 | Output

```text
output_dir/
├── 00_pipeline_info/
│   └── software_versions.yml                       # 软件版本与参数记录
├── {prefix}_deeptmhmm_summary.tsv                  # 主结果:干净汇总表
├── {prefix}_deeptmhmm_topologies.3line             # 每条蛋白的三行拓扑原始输出
├── {prefix}_deeptmhmm_tmr.gff3                     # 跨膜螺旋/信号肽的 GFF3 注释
├── {prefix}_deeptmhmm_results.md                   # DeepTMHMM 官方 Markdown 报告
└── 99_logs/
    └── deeptmhmm.log                               # 运行日志
```

- `{prefix}` 默认是输入文件名（如 proteins.fa → proteins）
- `{prefix}_deeptmhmm_summary.tsv` 是给用户看的主结果，列：`ID / Length / Protein_Type / Pred_TMHs / Signal_Peptide / TM_Regions`

## 结果解读 | Interpreting Results

### 1. 汇总表（`{prefix}_deeptmhmm_summary.tsv`）

**通俗理解|In plain words:** 每行一个蛋白，一看「跨膜螺旋个数」和「信号肽」两列就能定性。

- `Pred_TMHs`：预测的跨膜螺旋个数。0 个=不是跨膜蛋白；≥1 个=跨膜蛋白（膜蛋白）
- `Signal_Peptide`：`no` 或 `yes (起-止)`；有信号肽的蛋白通常是分泌蛋白或膜蛋白
- `Protein_Type`：DeepTMHMM 给出的蛋白类型标签（如 Globular、TM、SP 等）
- `TM_Regions`：各跨膜螺旋在序列上的坐标区间，多个用 `;` 分隔，无则 `-`

### 2. 好坏判据 | Judgment

- **跨膜螺旋数 1 个且 N 端有信号肽**：典型的单次跨膜受体模式；**跨膜螺旋 ≥2 个**：多次跨膜蛋白（如转运蛋白、通道）
- **只有信号肽、无跨膜区**：分泌蛋白的典型特征
- 跨膜区坐标和信号肽位置可对照 `{prefix}_deeptmhmm_tmr.gff3` 里的精确注释

## 参数选择建议 | Parameter Guidance

- **`--prefix`**：同一条输入在不同项目里复用、想区分结果时设置，否则用默认即可
- **`--conda-env` / `--deeptmhmm-dir`**：仅当部署环境不同时才动，换机器时按实际路径指定

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 输入蛋白质FASTA文件｜Input protein FASTA file |
| `-o, --output-dir` | 必填 | Path | 输出目录｜Output directory |
| `--prefix` | — |  | 输出文件前缀(默认输入文件名)｜Output file prefix |
| `--conda-env` | `deeptmhmm_v.1.0` |  | conda环境名｜conda env name |
| `--deeptmhmm-dir` | `~/software/deeptmhmm/DeepTMHMM-Academic-License-v1.0` |  | DeepTMHMM安装目录｜DeepTMHMM install directory |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | [FILE] 输入蛋白质FASTA文件｜Input protein FASTA file |
| `-o, --output-dir` | 必填 |  | [DIR] 输出目录｜Output directory |
| `--prefix` | — |  | [STR] 输出文件前缀(默认输入文件名)｜Output file prefix |
| `--conda-env` | `deeptmhmm_v.1.0` |  | [STR] conda环境名｜conda env name |
| `--deeptmhmm-dir` | `~/software/deeptmhmm/DeepTMHMM-Academic-License-v1.0` |  | [STR] DeepTMHMM安装目录｜DeepTMHMM install directory |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- DeepTMHMM 1.0 安装目录（默认 `~/software/deeptmhmm/DeepTMHMM-Academic-License-v1.0`，含 predict.py）
- conda 环境 `deeptmhmm_v.1.0`（提供运行 predict.py 的 Python 及依赖）
- Python 3 + PyYAML（记录软件版本）

## 常见问题 | FAQ

**Q1：为什么输出目录里会冒出临时目录？**
DeepTMHMM 官方 predict.py 要求输出目录「不能预先存在」，所以程序先在输出目录下建一个临时目录传给它，跑完再把结果搬出来、删除临时目录。若中途崩溃可能残留 `_deeptmhmm_run_*` 临时目录，可手动删除。

**Q2：换参数重跑，为什么结果没变？**
断点续传要求三个输出（summary.tsv、topologies.3line、tmr.gff3）都齐全才跳过。换输入或参数重跑前，先删掉这三个旧输出（或换输出目录）。

**Q3：报「conda 环境 python 不存在」？**
说明 conda 环境 `deeptmhmm_v.1.0` 没装或环境名不同。用 `--conda-env` 指定正确环境名，或在计算节点先初始化 conda。

**Q4：DeepTMHMM 目录在哪？**
默认 `~/software/deeptmhmm/DeepTMHMM-Academic-License-v1.0`（学术许可版本）。若装在不同位置，用 `--deeptmhmm-dir` 指定。