# Phobius - 跨膜拓扑与信号肽预测 | Phobius Transmembrane Topology & Signal Peptide Prediction

一句话理解：**给每条蛋白判断「有没有信号肽、有几个跨膜区、跨膜区在哪」**。输入一个蛋白质 FASTA，输出一张汇总表，用于区分分泌蛋白和膜蛋白、标注结构位置。

## 功能概述 | Overview

- 同时预测信号肽（signal peptide）与跨膜螺旋（TM helix），二者联合建模、互相校正
- 运行 Phobius 的 `-short` 与 `-long` 两种输出，合并解析成一张干净 TSV
- 输出每个蛋白的跨膜区个数与坐标、信号肽有无与位置、以及 Phobius 的预测类型
- 断点续传：short/long 两个原始输出各自独立判断，已存在则跳过对应步骤

## 快速开始 | Quick Start

```bash
biopytools phobius -i proteins.fa -o output_dir/
```

最小输入：一个蛋白质 FASTA 文件 + 输出目录。输出前缀默认取输入文件名。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| 信号肽 | 蛋白质 N 端的一小段「地址标签」，引导蛋白分泌或定位，通常随后被切掉 |
| 跨膜区 | 蛋白质「横穿」细胞膜的一段，像一根钉子钉穿木板 |
| 拓扑(topology) | 蛋白质相对细胞膜「怎么穿」的整体结构：哪段在膜内、哪段在膜外 |
| -short / -long | Phobius 的两种输出：short 给一句话结论，long 给逐特征的坐标表 |

## 输入 | Input

蛋白质 FASTA 文件（支持 .fa / .faa / .fasta / .fna / .ffn，含 .gz 压缩）：

```text
>protein_01
MKKLLIAAMMAAALAACSQEAKTEVFSKSADEGGAPK...
```

- 序列必须是氨基酸（蛋白质）
- 输入文件名自动作为默认输出前缀（剥离扩展名，空格/斜杠替换为下划线）

## 参数说明 | Parameters

### 必需参数 | Required

**通俗理解|In plain words:** 两个必填：输入蛋白文件、输出目录。没有默认值。

### 输出前缀与路径 | Prefix & path

**通俗理解|In plain words:** `--prefix` 决定输出文件名，默认取输入文件名，**一般不用动**；`--phobius-path` 指向 phobius.pl 的安装位置，部署时已配好，**普通用户不用动**。

## 分析流程 | Pipeline

**通俗理解|In plain words:** 分别跑 short（一句话结论）和 long（坐标特征表），各自可断点续传，再合并解析成汇总表。

```text
蛋白质 FASTA
    │
    ├─ 步骤1: phobius -short → {prefix}_phobius_short.txt  (已存在则跳过)
    ├─ 步骤2: phobius -long  → {prefix}_phobius_long.txt   (已存在则跳过)
    │
    ▼
解析 short(结论) + long(坐标) → 合并 → {prefix}_phobius.tsv
```

## 输出 | Output

```text
output_dir/
├── {prefix}_phobius_short.txt    # Phobius -short 原始输出
├── {prefix}_phobius_long.txt     # Phobius -long 原始输出(FT 特征表)
├── {prefix}_phobius.tsv          # 主结果:合并汇总表
└── 99_logs/
    └── phobius.log               # 运行日志
```

- `{prefix}` 默认是输入文件名
- 主结果 `{prefix}_phobius.tsv` 列：`ID / TM / SP / SP_region / TM_regions / Prediction`

## 结果解读 | Interpreting Results

### 1. 汇总表（`{prefix}_phobius.tsv`）

**通俗理解|In plain words:** 每行一个蛋白，看 `SP` 和 `TM` 两列即可定性。

- `TM`：跨膜区个数。0=非跨膜蛋白；≥1=跨膜（膜）蛋白
- `SP`：有无信号肽（Y/N）；`SP_region`：信号肽坐标区间
- `TM_regions`：各跨膜区坐标，多个用 `;` 分隔，无则 `-`
- `Prediction`：Phobius short 给出的一句话预测类型（如 SIGNAL、TM、CYT 等）

### 2. 好坏判据 | Judgment

- **有信号肽 + 无跨膜区**：分泌蛋白的典型特征
- **有信号肽 + 1 个跨膜区**：单次跨膜受体/膜锚定蛋白的常见模式；**跨膜区 ≥2**：多次跨膜蛋白
- Phobius 把信号肽和跨膜区联合建模，比单纯判断「N 端疏水」更准，尤其能区分二者

## 参数选择建议 | Parameter Guidance

- **`--prefix`**：同输入在多个项目复用想区分结果时设置，否则用默认
- **`--phobius-path`**：仅部署环境不同时才动
- 本工具参数极少，绝大多数场景一条最简命令即可

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input, -i` | 必填 |  | 输入蛋白质FASTA｜Input protein FASTA |
| `--output-dir, -o` | 必填 | Path | 输出目录｜Output directory |
| `--prefix` | — |  | 输出前缀(默认输入文件名)｜Output prefix (default: input filename) |
| `--phobius-path` | — |  | phobius.pl路径｜phobius.pl path |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | [FILE] 输入蛋白质FASTA｜Input protein FASTA |
| `-o, --output-dir` | 必填 |  | [DIR] 输出目录｜Output directory |
| `--prefix` | — |  | [STR] 输出前缀(默认输入文件名)｜Output prefix (default: input filename) |
| `--phobius-path` | `~/miniforge3/envs/protein/bin/phobius.pl` |  | phobius.pl路径｜phobius.pl path |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- Phobius（`phobius.pl`，默认 `~/miniforge3/envs/protein/bin/phobius.pl`，通过 conda 环境自动包装调用）
- Python 3

## 常见问题 | FAQ

**Q1：short 和 long 各跑一遍，能不能只跑一个？**
不能跳过——本工具用 short 拿「结论」、long 拿「坐标」，二者合并才得到完整 TSV。断点续传会分别检查这两个原始文件，缺哪个补跑哪个。

**Q2：换参数重跑，为什么结果没变？**
断点续传按 `_phobius_short.txt` / `_phobius_long.txt` 是否存在判断。想强制重跑，删掉这两个文件（或换输出目录）。

**Q3：报「phobius.pl 不存在」？**
默认路径 `~/miniforge3/envs/protein/bin/phobius.pl` 里没有该脚本。用 `--phobius-path` 指定实际位置，或确认 protein 环境已装 Phobius。