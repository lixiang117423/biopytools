# SignalP - 信号肽预测 | SignalP Signal Peptide Prediction

一句话理解：**判断每条蛋白「带不带信号肽、是不是会被分泌出去」，并标出信号肽在哪、切割位点在哪**。输入一个蛋白质 FASTA，输出预测结果，用于识别分泌蛋白、脂蛋白等。

## 功能概述 | Overview

- 基于 SignalP 6.0 深度学习模型预测信号肽及切割位点
- 支持真核（默认）与细菌（other）两种生物类型；输出 txt / png / eps / all / none 多种格式
- 三种速度模式：fast（默认）、slow、slow-sequential，兼顾速度与精度
- 自动把原始 `prediction_results.txt` 转成中文易读版与 R 友好纯 TSV（`signalp_summary.tsv`）
- 默认自动清理 plot 文件省空间（`--keep-plots` 保留）
- 断点续传：主输出 `prediction_results.txt` 已存在则跳过 signalp 运行，直接重做格式化

## 快速开始 | Quick Start

```bash
biopytools signalp -i proteins.faa -o output_dir
```

最小输入：一个蛋白质 FASTA 文件 + 输出目录。默认真核生物、fast 模式、txt 格式。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| 信号肽 | 蛋白质 N 端的一小段「地址标签」，引导蛋白走分泌通路，随后被信号肽酶切掉 |
| 切割位点(CS) | 信号肽被剪断的确切位置，常写成「20-21」表示第 20/21 位之间 |
| SP 概率 | 模型判断「这是信号肽」的把握，0–1，越接近 1 越像有信号肽 |
| Sec/SPI、Tat/SPI | 经典分泌途径；SPII 脂蛋白、SPIII 外膜蛋白（细菌特有） |
| organism | 生物类型：eukarya=真核（动植物真菌），other=原核/细菌 |

## 输入 | Input

蛋白质 FASTA 文件（氨基酸序列，.fa / .fasta / .faa / .ffn / .fna）：

```text
>protein_01
MKKLLIAAMMAAALAACSQEAKTEVFSKSADEGGAPK...
```

- 序列必须是氨基酸（蛋白质），不接受核酸
- 扩展名需为 .fa/.fasta/.faa/.ffn/.fna 之一（否则校验报错）

## 参数说明 | Parameters

### 必需参数 | Required

**通俗理解|In plain words:** 两个必填：输入蛋白文件、输出目录。没有默认值。

### 生物类型与模式 | Organism & mode

**通俗理解|In plain words:** `--organism` 告诉模型你的序列是哪种生物（真核还是细菌），选错会明显影响准确率。`--mode` 选速度：fast 最快（默认，够用），slow 更慢更准，slow-sequential 单线程最慢。**日常用 fast 即可。**

### 输出格式与性能 | Format & performance

**通俗理解|In plain words:** `--format` 决定生成哪些文件（默认 txt 文字结果）；`--bsize`/`--write-procs`/`--torch-num-threads` 是批次大小和并行线程数，**默认 12 一般不用动**，超大蛋白组或内存紧张时才调。

### 运行环境 | Runtime

**通俗理解|In plain words:** `--signalp-path`/`--model-dir` 指向程序与模型权重，部署时已配好，**普通用户不用动**；`--skip-resolve` 跳过结果冲突处理，`--keep-plots` 保留 plot 文件，日常都不用开。

## 分析流程 | Pipeline

**通俗理解|In plain words:** 先查主输出是否已存在（有就只重做格式化），否则运行 signalp6，再转出易读版和 R 友好 TSV，最后清理 plot 文件。

```text
蛋白质 FASTA
    │
    ▼
(断点续传) prediction_results.txt 已存在? → 是则跳过 signalp 运行
    │
    ▼
运行 signalp6 → prediction_results.txt + output.gff3 + output.json 等
    │
    ▼
格式化: 中文易读 prediction_results_readable.txt + R友好 signalp_summary.tsv
    │
    ▼
(默认) 清理 output_*_plot.txt + 记录软件版本
```

## 输出 | Output

```text
output_dir/
├── 00_pipeline_info/
│   └── software_versions.yml              # 软件版本与参数记录
├── prediction_results.txt                 # SignalP 原始预测(主输出,续传判据)
├── prediction_results_readable.txt        # 中文易读版(含统计与图例)
├── signalp_summary.tsv                    # R友好纯TSV(英文列名,无注释行)
├── processed_entries.fasta                # 处理后的序列
├── output.gff3                            # 预测的 GFF3 注释
├── region_output.gff3                     # 区域 GFF3 注释
├── output.json                            # 预测的 JSON
└── 99_logs/
    └── signalp.log                        # 运行日志
```

- `signalp_summary.tsv` 是给下游分析（R/Excel）用的主表：`蛋白质ID / 预测代码 / Has_SP / SP概率 / SP起止 / 切割位点起止 / 切割置信度`
- `output_*_plot.txt` 是逐蛋白绘图数据，默认被清理，用 `--keep-plots` 可保留

## 结果解读 | Interpreting Results

### 1. 预测类型（`signalp_summary.tsv` 的「预测代码」列）

**通俗理解|In plain words:** 每个蛋白会被归到下面某一类，看代码即可定性。

| 代码 | 含义 |
|------|------|
| SP / SPI | 分泌蛋白-经典信号肽(Sec/SPI) |
| SPII | 脂蛋白-脂蛋白信号肽(细菌) |
| SPIII | 外膜蛋白-细菌信号肽(细菌) |
| Tat | 双精氨酸转运蛋白-Tat信号肽(细菌) |
| OTHER | 非信号肽蛋白 |

### 2. 关键列 | Key columns

- `Has_SP`：TRUE/FALSE，是否有信号肽（OTHER 为 FALSE）
- `SP概率`：模型给出的信号肽概率（0–1），越高越可信
- `切割位点起止`：预测的信号肽切割位置；`切割置信度` 是该切割的把握

### 3. 好坏判据 | Judgment

- **SP 概率接近 1 且预测为 SP/SPI**：高度可信的分泌蛋白；概率 <0.5 的 SP 判读要谨慎
- **OTHER**：不带信号肽（胞内蛋白或经其他途径定位的蛋白）
- 预测类型为 SPII/SPIII/Tat 只对细菌有意义，真核数据出现这类应检查 `--organism` 是否设错

## 参数选择建议 | Parameter Guidance

- **`--organism`**：真核（动植物真菌）用默认 eukarya；细菌/古菌改用 other，选错会影响准确率
- **`--mode`**：批量初筛用 fast；对关键候选要最准结论时用 slow
- **`--format txt`**（默认）够用；需要可视化图时用 `png` 或 `all`
- **`--keep-plots`**：默认清理 plot 省空间，只有需要逐蛋白画图数据时才加

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input, -i` | 必填 |  | 输入FASTA文件(氨基酸序列)｜Input FASTA file (amino acid sequences) |
| `--output-dir, -o` | 必填 | Path | 输出目录｜Output directory |
| `--organism, -org` | `eukarya` | eukarya/other/euk | 生物类型｜Organism type (default: eukarya) |
| `--format, -fmt` | `txt` | txt/png/eps/all/none | 输出格式｜Output format (default: txt) |
| `--mode, -m` | `fast` | fast/slow/slow-sequential | 预测模式｜Prediction mode (default: fast) |
| `--bsize, -bs` | `12` | int | 批处理大小｜Batch size |
| `--write-procs, -wp` | `12` | int | 写入进程数｜Number of write processes |
| `--torch-num-threads, -tt` | `12` | int | PyTorch线程数｜PyTorch threads |
| `--signalp-path` | `~/miniforge3/envs/protein/bin/signalp6` |  | SignalP程序路径｜SignalP program path |
| `--model-dir, -md` | — |  | 模型权重目录｜Model weights directory |
| `--skip-resolve` | — |  | 跳过冲突解析｜Skip conflict resolution |
| `--keep-plots` | — |  | 保留plot文件(默认自动清理)｜Keep plot files (auto-deleted by default) |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 输入FASTA文件(氨基酸序列)｜Input FASTA file (amino acid sequences) |
| `-o, --output-dir` | 必填 |  | 输出目录｜Output directory |
| `--organism` | `eukarya` | eukarya/other/euk | 生物类型｜Organism type (default: eukarya) |
| `--format` | `txt` | txt/png/eps/all/none | 输出格式｜Output format (default: txt) |
| `--mode` | `fast` | fast/slow/slow-sequential | 预测模式｜Prediction mode (default: fast) |
| `--bsize` | `12` | int | 批处理大小｜Batch size (default: 12) |
| `--write-procs` | `12` | int | 写入进程数｜Number of write processes (default: 12) |
| `--torch-num-threads` | `12` | int | PyTorch线程数｜PyTorch threads (default: 12) |
| `--signalp-path` | `~/miniforge3/envs/protein/bin/signalp6` |  | SignalP程序路径｜SignalP program path |
| `--model-dir` | — |  | 模型权重目录｜Model weights directory |
| `--skip-resolve` | — | store_true | 跳过resolve步骤｜Skip resolve step |
| `--keep-plots` | — | store_true | 保留plot文件(默认自动清理plot文件)｜Keep plot files (plots are auto-deleted by default) |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- SignalP 6.0（`signalp6`，默认 `~/miniforge3/envs/protein/bin/signalp6`，通过 conda 环境自动包装调用）
- SignalP 模型权重（`--model-dir` 可指定）
- PyTorch（SignalP 6.0 为深度学习模型，运行时需要）
- Python 3 + PyYAML（记录版本）

## 常见问题 | FAQ

**Q1：换 `--organism` 或 `--mode` 重跑，为什么结果没变？**
断点续传按 `prediction_results.txt` 是否存在判断。换参数重跑前，先删掉输出目录里的 `prediction_results.txt`（或换输出目录），否则会复用旧预测只重做格式化。

**Q2：plot 文件去哪了？**
默认在格式化后自动删除 `output_*_plot.txt` 以省空间。需要时加 `--keep-plots` 保留。

**Q3：`--organism` 里 euk 和 eukarya 有什么区别？**
没区别，euk 是 eukarya 的别名，程序会自动归一到 eukarya。

**Q4：信号肽概率高但预测 OTHER？**
通常不会——概率与类型是一致的。若出现异常，检查输入是否有非标准序列，或换 slow 模式复核。