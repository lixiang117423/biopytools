# DeepBSA 批量 BSA 分析 | DeepBSA Batch BSA Analysis

一句话理解：**把 BSA(混池分离分析)里常用的 7 种找「性状关联区域」的方法一次性并行跑完，并把结果自动合并成一张总表和一套图**——输入一个 BSA 的 VCF(高表型池 vs 低表型池)，输出每种方法找到的候选区域。

## 功能概述 | Overview

- 封装 DeepBSA 工具，批量运行 7 种 BSA 方法：`DL`、`K`、`ED4`、`SNP`、`SmoothG`、`SmoothLOD`、`Ridit`
- 方法级并行执行，每个方法独立工作目录，互不冲突
- 自动结果合并：CSV 结果(标注来源方法)+ PNG 图片(加方法名前缀)+ 汇总报告
- 自动处理复杂 VCF 文件名、为 DL 方法自动创建 Models 符号链接
- 支持断点续传(某方法已有 `values.txt` 则跳过)，`--force` 强制重跑

## 快速开始 | Quick Start

```bash
biopytools deepbsa run -i variant.vcf -o bsa_results
```

最小输入：一个 BSA 的 VCF，默认并行跑全部 7 种方法并自动合并结果。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| BSA(混池分离分析) | 把一群「表现极端」的个体(如最高的几十株、最矮的几十株)分别混在一起测序，找两组之间差异最大的基因组区域 |
| 混池(pool) | 把多个个体的 DNA 混在一起测，省成本、看「整体倾向」 |
| QTL | 控制性状的基因组区域，BSA 要找的东西 |
| SNP-index | 每个位点上「来自高表型亲本的读数占比」，两组差得越大越像候选区 |
| 候选区域 | 曲线上冒尖的峰，通常是一段几 Mb 的区间 |
| DL / K / ED4 / SNP / SmoothG / SmoothLOD / Ridit | 7 种不同的统计/机器学习算法，各自从不同角度算「两组差多大」，互相印证 |

## 输入 | Input

- `run` 的输入是一个 BSA 的 VCF(或 CSV)。VCF 通常来自「高表型池 + 低表型池」的变异检测，含 AD(等位深度)信息最佳。
- 复杂文件名(如 `variation.filtered.snp.biallelic.vcf`)会自动处理(内部建符号链接)，无需改名。
- `vcf2csv` 子命令可把 VCF 转成 DeepBSA 用的 CSV(提取 FORMAT 里的 AD 信息)：`biopytools deepbsa vcf2csv -i input.vcf -o out.csv`。

## 参数说明 | Parameters

### run - 运行分析 | Run

**通俗理解|In plain words:** 最常用的子命令。`-m` 选跑哪些方法(默认全部 7 个，逗号分隔)；默认**并行**跑，`--no-parallel` 改串行；`--threads` 是每个方法的线程数(默认 6，并行时别设太大)；`--force` 强制重跑已完成的步骤。

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `-i, --input-file` | 必填 | 输入 VCF/CSV |
| `-o, --output-dir` | `deepbsa_results` | 输出目录 |
| `-m, --methods` | 全部 | 方法列表，逗号分隔(如 DL,K,ED4) |
| `-p/--parallel` | 开 | 并行运行(默认) |
| `--no-parallel` | 关 | 串行运行 |
| `--threads` | `6` | 每个方法线程数 |
| `--smooth-func` | `Tri-kernel-smooth` | 平滑函数 |
| `--skip-merge` | 关 | 跳过结果合并 |
| `-f, --force` | 关 | 强制重跑 |

### batch - 生成批量脚本 | Batch

**通俗理解|In plain words:** 生成一个 shell 脚本(和按方法拆分的作业脚本)，适合在集群上投递。生成后手动运行 `run.sh` 再合并。

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `-i, --input-file` | 必填 | 输入 VCF/CSV |
| `-o, --output-dir` | 必填 | 批量任务输出目录 |
| `-s, --script-name` | `run_deepbsa_methods.sh` | 生成的脚本名 |
| `--threads` | `88` | 每个方法线程数 |

### merge - 合并结果 | Merge

**通俗理解|In plain words:** 单独把某个输出目录的结果重新合并(高级用法，一般 run 已自动合并)。

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `-i, --input-dir` | 必填 | DeepBSA 输出目录 |
| `-o, --output-dir` | 必填 | 合并结果目录 |
| `-m, --methods` | 全部 | 要合并的方法 |

### vcf2csv - VCF 转 CSV | VCF to CSV

**通俗理解|In plain words:** 从 VCF 提取 AD(等位深度)转成 DeepBSA 认识的 CSV，做输入预处理。

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `-i, --input-file` | 必填 | 输入 VCF |
| `-o, --output-file` | 必填 | 输出 CSV |

## 分析流程 | Pipeline

**通俗理解|In plain words:** 先给输入(必要时 vcf2csv 转格式)→ 为每个方法建独立目录(并行跑)→ 各自出 QTL 结果和图 → 汇总成总表和总图。

```text
输入 VCF/CSV
    │ (可选 vcf2csv 转格式)
    ▼
为每个方法建独立目录(复杂文件名自动建符号链接;DL 自动建 Models 链接)
    │
    ▼
并行跑 7 种方法(各自独立日志)
    │
    ▼
自动合并: merged_results.csv + images/ + summary_report.txt
```

## 输出 | Output

```text
bsa_results/
├── deepbsa.log                    # 主日志
├── merged_results/                # 合并结果(各方法结果汇总)
│   ├── merged_results.csv         # 所有方法的 QTL 结果(含 Method 列)
│   ├── images/                    # 所有方法图片(带方法名前缀)
│   └── summary_report.txt         # 各方法统计汇总
├── DL/                            # DL 方法目录(Models 符号链接, Results/, DL.log)
├── K/  ├── ED4/  ├── SNP/  ├── SmoothG/  ├── SmoothLOD/  └── Ridit/
```

## 结果解读 | Interpreting Results

### 合并总表(`merged_results/merged_results.csv`)

**通俗理解|In plain words:** 一张表看全部 7 种方法找到的候选区。含 `Method`(来源方法)、`Chr/Left/Peak/Right`(候选区位置)、`Value`(方法打分)。**多种方法在同一个区域都出峰，这个区域最可信。**

### 各方法图(`merged_results/images/`)

**通俗理解|In plain words:** 每种方法一张曲线图，看「冒尖的峰」就是候选区域。

### 汇总报告(`summary_report.txt`)

**通俗理解|In plain words:** 每方法出了几个 QTL、最高分多少，一眼看各方法表现。

## 参数选择建议 | Parameter Guidance

| 场景 | 建议 |
|------|------|
| 常规分析 | `run` 全默认(并行 7 方法) |
| 只想跑某几个方法 | `-m ED4,SNP` |
| 集群投递 | 用 `batch` 生成脚本再投递 |
| 内存/核数有限 | `--no-parallel` 串行，`--threads 16` |
| 重跑某方法 | 删掉对应方法目录，或 `-f` 强制重跑 |
| VCF 含 AD 想转 CSV | 先 `vcf2csv` 预处理 |

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input-file` | 必填 |  | 输入文件路径（VCF/CSV格式）｜Input file path (VCF/CSV format) |
| `-m, --methods` | — |  | 要运行的方法，逗号分隔（默认：全部）｜Methods to run, comma-separated (default: all) |
| `-o, --output-dir` | `deepbsa_results` |  | 输出目录（默认：deepbsa_results）｜Output directory (default: deepbsa_results) |
| `--no-parallel` | — | store_false | 串行运行所有方法｜Run all methods serially |
| `-p, --parallel` | — | store_true | 并行运行所有方法（默认）｜Run all methods in parallel (default) |
| `-n, --no-auto-clean` | — | store_true | 不自动清理VCF注释行｜Don't auto clean VCF comment lines |
| `-k, --keep-clean` | — | store_true | 保留清理后的文件｜Keep cleaned file |
| `--deepbsa-path` | `~/software/DeepBSA/DeepBSA_multithread/bin/main.py` |  | DeepBSA主程序路径（多线程版本）｜DeepBSA main script path (multithread version) |
| `--conda-env` | `~/miniforge3/envs/DeepBSA` |  | Conda环境路径｜Conda environment path |
| `--threads` | `0` | int | DeepBSA线程数（0=自动检测，1=单线程，>1=多线程，默认：0）｜Number of threads for DeepBSA (0=auto, 1=single, >1=multi, default: 0) |
| `--smooth-func` | `Tri-kernel-smooth` | Tri-kernel-smooth/LOWESS/Moving Average | 平滑函数（默认：Tri-kernel-smooth）｜Smooth function (default: Tri-kernel-smooth) |
| `-v, --verbose` | — | store_true | 详细输出｜Verbose output |
| `-f, --force` | — | store_true | 强制重新运行所有步骤（不跳过已完成）｜Force rerun all steps (don't skip completed) |
| `--skip-merge` | — | store_true | 跳过合并步骤｜Skip merge step |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- DeepBSA 主程序：默认用**内置版本**(`deepbsa_builtin`)，也可 `--deepbsa-path` 指定
- conda 环境 `~/miniforge3/envs/DeepBSA`(DeepBSA 运行所需 Python 依赖)
- DL 方法需预训练模型(内置 Models 自动链接)

## 常见问题 | FAQ

**Q1：`biopytools deepbsa` 和 `python -m biopytools.deepbsa` 有什么区别？**
`biopytools deepbsa` 走子命令接口(batch/run/merge/vcf2csv，默认 `run --threads 6`、内置 DeepBSA)；`python -m biopytools.deepbsa` 是旧版扁平接口(见下方参数表，默认 `--threads 0`、外置 multithread 路径)。**日常用 `biopytools deepbsa` 即可**。

**Q2：DL 方法报错找不到 Models？**
本工具已自动为 DL 创建 Models 符号链接；若仍失败，检查内置 `deepbsa_builtin/Models` 是否存在。

**Q3：某方法失败怎么办？**
看 `输出目录/METHOD_NAME/METHOD_NAME.log`；常见原因是 VCF 格式不对或 SNP 太少。可 `-m` 单独重跑该方法。

**Q4：换参数重跑没生效？**
断点续传按 `values.txt` 存在性跳过。换参数后用 `-f` 强制重跑，或删除对应方法目录。

**Q5：能并行吗？内存够吗？**
默认并行 7 方法，优化目标是 64 核 + 300GB 集群；核数/内存有限时用 `--no-parallel` 串行，或 `--threads` 调小。
