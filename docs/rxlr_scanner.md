# RxLR 效应蛋白扫描 | RxLR Effector Protein Scanner

一句话理解：**在蛋白质序列里，把 N 端窗口内带 RxLR(或 QxLR/GxLR)和 EER 基序的候选效应蛋白批量扫出来，输出 Excel/TSV 表格。**

## 功能概述 | Overview { #overview }

- 扫描 RxLR / QxLR / GxLR 与 EER 基序
- 默认窗口为第 21~120 位氨基酸(可调)
- 输出全量结果表 + 候选蛋白清单
- 纯 Python 实现，无外部生信软件依赖

## 快速开始 | Quick Start { #quick-start }

```bash
biopytools rxlr-scanner -i proteins.fa -o rxlr_results
```

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语 | 通俗理解 |
|------|----------|
| 效应蛋白(effector) | 病原分泌到宿主里的蛋白 |
| RxLR / EER | 卵菌效应子的标志性基序，像「身份标签」 |
| 基序(motif) | 序列上固定的小段特征 |
| 正则表达式 | 模糊匹配模式，如 `R.LR` 表示 R-任意氨基酸-L-R |
| 窗口(window) | 只看 N 端某一段，因为 RxLR 基序通常出现在信号肽切割点之后 |

## 输入 | Input { #input }

一个蛋白质 FASTA 文件：

```text
>protein1
MSKTRKVL...
>protein2
MAAQ...
```

## 参数说明 | Parameters { #parameters }

### 必需参数 | Required { #parameters-required }

**通俗理解|In plain words:** `-i` 给输入 FASTA；`-o` 是「输出文件前缀」，产物名会是 `<前缀>.xlsx`、`<前缀>.tsv` 等。

### 搜索窗口 | Search window { #parameters-window }

**通俗理解|In plain words:** 只在这段区间里找基序。`--window-start` 默认 20(对应第 21 位氨基酸)、`--window-end` 默认 120，正好覆盖经典 RxLR 基序出现的位置。`--min-length` 是最短序列要求，短于它的序列会被标记为「长度无效」。一般不用动。

### 输出 | Output { #parameters-output }

**通俗理解|In plain words:** `--output-dir` 才是输出目录(默认 ./rxlr_scanner_output)，`-o` 是文件名前缀。`--no-excel` / `--no-tsv` 分别关掉对应格式(至少要留一种)。

### 日志 | Logging { #parameters-log }

**通俗理解|In plain words:** `-v` 打印更详细日志；`--log-file` 额外写一份 `<前缀>.log` 到输出目录。

## 分析流程 | Pipeline { #pipeline }

```text
读入蛋白 FASTA
    │
    ▼
逐条序列取 N 端窗口(默认 21~120 位)
    │
    ▼
正则扫描 RxLR / QxLR / GxLR / EER 基序
    │
    ▼
判定候选(命中任一基序即为候选)
    │
    ▼
输出全量表 + 候选清单
```

## 输出 | Output { #output }

```text
rxlr_scanner_output/
├── rxlr_results.xlsx                    # 全量结果(Excel)
├── rxlr_results.tsv                     # 全量结果(TSV)
├── rxlr_results_candidates_only.tsv     # 仅候选蛋白(核心结果)
└── rxlr_results.log                     # 仅 --log-file 时生成
```

结果表关键列：`Sequence_ID`、`Sequence_Length`、`Valid_Length_(≥120)`、`Has_RxLR_motif`、`Has_EER_motif`、`RxLR_Candidate`、`Motif_Types`、`Motif_Positions`、`Total_Motif_Count`。

## 结果解读 | Interpreting Results { #interpreting-results }

- `RxLR_Candidate=Yes` 即候选效应蛋白，集中在 `*_candidates_only.tsv` 里
- `Has_RxLR_motif=Yes` 表示命中 RxLR/QxLR/GxLR 基序，`Has_EER_motif=Yes` 表示命中 EER 基序；两者都有的更典型
- `Motif_Positions` 给出基序在蛋白里的绝对位置(1-based)，可核对是否落在预期的 N 端区域
- `Valid_Length=No` 表示序列短于 `--min-length`，通常不纳入候选判断
- 这是「基序级」初筛，不能替代完整效应子鉴定(如需信号肽 + HMM 证据，见 phyto_effector 模块)

## 参数选择建议 | Parameter Guidance { #parameter-guidance }

- 常规扫描：默认窗口(21~120)+ 默认最短长度(120)即可
- 序列整体偏短(如片段化蛋白库)：可把 `--min-length` 调小(如 50)，避免把有效蛋白全标成无效
- 想放宽/收紧基序位置：调整 `--window-start` / `--window-end`
- 只需候选清单、不关心全量表：直接用 `*_candidates_only.tsv`

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 输入FASTA文件｜Input FASTA file path |
| `-o, --output-prefix` | 必填 |  | 输出文件前缀｜Output file prefix |
| `--window-start` | `20` | int | 窗口起始位置(20对应第21位)｜Window start position (20 for position 21) |
| `--window-end` | `120` | int | 窗口结束位置｜Window end position |
| `--min-length` | `120` | int | 最小序列长度｜Minimum sequence length |
| `--output-dir` | `./rxlr_scanner_output` | Path | 输出目录｜Output directory |
| `--no-excel` | — |  | 不生成Excel输出｜Do not generate Excel output |
| `--no-tsv` | — |  | 不生成TSV输出｜Do not generate TSV output |
| `-v, --verbose` | — |  | 详细输出模式｜Verbose output mode |
| `--log-file` | — |  | 生成日志文件｜Generate log file |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 输入FASTA文件｜Input FASTA file path |
| `-o, --output-prefix` | 必填 |  | 输出文件前缀｜Output file prefix |
| `--window-start` | `20` | int | 窗口起始位置(默认20对应第21位)｜Window start position (default: 20 for position 21) |
| `--window-end` | `120` | int | 窗口结束位置(默认120)｜Window end position (default: 120) |
| `--min-length` | `120` | int | 最小序列长度(默认120)｜Minimum sequence length (default: 120) |
| `--output-dir` | `./rxlr_scanner_output` |  | 输出目录｜Output directory (default: ./rxlr_scanner_output) |
| `--no-excel` | — | store_true | 不生成Excel输出｜Do not generate Excel output |
| `--no-tsv` | — | store_true | 不生成TSV输出｜Do not generate TSV output |
| `-v, --verbose` | — | store_true | 详细输出模式｜Verbose output mode |
| `--log-file` | — | store_true | 生成日志文件｜Generate log file |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies { #dependencies }

- 纯 Python 3 实现；输出 Excel 需要 pandas + openpyxl

## 常见问题 | FAQ { #faq }

**Q1：支持断点续传吗？**
不支持。单次扫描很快，重跑会覆盖同名输出。

**Q2：`-o` 和 `--output-dir` 什么区别？**
`-o` 是输出文件前缀(决定文件名)，`--output-dir` 才是输出目录(默认 ./rxlr_scanner_output)。

**Q3：序列太短会怎样？**
短于 `--min-length` 的序列会被标记为「长度无效」，但仍保留在结果表里，只是不参与候选判断。

**Q4：`--window-start 20` 为什么是第 21 位？**
内部用 0-based 坐标，20 对应第 21 个氨基酸。默认窗口即第 21~120 位，覆盖经典 RxLR 基序区。

**Q5：这和 phyto_effector 模块的 RxLR 鉴定有什么区别？**
本工具只做「基序扫描」初筛；phyto_effector 是完整鉴定流程(信号肽 + HMM + BLASTP + TMHMM 过滤)。需要严谨效应子清单时用 phyto_effector。