# Primer3 批量引物设计 | Primer3 Batch Primer Design

一句话：给它一批 DNA 序列（FASTA 文件），它自动调用 Primer3 为每条序列设计一对（或多对）PCR 引物，并整理成一张带 Tm、GC 含量、二聚体评估的结果表格，可直接用于合成订单或实验记录。

## 功能概述 | Overview

- 批量输入 FASTA，单次运行完成全部引物设计
- 两种设计策略：`all`（默认，引物锚定在序列两端，适合扩增完整插入片段）与 `random`（在序列内部自由选位）
- 引物长度、退火温度（Tm）、产物大小、引物对数量等参数全部可调
- 结果输出 CSV/TSV/XLSX，中英文表头可选，含二聚体与发夹结构热力学评估值
- 支持断点续传：primer3 原始输出已存在时自动跳过重新计算
- 自动记录软件版本到 `00_pipeline_info/software_versions.yml`

## 快速开始 | Quick Start

```bash
biopytools primer3 -i sequences.fasta -o primer3_output
```

## 零基础概念速览 | Concepts in plain words

| 术语<br>Term | 通俗解释<br>Plain explanation |
|---|---|
| 引物 primer | 一小段（约 20 个碱基）人工合成的 DNA"起跑线"。DNA 聚合酶只能从这段起跑线开始复制，一对引物（一前一后）夹住的区段就是 PCR 会扩增的目标 |
| Tm（解链温度） | 引物和模板"粘住"的温度门槛，可以理解为双链DNA"融化"一半时的温度。两条引物的 Tm 要接近，PCR 才能同时工作。一般设在 58-62°C |
| GC 含量 | 引物里 G+C 碱基的比例。GC 之间有 3 个氢键，粘得更牢。40-60% 比较合适，太高太低都会降低扩增效率 |
| GC clamp | 引物 3' 末端（起跑线的"踩踏端"）放 1-2 个 G/C，像锚一样抓牢模板，提高扩增特异性 |
| 二聚体 dimer | 两条引物互相粘在一起而不是粘模板，PCR 就白做了。工具给出量化分数（TH 值），分数越低越好 |
| 产物大小 product size | 一对引物夹住的片段长度。常规检测 100-300 bp；要测完整基因或做测序则更长 |
| Penalty（惩罚值） | Primer3 给每对引物的综合打分，综合了长度、Tm、GC、二聚体等所有指标。越低越好，低于 1 通常令人满意 |

## 输入 | Input

一个 FASTA 文件（`-i`），每条序列为一行标题（`>` 开头）加若干行序列：

```text
>gene_alpha
GCTAAAGACAATTACATAACATACACGTCAGC...
>gene_beta
TTTGGGGCCCCAAAATTTTGGGGCCCCAAAA...
```

要求与说明：

- 序列只含 A/C/G/T/N（其他字符会被 Primer3 拒绝）
- 序列 ID 中的空格、`|`、`:` 等特殊字符会自动替换为下划线（Primer3 对 ID 有格式限制），结果表中以替换后的 ID 呈现
- `method=all` 模式下序列长度应大于约 150 bp（两端各留出引物结合位点的空间），过短序列会设计失败

## 参数说明 | Parameters

### 输入输出 | Input/Output

**通俗理解|In plain words:** 这组管"文件从哪来、结果到哪去"。一般只需要动 `-i` 和 `-o`；表头语言按实验记录本习惯选，给国内同事看选 `zh`，写英文论文附录选 `en`。

### 引物与产物 | Primer and product

**通俗理解|In plain words:** 这组管"引物长什么样、扩增出多长的片段"。默认值是常规检测场景的成熟起点，一般不用动。引物长度加长（如 22-25）特异性更好但合成更贵；Tm 范围跟着你实验室的 PCR 退火温度走（退火温度通常设为 Tm 下限减 3-5°C）；产物大小在 `random` 模式下生效，`all` 模式会被自动产物范围覆盖。

### 设计策略 | Design strategy

**通俗理解|In plain words:** 这组管"引物允许站在序列的哪些位置"。`all`（默认）把两条引物分别钉在序列头尾 `--primer-end-margin` bp 的范围内，适合"我要扩增这整段序列"的场景；`random` 允许引物落在序列任何位置，适合"只要在这个序列里扩增出一段就行"的场景。`--auto-product-size` 在 `all` 模式下让产物覆盖整条序列，无需手算。

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input-fasta, -i` | 必填 |  | 输入FASTA文件｜Input FASTA file path |
| `--output-dir, -o` | 必填 | Path | 输出目录｜Output directory |
| `--primer3-core-path` | `~/miniforge3/envs/misc/bin/primer3_core` |  | Primer3核心程序路径｜Primer3 core program path |
| `--primer-min-size` | `20` | int | 最小引物长度｜Minimum primer size |
| `--primer-opt-size` | `20` | int | 最优引物长度｜Optimal primer size |
| `--primer-max-size` | `22` | int | 最大引物长度｜Maximum primer size |
| `--primer-min-tm` | `53.0` | float | 最小退火温度(°C)｜Minimum annealing temperature (°C) |
| `--primer-opt-tm` | `58.0` | float | 最优退火温度(°C)｜Optimal annealing temperature (°C) |
| `--primer-max-tm` | `63.0` | float | 最大退火温度(°C)｜Maximum annealing temperature (°C) |
| `--product-min-size` | `100` | int | 最小产物大小(bp)｜Minimum product size (bp) |
| `--product-max-size` | `300` | int | 最大产物大小(bp)｜Maximum product size (bp) |
| `--primer-num-return` | `5` | int | 返回引物对数量｜Number of primer pairs to return |
| `--primer-max-ns` | `0` | int | 允许的N碱基数量｜Number of N bases accepted |
| `--primer-gc-clamp` | `1` | int | GC clamp数量｜GC clamp count |
| `--output-format` | `csv` | csv/tsv/xlsx | 输出文件格式｜Output file format |
| `--output-header-lang` | `zh` | zh/en | 输出表头语言(zh:中文, en:英文)｜Output header language (zh: Chinese, en: English) |
| `--method, -m` | `all` | all/random | 引物设计策略: all=覆盖头尾, random=随机设计｜Primer design strategy: all=cover ends, random=random design |
| `--primer-end-margin` | `200` | int | 两端允许的引物位置范围bp(仅用于method=all)｜Allowed margin at ends in bp (only for method=all) |
| `--auto-product-size/--no-auto-product-size` | `True` |  | 自动根据序列长度设置产物大小范围(默认开启)｜Auto set product size range based on sequence length (enabled by default) |
| `--product-size-min-ratio` | `0.5` | float | 产物最小长度占序列长度的比例｜Min product size ratio to sequence length (default: 0.5) |
| `--product-size-max-ratio` | `1.0` | float | 产物最大长度占序列长度的比例｜Max product size ratio to sequence length (default: 1.0) |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input-fasta` | 必填 |  | 输入FASTA文件路径｜Input FASTA file path |
| `-o, --output-dir` | 必填 |  | 输出目录路径｜Output directory path |
| `--primer3-core-path` | `~/miniforge3/envs/misc/bin/primer3_core` |  | Primer3核心程序路径｜Primer3 core program path |
| `--primer-min-size` | `20` | int | 最小引物长度｜Minimum primer size |
| `--primer-opt-size` | `20` | int | 最优引物长度｜Optimal primer size |
| `--primer-max-size` | `22` | int | 最大引物长度｜Maximum primer size |
| `--primer-min-tm` | `53.0` | float | 最小退火温度(°C)｜Minimum annealing temperature (°C) |
| `--primer-opt-tm` | `58.0` | float | 最优退火温度(°C)｜Optimal annealing temperature (°C) |
| `--primer-max-tm` | `63.0` | float | 最大退火温度(°C)｜Maximum annealing temperature (°C) |
| `--product-min-size` | `100` | int | 最小产物大小(bp)｜Minimum product size (bp) |
| `--product-max-size` | `300` | int | 最大产物大小(bp)｜Maximum product size (bp) |
| `--primer-num-return` | `5` | int | 返回引物对数量｜Number of primer pairs to return |
| `--primer-max-ns` | `0` | int | 允许的N碱基数量｜Number of N bases accepted |
| `--primer-gc-clamp` | `1` | int | GC clamp数量｜GC clamp count |
| `--output-format` | `csv` | csv/tsv/xlsx | 输出文件格式｜Output file format |
| `--output-header-lang` | `zh` | zh/en | 输出表头语言(zh:中文, en:英文)｜Output header language (zh: Chinese, en: English) |
| `--method, -m` | `all` | all/random | 引物设计策略: all=覆盖头尾(默认), random=随机设计｜Primer design strategy: all=cover ends (default), random=random design |
| `--primer-end-margin` | `200` | int | 两端允许的引物位置范围bp,仅用于method=all｜Allowed margin at ends in bp (only for method=all) |
| `--auto-product-size` | `True` | store_true | 自动根据序列长度设置产物大小范围(默认开启)｜Auto set product size range based on sequence length (enabled by default) |
| `--no-auto-product-size` | — | store_false | 禁用自动产物大小范围｜Disable automatic product size range |
| `--product-size-min-ratio` | `0.5` | float | 产物最小长度占序列长度的比例｜Min product size ratio to sequence length (default: 0.5) |
| `--product-size-max-ratio` | `1.0` | float | 产物最大长度占序列长度的比例｜Max product size ratio to sequence length (default: 1.0) |

<!-- END PARAMS:auto -->

## 输出 | Output

```text
primer3_output/
├── 00_pipeline_info/
│   └── software_versions.yml      # primer3 版本与本次关键参数
├── 01_primer_design/
│   ├── primer3_output.txt         # Primer3 原始输出(断点续传依据)
│   └── primers_result.csv         # 设计结果总表(--output-format 决定后缀)
├── 99_logs/
│   └── primer3_design.log         # 运行日志(含全部执行命令)
└── tmp/                           # 临时目录(运行结束自动清空)
```

- **`primers_result.csv`**：核心结果，每行一对引物，含正/反向引物序列、长度、Tm、GC 含量、退火温度（取两条引物 Tm 的较小值）、产物大小、惩罚值与各项二聚体评分
- **`primer3_output.txt`**：Primer3 原始输出。重跑时若此文件存在会跳过重新计算（断点续传）；它也是排查"某序列为什么没有引物"的第一手材料
- **`software_versions.yml`**：论文 Methods 与问题复现时使用

## 结果解读 | Interpreting Results

结果表逐列怎么读、好坏判据：

| 列<br>Column | 怎么读<br>How to read |
|---|---|
| 惩罚值 Penalty | 越低越好。`<1` 可放心用；`1-5` 可用但建议对照其他指标；`>10` 建议放弃或调整参数重新设计 |
| 退火温度 Annealing_Temp | 两条引物 Tm 的较小值。PCR 退火温度建议设为该值减 3-5°C；同一反应中多条引物对的退火温度差最好不超过 5°C |
| 正/反向引物 Tm 差 | 同一对引物两条 Tm 相差 `<5°C` 理想；`>10°C` 扩增可能偏斜 |
| 正/反向 GC 含量 | `40-60%` 合适；`>70%` 或 `<30%` 扩增效率通常下降 |
| 产物大小 | 按实验目的判断：酶切鉴定 `200-500 bp`、qPCR `70-200 bp`、测序 `400-800 bp` |
| 自身二聚体 Any/End | 引物自己粘自己的评分。End（3'端）比 Any 更致命，经验阈值 `End < 9`、`Any < 12` |
| 引物间二聚体 Any/End | 正反向引物互粘的评分，判据同上；引物间二聚体比自身二聚体危害更大 |

没有引物输出的序列（结果表中无对应行）：先查 `01_primer_design/primer3_output.txt` 中该序列的 `PRIMER_PAIR_NUM_RETURNED`，为 0 说明 Primer3 判定无合格引物，常见原因是序列过短、GC 极端或重复序列。

## 参数选择建议 | Parameter Guidance

- **常规基因检测/克隆筛验证**：默认参数即可
- **qPCR 引物**：`--product-min-size 70 --product-max-size 200`，Tm 提高到 `--primer-min-tm 60 --primer-opt-tm 62 --primer-max-tm 65`，`--primer-num-return 3` 逐对评估
- **扩增完整 ORF/片段（几百 bp 到数 kb）**：保持 `--method all`，加大 `--primer-end-margin`；产物范围自动覆盖全序列，无需手算
- **在长序列内任意找一段扩增**：`--method random`，并用 `--product-min-size/--product-max-size` 指定期望产物区间
- **含大量 N 的草图序列**：`--primer-max-ns 1` 放宽容忍度，但需自行评估风险
- **序列很短（<150 bp）**：改用 `--method random` 并调小 `--product-min-size`，或考虑合成寡核苷酸而非 PCR 扩增

## 依赖 | Dependencies

- Python 环境随 biopytools 主环境（Biopython、pandas）
- `primer3_core`：默认 `~/miniforge3/envs/misc/bin/primer3_core`（misc 环境），可用 `--primer3-core-path` 指定其他位置或设 `PRIMER3_PATH` 环境变量
- 可选：输出 xlsx 需要 `openpyxl`

## 常见问题 | FAQ

- **改了参数重跑，结果没变？** 断点续传在起作用：`01_primer_design/primer3_output.txt` 已存在时跳过重新计算。换参数重跑请先删除旧输出目录（或仅删该文件）。
- **序列 ID 里的 `|` 或空格去哪了？** 被 `format_sequence_id` 替换为下划线（Primer3 对 ID 格式有限制），结果表以替换后的 ID 呈现。
- **某条序列结果表里没有？** 该序列设计失败。查 `primer3_output.txt` 对应记录的 `PRIMER_PAIR_NUM_RETURNED` 是否为 0，并按"结果解读"末段排查。
- **xlsx 打不开？** 当前 Python 环境缺 `openpyxl`，安装后重跑，或改用 `--output-format csv`。
- **旧版本（v1.0.0）的输出目录还能用吗？** v1.1.0 起结果移入 `01_primer_design/` 子目录、日志移入 `99_logs/`，与旧版平铺结构不同；旧目录重跑会自动补建新结构并重新计算（旧位置的结果文件不会被读取）。
- **产物大小为什么被自动改了？** `--auto-product-size` 默认开启，会按序列长度自动设定产物范围（`all` 模式强制覆盖全序列，`random` 模式受全局上限约束）。要手动控制请加 `--no-auto-product-size`。
