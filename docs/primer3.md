# primer3 - 批量 PCR 引物设计 | Primer3 Primer Design

一句话理解：**给一批序列批量设计 PCR 引物——自动为每条序列找出一组「正向+反向」引物对，连同长度、退火温度、GC 含量、产物大小等指标一起导出成表格，供实验直接下单合成。**

## 功能概述 | Overview { #overview }

- 用 Primer3 core 批量设计 PCR 引物，一条序列可返回多对（默认 5 对）
- 支持两种策略：all（覆盖序列头尾）与 random（随机设计）
- 自动根据序列长度调整产物大小范围（可关）
- 输出 CSV / TSV / XLSX，表头可选中英文
- 引物长度、退火温度、产物大小等参数全可定制

## 快速开始 | Quick Start { #quick-start }

```bash
biopytools primer3 -i sequences.fasta -o primer3_output
```

最小输入：一个含目标序列的 FASTA 文件。

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语 | 通俗理解 |
|------|----------|
| PCR | 用酶把一段 DNA 大量复制出来，分子生物学的「复印机」 |
| 引物(primer) | 十几到二十几个碱基的短片段，引导复制从哪开始，PCR 成败的关键 |
| 正向/反向引物(F/R) | 一条在目标区间起点、一条在终点，两把「夹子」夹出要扩增的片段 |
| 退火温度(Tm) | 引物与模板结合的温度；太高结合不上、太低乱结合 |
| 产物(amplicon/product) | 两条引物之间被扩增出来的那段 DNA |
| GC 含量 | 序列里 G+C 的比例；过高过低都不利于稳定结合 |
| 二聚体(dimer) | 引物自己或两条引物之间配对，会消耗引物、降低效率 |
| 惩罚值(penalty) | Primer3 给每对引物打的综合分，越低越优 |

## 输入 | Input { #input }

- FASTA 文件(-i)：含要设计引物的序列（DNA）。序列 ID 中的空格和特殊字符会自动替换为下划线（Primer3 要求）。

示例：

```text
>gene1
ATGGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
>gene2
ATGGCTAGCGAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
```

## 参数说明 | Parameters { #parameters }

### 必需参数 | Required

**通俗理解|In plain words:** -i 输入 FASTA，-o 输出目录，两者必填。其余参数都有合理默认，新手可直接用默认。

### 引物长度 | Primer size

**通俗理解|In plain words:** --primer-min-size / --primer-opt-size / --primer-max-size（默认 20/20/22）设定引物长度范围，Primer3 优先取「最优值」。长度太短特异性差、太长容易出错；20 bp 左右是通用值。一般不用动，除非实验有特殊要求。

### 退火温度 | Annealing temperature

**通俗理解|In plain words:** --primer-min-tm / --primer-opt-tm / --primer-max-tm（默认 53/58/63）设定引物退火温度范围，Primer3 优先取最优值 58。想更特异可整体调高，想更容易扩增可调低。一般用默认。

### 产物大小 | Product size

**通俗理解|In plain words:** --product-min-size / --product-max-size（默认 100/300）设定扩增片段长度范围。--auto-product-size 默认开启，会按序列长度自动重新设定产物范围（--product-size-min-ratio 0.5 / --product-size-max-ratio 1.0 是产物占序列长度的比例）；想严格用手动范围时用 --no-auto-product-size 关闭。一般不用动。

### 设计策略 | Design strategy

**通俗理解|In plain words:** --method 选 all（默认，把引物强制设计在序列两端，产物覆盖整条序列，适合「测全长/拿到整段」）或 random（整条序列里随便找最优引物对，适合「只扩增其中一段」）。--primer-end-margin（默认 200）只在 all 模式生效，控制两端允许放引物的范围。一般 all 就够。

### 输出控制 | Output control

**通俗理解|In plain words:** --primer-num-return（默认 5）是每条序列返回的引物对数量，多给几对方便挑；--output-format 选 csv/tsv/xlsx（默认 csv）；--output-header-lang 选表头中文 zh 或英文 en（默认 zh）。按需改。

## 分析流程 | Pipeline { #pipeline }

```text
解析 FASTA（ID 清洗 + 序列大写）
    |
    v
生成 Primer3 输入格式（含全局与每条序列的设定）
    |
    v
conda run primer3_core 批量设计引物
    |
    v
解析 Primer3 输出（逐引物对提取序列/Tm/GC/产物/惩罚值）
    |
    v
格式化并保存 primers_result.<csv|tsv|xlsx>
```

## 输出 | Output { #output }

```text
primer3_output/
|-- primers_result.csv        # 引物结果表（核心，格式由 --output-format 决定）
-- primer3_design.log         # 运行日志
```

结果表主要列（中文表头）：序列ID、输入序列、引物对编号、正向引物、反向引物、正向/反向引物长度、正向/反向引物 Tm 值、退火温度、正向/反向 GC 含量(%)、产物大小(bp)、产物 Tm 值、惩罚值、各类二聚体指标。

## 结果解读 | Interpreting Results { #interpreting-results }

**通俗理解|In plain words:** 打开结果表，每行是一对引物。先看有没有行（没行=没设计出来），再看惩罚值（越低越优）和各项指标是否落在合理范围。

- 正向引物 / 反向引物：直接拿去合成的序列
- 退火温度：取正反引物 Tm 的较小值，做 PCR 时参考
- GC 含量：通常在 40-60% 之间比较理想，偏离太多扩增困难
- 产物大小(bp)：扩增片段长度，应符合预期
- 惩罚值：Primer3 综合打分，越小越优，同一序列内选惩罚值最小的那对
- 二聚体指标：数值越低越好，高说明引物容易自己配对

好坏判据：每条序列能返回若干对引物、惩罚值较低、GC 与 Tm 在常规范围即算成功；完全没结果多半是序列太短或含 N 过多（见 FAQ）。

## 参数选择建议 | Parameter Guidance { #parameter-guidance }

- 常规 PCR：全部默认（引物 20-22 bp、Tm 53-63、产物 100-300 bp）
- 想扩增整段序列：默认 method=all 即可
- 只扩增其中一段：--method random 并手动设 --product-min-size / --product-max-size
- 想要更多候选：调大 --primer-num-return
- 给 Excel 用：--output-format xlsx；给非中文同事：--output-header-lang en

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input-fasta, -i` | 必填 |  | 输入FASTA文件｜Input FASTA file path |
| `--output-dir, -o` | 必填 | Path | 输出目录｜Output directory |
| `--primer-min-size` | `20` | int | 最小引物长度｜Minimum primer size |
| `--primer-opt-size` | `20` | int | 最优引物长度｜Optimal primer size |
| `--primer-max-size` | `22` | int | 最大引物长度｜Maximum primer size |
| `--primer-min-tm` | `53.0` | float | 最小退火温度(°C)｜Minimum annealing temperature (°C) |
| `--primer-opt-tm` | `58.0` | float | 最优退火温度(°C)｜Optimal annealing temperature (°C) |
| `--primer-max-tm` | `63.0` | float | 最大退火温度(°C)｜Maximum annealing temperature (°C) |
| `--product-min-size` | `100` | int | 最小产物大小(bp)｜Minimum product size (bp) |
| `--product-max-size` | `300` | int | 最大产物大小(bp)｜Maximum product size (bp) |
| `--primer-num-return` | `5` | int | 返回引物对数量｜Number of primer pairs to return |
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

## 依赖 | Dependencies { #dependencies }

- Primer3（命令行 primer3_core），默认 conda 环境 misc（~/miniforge3/envs/misc/bin/primer3_core），通过 conda run 自动检测调用
- Biopython（解析 FASTA）
- pandas（结果格式化）；输出 xlsx 额外需要 openpyxl
- Python 3

## 常见问题 | FAQ { #faq }

Q1：支持断点续传吗？
不支持。每次运行都全量重新解析、重新设计并覆盖输出文件。

Q2：method=all 和 random 有什么区别？
all 把引物强制设计在序列两端（用 SEQUENCE_PRIMER_PAIR_OK_REGION_LIST 限制），产物覆盖整条序列；random 不限制位置，整条序列里找最优引物对。默认 all。

Q3：为什么某条序列没设计出引物？
可能原因：序列太短（容不下引物+产物）、含 N 过多（默认不接受 N，PRIMER_MAX_NS_ACCEPTED=0）、GC 极端、或产物大小范围与序列长度冲突。检查该序列长度与内容。

Q4：auto-product-size 是什么？
默认开启，按序列长度自动设定产物大小范围（最小 = max(序列长度×0.5, 全局最小值)，method=all 时最大 = 序列长度×1.0）。想要严格手动范围就加 --no-auto-product-size。

Q5：primer3_core 找不到怎么办？
默认按「PRIMER3_PATH 环境变量、配置文件、~/miniforge3/envs/misc/bin/primer3_core」查找；找不到会直接报错，需确认装了 Primer3 并调整路径。

Q6：输出 xlsx 报错？
需要 openpyxl 库；未安装时改用 --output-format csv 或 tsv。