# 基因组装配统计 | Genome Assembly Statistics

一句话理解：**快速摸清一个 FASTA/FASTQ 文件里「有多少条序列、总共多长、N50 是多少」**，一口气给出组装常用的长度统计指标，还能批量处理一个目录。

## 功能概述 | Overview { #overview }

- 对 FASTA / FASTQ 计算：总长（sum）、序列数（n）、平均长（ave）、最长（largest）、N50–N100（含各 Nx 对应序列数）
- 统计 N 碱基总数（N_count）与连续 N 段数量（Gaps）
- 支持单文件或整个目录（自动发现 `.fa/.fasta/.fna/.fq/.fastq` 及大写扩展名）
- 结果同时打印到屏幕，并生成每个文件的 CSV / Excel 报告 + 一份汇总表
- 纯 Python + pandas 实现，无外部生信软件依赖

## 快速开始 | Quick Start { #quick-start }

```bash
biopytools assembly-stats -i genome.fa
```

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语<br>Term | 通俗理解 |
|---|---|
| N50 | 把序列从长到短排，累加到总长 50% 时最后加进来的那条序列的长度；N50 越大组装越「连贯」 |
| Nx（N50/N60/…/N100） | 同一套玩法，只是累加目标从 50% 换成 60%…100% |
| contig / scaffold | 没有缺口的连续序列 / 由 contig 加缺口拼成的大段序列 |
| Gap（缺口） | 序列里连续的一段 N（不确定碱基） |
| 平均长度 | 总长除以序列条数，粗略反映每条序列的典型长度 |

## 输入 | Input { #input }

- `-i` 可以是单个文件，也可以是一个目录（目录下按扩展名自动收集所有 FASTA/FASTQ）
- FASTA 取 `>` 后第一个词作序列名；FASTQ 按四行一条读取
- 长度小于 `-l/--min-length` 的序列会被过滤掉（默认 1，即不过滤）

```text
>chr1
ACGTACGTACGTNNNACGT
>chr2
TTTTGGGGCCCC
```

## 参数说明 | Parameters { #parameters }

### 长度过滤 | Length filter

**通俗理解|In plain words:** `-l, --min-length` 设一个「最低长度门槛」，短于它的序列直接不参与统计。**默认 1 即不过滤；想只看大 contig、排除碎序列时再调大**，比如设 1000 只看 ≥1kb 的序列。

相关参数：`-l, --min-length`。

### 输出格式 | Output format

**通俗理解|In plain words:** 控制屏幕打印样式，不影响生成的 CSV/Excel 报告。`-s` 输出「grep 友好」格式（一行一个键值，方便 `grep N50` 抓取）；`-t` 改成 Tab 分隔的一行；`-u` 是 `-t` 且不带表头。**只想在终端看结果就不加；要写脚本抓字段用 `-s` 或 `-t`。**

相关参数：`-s`、`-t`、`-u`。

## 输出 | Output { #output }

```text
assembly_stats_output/
├── {文件名}_stats.csv            # 每个文件单独的竖版报告（键值两列）
├── {文件名}_stats.xlsx           # 每个文件单独的 Excel 报告
├── assembly_stats_summary.csv    # 所有文件汇总表（每行一个文件）
├── assembly_stats_summary.xlsx   # 汇总 Excel
└── assembly_stats.log            # 运行日志
```

屏幕输出示例（竖版）：

```text
File|文件	genome.fa
Sum|总和	23328019
N|序列数	16
Average|平均长度	1458001
Largest|最大长度	3291936
N50	1687656
N50_count|N50_序列数	5
...
N_count|N碱基数	0
Gaps|Gap数	0
```

## 结果解读 | Interpreting Results { #interpreting-results }

- `N50` 越大，说明一半以上的组装被少数几条长序列覆盖，组装越连贯（同一物种 N50 翻倍通常意味着明显更好的组装）
- `N|序列数` 越小越好：在总长相近时，序列条数少 = 更少碎片
- `Gaps` 与 `N_count` 反映缺口情况：两者都为 0 说明没有不确定碱基；缺口多则组装有未填的洞
- `Largest` 是最长序列长度，约等于最长染色体的规模

## 参数选择建议 | Parameter Guidance { #parameter-guidance }

- 快速看一个文件：`-i genome.fa`
- 批量统计一个目录：`-i assemblies/`（自动收集并生成汇总表）
- 只想看大序列：`-l 1000`
- 写脚本抓字段：`-s`（grep 友好）或 `-u`（Tab 无表头）

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 输入文件或文件夹｜Input file or directory path |
| `-l, --min-length` | `1` | int | 最小序列长度过滤｜Minimum sequence length cutoff |
| `-s` | — |  | Grep友好输出格式｜Print grep-friendly output |
| `-t` | — |  | Tab分隔输出｜Print tab-delimited output |
| `-u` | — |  | Tab分隔输出且无header｜Print tab-delimited output without header |
| `-o, --output-dir` | `./assembly_stats_output` |  | 输出目录｜Output directory |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 输入文件或文件夹｜Input file or directory path |
| `-l, --min-length` | `1` | int | 最小序列长度过滤｜Minimum sequence length cutoff |
| `-s` | — | store_true | Grep友好输出格式｜Print grep-friendly output |
| `-t` | — | store_true | Tab分隔输出｜Print tab-delimited output |
| `-u` | — | store_true | Tab分隔输出且无header｜Print tab-delimited output without header |
| `-o, --output-dir` | `./assembly_stats_output` |  | 输出目录｜Output directory |
| `--version` | — | version |  |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies { #dependencies }

- Python 3
- pandas（报告生成）
- openpyxl（生成 `xlsx` 报告时需要）
- 无 conda 环境、无外部生信软件依赖

## 常见问题 | FAQ { #faq }

**Q1：目录里没输出任何文件？**
目录模式下只认 `.fa/.fasta/.fna/.fq/.fastq` 及对应大写扩展名。确认文件扩展名在列表内，或改成直接传文件路径。

**Q2：为什么统计的序列数比文件里少？**
`-l/--min-length` 过滤掉了短于门槛的序列，默认 1 本应不过滤，若设大了会变少。检查该参数。

**Q3：N50 和 Largest 有什么区别？**
Largest 是最长一条序列的长度（单条记录）；N50 是「累加到总长 50% 时那条序列的长度」，是整体连续度的指标，通常小于 Largest。

**Q4：会断点续传吗？**
不会。统计是一次性完成的，重跑会覆盖报告文件。
