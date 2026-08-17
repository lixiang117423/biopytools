# vcf-sample-hete - VCF基因型统计 | VCF Genotype Statistics

一句话理解：**对 VCF 里每个样本的基因型做一次「体检」，算出杂合率、纯合率、缺失率等指标**，帮你快速揪出测序质量差（缺失多）或异常（杂合率过高）的样本。

## 功能概述 | Overview { #overview }

- 逐样本统计基因型分类：参考纯合(0/0)、杂合(0/1)、变异纯合(1/1)、缺失(./.)
- 支持定相(0|1)与未定相(0/1)基因型，以及多等位(0/2、1/2 等)
- 按深度、质量、缺失三类条件过滤
- 输出汇总表、简化比率表、详细表、逐样本报告多级结果
- 自动标记潜在异常样本（调用率 <0.9、杂合率 >0.7）
- 纯 Python 实现，不调用任何外部命令行工具

## 快速开始 | Quick Start { #quick-start }

```bash
biopytools vcf-sample-hete -i variants.vcf -o vcf_stats_output
```

统计 `variants.vcf` 里每个样本的杂合/纯合/缺失等指标，结果写到 `vcf_stats_output/`。

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语 | 通俗理解 |
|------|----------|
| 基因型 | 一个位点上两条同源染色体的构成，0 代表参考、1 代表变异 |
| 纯合 | 两条染色体一样（0/0 或 1/1），像两只手都是左手或都是右手 |
| 杂合 | 两条染色体不一样（0/1），像一左一右 |
| 缺失 | 这个位点没测出来（./.），数据里是空白 |
| 定相/未定相 | 用 \| 还是 / 分隔两个等位基因，只是写法不同，本工具一视同仁 |
| QUAL | 每个位点的「质量分」，分数越低越可能是测序错误 |
| DP | 覆盖深度，这个位点被读了多少次，太浅则不可信 |

## 输入 | Input { #input }

输入一个 VCF 文件，支持 `.vcf` 或 `.vcf.gz`（gzip 自动识别解压）。基因型可带 FORMAT 字段（如 `0/1:20:10,10`），本工具会自动只取 GT 部分。

## 参数说明 | Parameters { #parameters }

### 必需参数 | Required

**通俗理解|In plain words:** 指定输入 VCF，每次必填。

### 输出 | Output

**通俗理解|In plain words:** 决定结果写到哪个目录，默认 `vcf_stats_output`。

### 过滤 | Filtering

**通俗理解|In plain words:** 决定「什么样的数据信不过、直接不计数」。`-d` 按覆盖深度滤掉太浅的基因型（需要 FORMAT 里有 DP 才生效），`-q` 按质量分滤掉低质量位点，`-e` 把缺失基因型从统计里剔除。**三项默认都不过滤，绝大多数情况直接用默认即可。**

### 输出控制 | Output control

**通俗理解|In plain words:** 控制生成哪些报告文件。默认全生成；`-D` 省掉详细统计表、`-S` 省掉汇总表（适合只要快速看一眼比率）。**一般不用动。**

### 日志 | Logging

**通俗理解|In plain words:** 控制日志写到哪、多啰嗦。`--log-file` 指定日志文件，`--log-level` 指定级别（DEBUG/INFO/WARNING/ERROR/CRITICAL）。**一般不用动。**

## 分析流程 | Pipeline { #pipeline }

```text
输入 VCF
    │
    ▼
步骤1: 解析表头 -> 得到样本名列表
    │
    ▼
步骤2: 逐位点解析基因型，套用深度/质量/缺失过滤
    │
    ▼
步骤3: 计算每个样本的计数与比率(杂合率/纯合率/缺失率/调用率)
    │
    ▼
步骤4: 导出多级报告 + 生成分析总结(含异常样本标记)
```

## 输出 | Output { #output }

```text
vcf_stats_output/
├── 00_pipeline_info/
│   └── software_versions.yml            # 运行环境版本与参数记录
├── 99_logs/
│   └── vcf_stats_YYYYMMDD_HHMMSS.log    # 运行日志(带时间戳)
├── genotype_summary_statistics.txt      # 汇总统计表(每个样本一行,含全部比率)
├── genotype_rates_simple.txt            # 简化比率表(只看主要比率)
├── genotype_detailed_statistics.txt     # 详细基因型分布表(用 -D 时不生成)
├── per_sample_stats/                    # 每个样本一份报告
│   ├── {sample1}_genotype_stats.txt
│   └── ...
└── analysis_summary.txt                 # 分析总结 + 异常样本标记
```

## 结果解读 | Interpreting Results { #interpreting-results }

- **`genotype_summary_statistics.txt`**：核心结果表，每行一个样本，列含总位点数、有效/缺失调用数、调用率(Call_Rate)、缺失率(Missing_Rate)、杂合率(Heterozygosity_Rate)、纯合率(Homozygosity_Rate)等。
- **`genotype_rates_simple.txt`**：只保留 0/0、0/1、1/1、./. 四类计数与杂合率/纯合率，便于快速扫一眼。
- **`analysis_summary.txt`**：末尾列出「潜在异常样本」——调用率 <0.9（数据缺太多）或杂合率 >0.7（疑似污染/异常）的样本会点名。
- 好坏判据：健康的样本调用率通常 >0.9、缺失率 <0.1；杂合率取决于物种与群体（自然群体通常明显低于近交系），**同一批样本里明显偏离主流的个体最值得检查**。

## 参数选择建议 | Parameter Guidance { #parameter-guidance }

- 先不设任何过滤跑一遍看全貌，再决定是否过滤。
- **`-d/--min-depth`**：只有 FORMAT 含 DP 字段时才有效；想滤掉低覆盖噪声时给一个合理阈值（如 5）。
- **`-q/--min-qual`**：想排除低质量位点干扰时设置，按你 VCF 的 QUAL 分布取值。
- **`-e/--exclude-missing`**：计算比率时想把「没测到」的部分剔除出分母时使用；默认缺失计入分母。

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 | Path | 输入VCF文件路径｜Input VCF file path |
| `--output, -o` | `vcf_stats_output` | str | 输出目录｜Output directory |
| `--min-depth, -d` | `0` | int | 最小深度过滤阈值｜Minimum depth filter threshold |
| `--min-qual, -q` | `0.0` | float | 最小质量分数过滤阈值｜Minimum quality score filter threshold |
| `--exclude-missing, -e` | — |  | 排除缺失基因型｜Exclude missing genotypes |
| `--no-detailed, -D` | — |  | 不输出详细统计结果｜Do not output detailed statistics |
| `--no-summary, -S` | — |  | 不输出汇总统计结果｜Do not output summary statistics |
| `--verbose` | — |  | 增加输出详细程度｜Increase output verbosity |
| `--quiet` | — |  | 静默模式｜Quiet mode |
| `--log-file` | — | str | 日志文件路径｜Log file path |
| `--log-level` | `INFO` | DEBUG/INFO/WARNING/ERROR/CRITICAL | 日志级别｜Log level |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 输入VCF文件路径｜Input VCF file path |
| `-o, --output` | `vcf_stats_output` |  | 输出目录｜Output directory |
| `-d, --min-depth` | `0` | int | 最小深度过滤阈值｜Minimum depth filter threshold |
| `-q, --min-qual` | `0.0` | float | 最小质量分数过滤阈值｜Minimum quality score filter threshold |
| `-e, --exclude-missing` | — | store_true | 排除缺失基因型｜Exclude missing genotypes |
| `-D, --no-detailed` | — | store_true | 不输出详细统计结果｜Do not output detailed statistics |
| `-S, --no-summary` | — | store_true | 不输出汇总统计结果｜Do not output summary statistics |
| `--verbose, -v` | `0` | count | 增加输出详细程度｜Increase output verbosity |
| `--quiet` | — | store_true | 静默模式｜Quiet mode |
| `--log-file` | — | str | 日志文件路径｜Log file path |
| `--log-level` | `INFO` | DEBUG/INFO/WARNING/ERROR/CRITICAL | 日志级别｜Log level |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies { #dependencies }

- Python 3
- pandas（纯 Python 统计，仅用于结果表格输出；无外部命令行工具依赖）

## 常见问题 | FAQ { #faq }

**Q1：有断点续传吗？**
没有。每次运行都会重新读取整个 VCF 并覆盖输出目录里的报告文件。

**Q2：设置了 `-d` 深度过滤，为什么没效果？**
深度过滤需要每个样本的 FORMAT 里有 DP 字段。若你的 VCF 没有 DP 信息，`-d` 不会起任何作用。

**Q3：`-e` 排除缺失会改变什么？**
默认缺失基因型计入统计总数；加 `-e` 后，缺失的位点不进入分子也不进入分母，杂合率/纯合率的分母变小、数值会变。比较不同批次结果时要保持一致口径。

**Q4：定相(0|1)和多等位(0/2、1/2)能处理吗？**
能。定相与未定相统一处理；多等位里两个等位基因相同归纯合、不同归杂合，会归入对应类别。
