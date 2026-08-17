# HMMsearch - 蛋白质结构域搜索与结果处理 | HMMsearch Domain Search & Result Processing

一句话理解：**拿一个「结构域模型」去一锅蛋白里捞人，找出哪些蛋白带这个结构域**。输入一个 HMM 模型和一组蛋白序列（或直接给已有的搜索结果），输出命中表格、统计摘要，还能把命中的蛋白和结构域片段单独提取出来。

## 功能概述 | Overview

- 两种模式：模式 2 运行 hmmsearch 搜索 + 处理结果（给 HMM + 蛋白）；模式 1 只处理已有的 domtblout 结果文件
- 蛋白序列支持单个 FASTA 文件，也支持一个目录（自动逐文件搜索并合并、生成物种汇总表）
- 解析 domtblout 的 23 列标准格式，输出 CSV / Excel 表格
- 可按 domain E-value / 分数过滤，并提取命中蛋白序列与 domain 片段 FASTA
- 支持 hmmsearch 的 TC/GA/NC 三种模型内置阈值（cut_tc/cut_ga/cut_nc，互斥）
- 断点续传：最终 domtblout 已存在则跳过搜索；多文件时逐文件续传

## 快速开始 | Quick Start

```bash
biopytools hmmsearch -i NB-ARC.hmm -p proteins.fa
```

最小输入：一个 HMM 模型文件（`-i`）+ 一个蛋白 FASTA 文件（`-p`）。结果默认写到 `./hmmsearch_output/`。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| HMM 模型 | 一个「结构域的标准画像」，记录了这类结构域在序列上的保守模式，搜索时用它当模板比对 |
| domtblout | hmmsearch 输出的制表符分隔结果文件，一行一个结构域命中，是本工具解析的核心 |
| E-value | 期望值：这个命中「纯属巧合」的概率，越小越可信（1e-10 比 0.01 强得多） |
| score | 比对得分，越高越像真的；和 E-value 互为表里 |
| TC/GA/NC | 模型作者预设的三档阈值（可信/集合/噪声），比自定 E-value 更省心 |
| domain | 结构域：蛋白上独立折叠、独立功能的一段，一条蛋白可含多个 |

## 输入 | Input

### 模式 2（运行搜索）

- `-i`：HMM 模型文件（如 Pfam、InterPro 导出的 .hmm）
- `-p`：蛋白序列 FASTA 文件或目录（目录下所有 .fa/.faa/.fasta/.fa.gz 等会被逐文件搜索）

### 模式 1（处理已有结果）

- `-i`：hmmsearch 的 `--domtblout` 输出文件（不指定 `-p` 即进入此模式）
- 若要提取序列，仍需 `-p` 提供蛋白 FASTA

domtblout 示例（前几列）：

```text
# target name  accession  tlen  query name  accession  qlen  E-value  score  bias  ...
protein_001  -  500  NB-ARC  PF00931  280  1.2e-40  145.3  0.1  ...
```

## 参数说明 | Parameters

### 主输入与模式 | Main input & mode

**通俗理解|In plain words:** 程序看有没有 `-p` 来自动分模式：给了 `-p` 就是「跑搜索」（`-i` 当 HMM 用），不给就是「只处理已有结果」（`-i` 当 domtblout 用）。

### 输出配置 | Output

**通俗理解|In plain words:** 决定结果写到哪、文件名前缀是什么。默认写到 `./hmmsearch_output/`、前缀 `hmmsearch_results`，**一般不用动**。

### hmmsearch 运行参数 | hmmsearch run options

**通俗理解|In plain words:** 这三个 `--cut-*` 开关让 hmmsearch 直接用模型自带的阈值，比手动设 E-value 更省心；三者互斥，**同时只能用一个**。`--evalue-cutoff`/`--score-cutoff` 是传给 hmmsearch 的报告阈值（决定「报哪些」），一般留空用默认即可。

### 结果过滤 | Result filtering

**通俗理解|In plain words:** 这两个是「搜完之后再筛一遍」：`-e` 按 domain E-value 保留 ≤阈值 的，`-s` 按 domain 分数保留 ≥阈值 的。**默认不筛（全部保留）**；命中太多想收紧时才用。

### 序列提取与输出格式 | Extraction & format

**通俗理解|In plain words:** 默认会额外提取命中蛋白序列和 domain 片段 FASTA、输出 CSV 和 Excel。想只拿表格、省磁盘，可用 `--no-extract-proteins`、`--no-extract-domains`、`--no-csv`、`--no-excel` 关掉。

## 分析流程 | Pipeline

**通俗理解|In plain words:** 模式 2 先跑 hmmsearch 生成 domtblout，再解析、过滤、生成表格、提取序列；模式 1 跳过第一步直接解析。

```text
HMM 模型 + 蛋白 FASTA(或目录)
    │
    ▼
[模式2] 逐文件运行 hmmsearch → 各文件 domtblout → 合并
    │  (断点续传: 最终/单个 domtblout 已存在则跳过)
    ▼
解析 domtblout(23列) → domain 命中列表
    │
    ▼
(可选) 按 E-value / score 过滤
    │
    ▼
生成统计摘要 + CSV/Excel 表格
    │
    ▼
(可选) 提取蛋白序列 + domain 片段 FASTA
```

## 输出 | Output

```text
hmmsearch_output/
├── hmmsearch_results.domtblout            # 合并后的原始 domtblout(模式2)
├── hmmsearch_results_<物种>.domtblout     # 各文件的单独 domtblout(目录模式)
├── hmmsearch_results.csv                  # 命中结果主表格
├── hmmsearch_results.xlsx                 # 命中结果 Excel
├── hmmsearch_results_proteins.fa          # 命中的蛋白序列
├── hmmsearch_results_domains.fa           # 命中的 domain 片段序列
├── hmmsearch_results_species_summary.csv  # 物种汇总表(仅目录/多文件模式)
├── hmmsearch_results_species_summary.xlsx # 物种汇总 Excel(同上)
└── hmmsearch_analysis.log                 # 运行日志(输出目录根下)
```

- 主表格 CSV 含 target_name、query_name、full_evalue、domain_evalue、domain_score、坐标等约 20 列
- `_proteins.fa` 按唯一基因提取；`_domains.fa` 每个 domain 一条，fasta 头含 `query名|坐标|E-value` 信息

## 结果解读 | Interpreting Results

### 1. 主表格（CSV/Excel）

**通俗理解|In plain words:** 一行=一个结构域命中。看两列就能定论：`domain_evalue`（越小越可信）和 `target_name`（命中在哪个蛋白上）。

- `domain_evalue` < 1e-5 通常视为可信命中；`full_evalue` 是该蛋白整体比对的 E-value
- `domain_number` / `domain_total`：该蛋白上这是第几个结构域 / 共几个（同一蛋白可命中多次）
- `env_from`/`env_to` 与 `ali_from`/`ali_to`：结构域在蛋白上的坐标范围

### 2. 统计摘要（运行日志中打印）

**通俗理解|In plain words:** 告诉你有多少命中、多少唯一基因、每个基因带几个结构域。

- 总命中数 > 唯一基因数，说明有的蛋白带了多个该结构域，属正常
- 唯一基因数为 0 表示没有可信命中，检查 HMM 与物种是否匹配

## 参数选择建议 | Parameter Guidance

- **`--cut-ga`**：做「某结构域在哪些蛋白中存在」的定性注释时，优先用它（模型作者校准过的阈值，最省心）
- **`-e`**：结果太松、命中爆炸时设一个 domain E-value 阈值（如 1e-5）收紧；不确定就先不设
- **`--no-extract-domains`**：只想要命中清单、不关心片段序列时关掉，省磁盘
- **目录模式**：`-p` 给目录可一次搜多个物种的蛋白组，自动出物种汇总表

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | — |  | 输入文件：domtblout文件（模式1）或HMM文件（模式2，需同时指定-p）｜Input file: domtblout (mode 1) or HMM file (mode 2, requires -p) |
| `-p, --protein-fasta` | — |  | 蛋白序列FASTA文件或目录（模式2必需，模式1提取序列时需要）｜Protein FASTA file or directory (required for mode 2, needed for mode 1 if extracting sequences) |
| `-o, --output-dir` | `./hmmsearch_output` |  | 输出目录｜Output directory |
| `--output-prefix` | `hmmsearch_results` |  | 输出文件前缀｜Output file prefix |
| `--hmmsearch-path` | `~/miniforge3/envs/protein/bin/hmmsearch` |  | hmmsearch程序路径｜hmmsearch program path |
| `-t, --threads` | `12` | int | 线程数｜Number of threads |
| `--evalue-cutoff` | — | float | hmmsearch报告E-value阈值｜hmmsearch reporting E-value threshold |
| `--score-cutoff` | — | float | hmmsearch报告分数阈值｜hmmsearch reporting score threshold |
| `--cut-tc` | — |  | 使用模型的TC trusted cutoff｜Use model TC trusted cutoff |
| `--cut-ga` | — |  | 使用模型的GA gathering cutoff｜Use model GA gathering cutoff |
| `--cut-nc` | — |  | 使用模型的NC noise cutoff｜Use model NC noise cutoff |
| `-e, --evalue-threshold` | — | float | Domain E-value阈值(保留小于等于该值的)｜Domain E-value threshold (keep <= this value) |
| `-s, --score-threshold` | — | float | Domain分数阈值(保留大于等于该值的)｜Domain score threshold (keep >= this value) |
| `--no-extract-proteins` | — |  | 不提取蛋白序列｜Do not extract protein sequences |
| `--no-extract-domains` | — |  | 不提取domain序列｜Do not extract domain sequences |
| `--no-csv` | — |  | 不输出CSV文件｜Do not output CSV file |
| `--no-excel` | — |  | 不输出Excel文件｜Do not output Excel file |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | — |  | 输入文件：domtblout文件（模式1）或HMM文件（模式2，需同时指定-p）｜Input file: domtblout (mode 1) or HMM file (mode 2, requires -p) |
| `-p, --protein-fasta` | — |  | 蛋白序列FASTA文件或目录（模式2必需，模式1提取序列时需要）｜Protein FASTA file or directory (required for mode 2, needed for mode 1 if extracting sequences) |
| `-o, --output-dir` | `./hmmsearch_output` |  | 输出目录｜Output directory |
| `--output-prefix` | `hmmsearch_results` |  | 输出文件前缀｜Output file prefix |
| `--hmmsearch-path` | `~/miniforge3/envs/protein/bin/hmmsearch` |  | hmmsearch程序路径｜hmmsearch program path |
| `-t, --threads` | `12` | int | 线程数｜Number of threads |
| `--evalue-cutoff` | — | float | hmmsearch报告E-value阈值｜hmmsearch reporting E-value threshold |
| `--score-cutoff` | — | float | hmmsearch报告分数阈值｜hmmsearch reporting score threshold |
| `--cut-tc` | — | store_true | 使用模型的TC trusted cutoff｜Use model TC trusted cutoff |
| `--cut-ga` | — | store_true | 使用模型的GA gathering cutoff｜Use model GA gathering cutoff |
| `--cut-nc` | — | store_true | 使用模型的NC noise cutoff｜Use model NC noise cutoff |
| `-e, --evalue-threshold` | — | float | Domain E-value阈值(保留小于等于该值的)｜Domain E-value threshold (keep <= this value) |
| `-s, --score-threshold` | — | float | Domain分数阈值(保留大于等于该值的)｜Domain score threshold (keep >= this value) |
| `--extract-proteins` | `True` | store_true | 提取匹配的蛋白序列｜Extract matched protein sequences (default: True) |
| `--no-extract-proteins` | — | store_false | 不提取蛋白序列｜Do not extract protein sequences |
| `--extract-domains` | `True` | store_true | 提取domain序列｜Extract domain sequences (default: True) |
| `--no-extract-domains` | — | store_false | 不提取domain序列｜Do not extract domain sequences |
| `--no-csv` | — | store_true | 不输出CSV文件｜Do not output CSV file |
| `--no-excel` | — | store_true | 不输出Excel文件｜Do not output Excel file |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- HMMER 的 hmmsearch（默认 `~/miniforge3/envs/protein/bin/hmmsearch`，通过 conda 环境自动包装调用）
- Python 3：pandas、biopython（SeqIO）、openpyxl（Excel）

## 常见问题 | FAQ

**Q1：`-i` 到底传 HMM 还是 domtblout？**
看有没有 `-p`：给了 `-p` 则 `-i` 必须是 HMM（跑搜索）；不给 `-p` 则 `-i` 必须是已有的 domtblout（只处理结果）。

**Q2：三个 `--cut-*` 能一起用吗？**
不能，三者互斥，同时给会校验报错。选一个即可，通常 `--cut-ga` 最常用。

**Q3：搜索为什么没重跑？**
断点续传按最终 domtblout 是否存在判断。想强制重搜，删掉 `hmmsearch_results.domtblout`（目录模式下还要删对应 `_<物种>.domtblout`）再跑。

**Q4：提取的序列里为什么有重复/缺失？**
蛋白序列按「唯一基因」去重提取；若某个命中的序列 ID 在 FASTA 里找不到，会告警并跳过。请确认 `-p` 提供的 FASTA 与 domtblout 里 target_name 一致。