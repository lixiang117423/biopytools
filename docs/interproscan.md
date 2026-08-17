# interproscan — 蛋白质功能注释 | InterProScan Protein Function Annotation

一句话理解：**给每条蛋白质序列「贴标签」**——自动调用 InterProScan 识别蛋白上的结构域、家族与功能位点，并把零散的结果整理成一张张可直接看的 Excel / CSV 报告（可选带 GO 术语与代谢通路）。

## 功能概述 | Overview

- 封装 InterProScan，一键对蛋白质 FASTA 做功能结构域注释
- 支持同时输出多种格式（默认 TSV + XML，TSV 解析更可靠）
- 自动把 InterProScan 原始结果整理成 Excel / CSV 报告（详细注释、蛋白汇总、数据库统计等）
- 可选获取 GO 术语与 Pathway 信息，并内置 GO 数据库填充名称与 ontology
- 默认禁用预计算查找服务（`-dp`），避免离线环境卡网络

## 快速开始 | Quick Start

```bash
biopytools interproscan -i proteins.fa -o results
```

最小输入：一个蛋白质 FASTA 文件 + 一个输出前缀（`-o` 不带扩展名，程序会据此派生各类输出文件名）。

## 零基础概念速览 | Concepts in plain words

不熟悉生信术语的话，先花两分钟看这张表：

| 术语 | 通俗理解 |
|------|----------|
| 结构域 (domain) | 蛋白上能独立折叠、承担某种功能的一段「零件」，同一零件可装在不同蛋白上 |
| InterPro | 一个汇总了多个蛋白家族/结构域数据库的「大目录」，给每个零件一个统一编号 |
| 签名 (signature) | 某个家族/结构域的「指纹特征」，比对上了就说明蛋白含这个零件 |
| GO 术语 | 用一套标准词汇描述「这个蛋白在干什么」（分子功能 / 生物过程 / 细胞定位） |
| Pathway | 代谢通路，即蛋白参与的一连串生化反应 |
| 预计算查找 (precalc) | InterProScan 把常见蛋白的结果缓存到服务器，联网时可直接查而不用重算 |

## 输入 | Input

### 蛋白质 FASTA 文件（-i, --input）

氨基酸序列文件，标准 FASTA 格式。本工具不适用于 DNA/RNA 序列——那是另一个流程。

```text
>protein1
MKTLSVKKL...
>protein2
MAEIQTER...
```

## 参数说明 | Parameters

### 必需参数 | Required

**通俗理解|In plain words:** 两个必填项：输入蛋白文件和输出前缀。输出前缀不带扩展名，程序会在它后面拼出 `.tsv`、`.xml`、`_formatted_report.xlsx` 等一系列文件。

### 软件与运行 | Software & runtime

**通俗理解|In plain words:** `--interproscan-path` 指向 InterProScan 安装脚本（一般用默认路径，装在不同位置才改）；`--threads` 越大跑得越快但越占 CPU；`--format` 决定原始输出格式（默认 TSV,XML 就够，TSV 用于解析报告）；`--python-path` 是给 InterProScan 指定兼容的 Python 解释器（见 FAQ，Python 3.12+ 时才需要动）。

### 注释内容开关 | Annotation toggles

**通俗理解|In plain words:** `--goterms` / `--pathways` 控制要不要额外导出 GO 术语和代谢通路信息，默认都关（省时间）。需要下游做富集分析时再开。`--applications` 可以只跑指定的数据库（如只跑 Pfam）以提速，不给则全跑。

### 报告整理 | Report formatting

**通俗理解|In plain words:** 这些管「结果怎么整理给你」。默认生成 Excel + CSV 双格式报告并含汇总统计表，一般不用动；只想看原始输出时可用 `--no-report` 关闭整理。

## 分析流程 | Pipeline

**通俗理解|In plain words:** 三步：先调 InterProScan 对蛋白做注释，再解析它的输出，最后整理成报告。

```text
输入蛋白质 FASTA
    │
    ▼
1. 调用 InterProScan（-dp 禁用预计算，输出 TSV + XML）
    │
    ▼
2. 解析结果（优先 TSV，可选 GO/Pathway）
    │
    ▼
3. 生成 Excel / CSV 报告（详细注释 / 蛋白汇总 / 数据库统计）
```

## 输出 | Output

以 `-o results` 为例：

```text
results.tsv                             # InterProScan 原始 TSV
results.xml                             # InterProScan 原始 XML
results_formatted_report.xlsx           # Excel 报告（核心）
results_detailed_results.csv            # 详细注释 CSV
results_protein_summary.csv             # 蛋白注释汇总 CSV
results_database_stats.csv              # 数据库注释统计 CSV
results_go_annotations.csv              # GO 注释 CSV（仅 --goterms）
results_pathway_annotations.csv         # Pathway 注释 CSV（仅 --pathways）
```

### 关键文件说明 | Key files

- `results_formatted_report.xlsx`：多 sheet 的 Excel 报告——Detailed Results（每条匹配）、Protein Summary（每蛋白汇总）、Database Statistics（各数据库命中数）、以及开启对应开关后的 GO / Pathway sheet
- `results_detailed_results.csv`：每条蛋白匹配的明细，列含 Protein ID、Database、Signature ID、Description、Start、End、Score、InterPro ID、GO Terms、Pathways 等
- `results.tsv` / `results.xml`：InterProScan 的原始输出，方便用其它工具二次处理

## 结果解读 | Interpreting Results

**通俗理解|In plain words:** 核心看两件事——「每条蛋白匹配到了哪些结构域/家族」和「整体上哪个数据库命中最多」。匹配越多、得分越高，说明该蛋白的功能注释越丰富可信。

- **Detailed Results**：一行 = 一条蛋白在某个数据库的一条匹配。Score 越高越可信；Start / End 是匹配落在蛋白上的区段
- **Protein Summary**：按 Total Matches 降序排列，一眼看出哪些蛋白注释最丰富；注释全空的蛋白（0 匹配）说明没找到已知结构域
- **Database Statistics**：各数据库（Pfam、SMART、Gene3D 等）的命中数，可评估注释覆盖面
- **GO / Pathway sheet**：开启对应开关后才有，每个 GO / Pathway 一行展开，便于下游富集分析

## 参数选择建议 | Parameter Guidance

- **常规注释**：默认参数直接跑，`-t` 按机器核数给（如 32）
- **要快**：用 `--applications Pfam` 只跑最常用的库；或关掉报告整理 `--no-report`
- **要下游富集分析**：加 `--goterms --pathways`
- **Python 3.12+ 环境**：用 `--python-path` 指向 3.8–3.11 的解释器（见 FAQ）

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input, -i` | 必填 |  | 蛋白质FASTA文件｜Protein FASTA file path |
| `--output-prefix, -o` | 必填 |  | 输出文件前缀｜Output file prefix (without extension) |
| `--interproscan-path, -a` | `~/software/InterProScan/v.5.75-106.0/interproscan-5.75-106.0/interproscan.sh` |  | InterProScan软件路径｜InterProScan software installation path |
| `--format, -f` | `TSV,XML` |  | 输出格式,支持逗号分隔多格式(如 TSV,XML)｜Output format, supports comma-separated multi-format (e.g. TSV,XML) |
| `--threads, -t` | `12` | int | 线程数｜Number of threads |
| `--disable-precalc` | `True` |  | 禁用预计算查找服务(默认)｜Disable precalc lookup service (default) |
| `--enable-precalc` | — |  | 启用预计算查找服务(需要网络)｜Enable precalc lookup service (requires network) |
| `--goterms` | — |  | 获取GO术语(默认不获取)｜Get GO terms (disabled by default) |
| `--pathways` | — |  | 获取Pathway信息(默认不获取)｜Get pathway information (disabled by default) |
| `--applications, -appl` | — |  | 指定运行的应用，逗号分隔｜Specify applications to run, comma-separated |
| `--temp-dir` | — | Path | 临时目录路径｜Temporary directory path |
| `--python-path` | — | Path | Python解释器路径（兼容Python 3.8-3.11）｜Python interpreter path (compatible with Python 3.8-3.11) |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 输入蛋白质FASTA文件路径｜Input protein FASTA file path |
| `-o, --output-prefix` | 必填 |  | 输出文件前缀(不含扩展名)｜Output file prefix (without extension) |
| `--interproscan-path` | `~/software/InterProScan/v.5.75-106.0/interproscan-5.75-106.0/interproscan.sh` |  | InterProScan软件路径｜InterProScan software path |
| `-f, --format` | `TSV,XML` |  | 输出格式，支持多个格式逗号分隔 (例如: TSV,XML,JSON)｜Output format, support multiple formats comma-separated (e.g., TSV,XML,JSON) |
| `-t, --threads` | `24` | int | 线程数｜Number of threads |
| `--disable-precalc` | `True` | store_true | 禁用预计算查找服务(默认启用，解决网络问题)｜Disable precalc lookup service (enabled by default, solves network issues) |
| `--enable-precalc` | — | store_true | 启用预计算查找服务(需要网络连接)｜Enable precalc lookup service (requires network connection) |
| `--goterms` | — | store_true | 获取GO术语(默认不获取)｜Get GO terms (disabled by default) |
| `--pathways` | — | store_true | 获取Pathway信息(默认不获取)｜Get pathway information (disabled by default) |
| `-appl, --applications` | — |  | 指定运行的应用，逗号分隔(例如: Pfam,SMART,Gene3D)｜Specify applications to run, comma-separated (e.g., Pfam,SMART,Gene3D) |
| `--temp-dir` | — |  | 临时目录路径｜Temporary directory path |
| `--no-report` | — | store_true | 不生成整理后的报告｜Do not generate formatted report |
| `--report-format` | `both` | excel/csv/both | 报告格式 (默认: both，同时生成Excel和CSV)｜Report format (default: both, generate both Excel and CSV) |
| `--no-summary` | — | store_true | 不包含汇总统计表｜Do not include summary statistics |
| `--go-database` | — |  | GO术语JSON数据库文件路径（用于填充GO名称和ontology）｜GO term JSON database file path (for filling GO names and ontologies) |
| `--python-path` | — |  | Python解释器路径（用于指定兼容的Python 3.8-3.11版本）｜Python interpreter path (for specifying compatible Python 3.8-3.11) |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- InterProScan（默认路径 `~/software/InterProScan/v.5.75-106.0/interproscan-5.75-106.0/interproscan.sh`）
- Python 3.8–3.11（InterProScan 5.70+ 不兼容 3.12+）
- pandas + openpyxl（生成报告用；openpyxl 缺失时自动跳过 Excel、只出 CSV）
- 内置 GO 术语数据库（可用 `--go-database` 指定外部 JSON 覆盖）

## 常见问题 | FAQ

**Q1：报错提示 Python 版本不兼容？**
InterProScan 5.70 以上要求 Python 3.8–3.11。若系统是 Python 3.12+，程序会尝试自动查找兼容版本；找不到就用 `--python-path /path/to/python3.11` 显式指定。

**Q2：为什么默认禁用预计算查找服务（-dp）？**
预计算服务需要联网查询 InterProScan 的服务器，离线/受限环境会卡住或报错。默认 `-dp` 跳过联网，直接本地计算，最稳。

**Q3：Excel 报告没生成，只有 CSV？**
说明环境里没装 openpyxl。装一下（`pip install openpyxl`）即可恢复 Excel 输出，CSV 内容不受影响。

**Q4：能不能只跑某几个数据库？**
可以，用 `--applications Pfam,SMART` 逗号分隔指定，能大幅缩短运行时间。

**Q5：GO / Pathway 结果为什么是空的？**
这两个信息默认不导出，需要显式加 `--goterms` / `--pathways` 才会出现在输出里。
