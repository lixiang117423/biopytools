# blast - BLAST 序列比对分析 | BLAST Sequence Alignment Analysis

一句话理解：**把你的一批序列（样品）拿去和一个"目标库"（如候选基因）做 BLAST 比对，自动完成建库、比对、排序、筛选、统计和可视化，告诉你每条查询序列最像目标库里的谁、像到什么程度。**

## 功能概述 | Overview

- 封装 BLAST+ 完整流程：建库(makeblastdb) → 逐样本比对 → 合并 → 排序 → 高质量筛选 → 统计报告 → 可视化
- 支持 5 种 BLAST 程序（`blastn/blastp/blastx/tblastn/tblastx`），**默认自动检测**（根据查询与参考的序列类型）
- 输入支持单文件、目录、或样品映射文件（`-s`）
- 输出分层目录（`00_pipeline_info / 01_database / 02_blast / 03_alignments / 99_logs`），含 TSV、Excel、HTML/文本比对可视化
- 断点续传：已建库、已比对的样品默认跳过，`--force` 强制重算
- 自动按"覆盖度降序、相似度降序、E-value 升序"排序，并筛出高质量比对

## 快速开始 | Quick Start

~~~bash
biopytools blast -i sequences/ -r nlr_genes.fa -o blast_results
~~~

最小输入：一个查询目录（`-i`）+ 一个目标序列文件（`-r`，必需）。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| 查询(query) / 目标(subject) | 你想查的序列 / 拿来"对答案"的库；结果里 q 开头是查询侧，s 开头是目标侧 |
| 相似度(identity) | 两条序列对齐后有多少比例完全一样，越高越像 |
| 覆盖度(coverage) | 查询序列有多长被这段比对"盖住"（比对跨度÷查询序列长度），越满说明整条查询都对上了 |
| E-value | "纯靠随机也能碰到这么像"的概率，越小越可信 |
| Bit score | 比对的"力度分"，越大越强 |
| blastn/blastp/blastx/tblastn/tblastx | 核酸比核酸 / 蛋白比蛋白 / 核酸翻译后比蛋白 / 蛋白比核酸翻译 / 两边都翻译后比 |
| word_size | BLAST 找种子的"词长"，小=更敏感但慢，程序会自动设好 |

## 输入 | Input

- 查询：单个 FASTA 文件，或目录（按 `--input-suffix` 匹配，默认 `*.fa`），或用 `-s 样品映射文件`
- 目标库：单个序列 FASTA（`-r`，必需），工具会自动 `makeblastdb` 建库
- 样品映射文件（`-s`）两列：文件路径 + 样品名，制表符分隔

~~~text
# 样品映射文件 sample_map.tsv
/path/to/sample1.fa    sample1
/path/to/sample2.fa    sample2
~~~

- 序列类型（核酸/蛋白）自动检测；也可用 `--blast-type` 和 `--target-db-type` 显式指定

## 参数说明 | Parameters

### 输入与输出 | Input & output

**通俗理解|In plain words:** `-i`（查询）、`-s`（样品映射，与 `-i` 二选一）、`-r`（目标库，必填）、`-o`（输出目录，默认 `./blast_output`）。样品名默认从文件名自动提取（`--sample-name-pattern` 控制），单文件可用 `--sample-name` 手动命名。

### BLAST 类型与阈值 | BLAST type & thresholds

**通俗理解|In plain words:** `--blast-type` 决定用哪个 BLAST 程序，**默认自动检测**，一般不用管。`-e`（默认 1e-5）是 E-value 门槛，调小=只要更可信的；`--max-target-seqs`（默认 10）每条查询最多保留几个目标；`--word-size` 按程序自动设，一般不用动。`--task` 只对 blastn 生效，控制搜索"灵敏度档位"：`blastn`（默认）最敏感，能报出分歧大的同源（如 60–70% 相似度）；`megablast` 最快但只适合找高度相似的序列（实测报不出 <70% 的比对）；`dc-megablast` 居中，适合跨物种。`--min-identity`（默认 30%）会传给 blastn 的 `-perc_identity` 在比对阶段直接过滤，并同时过滤汇总表。

### 高质量筛选 | High quality filter

**通俗理解|In plain words:** `--high-quality-evalue`（默认 1e-10）配合 `--min-identity`、`--min-coverage` 一起，从排序结果里再挑出"既像又满又可信"的顶级命中（查询覆盖度 ≥50% 等）。阈值收紧=结果更精但可能漏掉边缘的弱同源。`--min-coverage` 按查询覆盖度算（整条查询被盖住多少），衡量"是不是全长同源"。

### 比对可视化 | Alignment visualization

**通俗理解|In plain words:** `--alignment-output`（默认 `both`）控制是否生成比对可视化（`none/text/html/both`）。`--alignment-width` 是每行显示的字符数，`--alignment-max-per-query`（默认 5）**每条查询序列**只画相似度最高的 5 条，保证每个基因都会出现、不被高相似的头部基因挤掉；`--alignment-max-per-sample`（默认 2000）是整个样品的总上限，防止 HTML 过大。`--html-theme` 换 HTML 配色。`--merge-html`（默认开启）把 HTML 合并成**单个自包含文件**（概览 + 每样品一个 tab），一个文件就能发给别人看；`--no-merge-html` 恢复旧的 index + 分样品多文件。**只想拿表格结果时设 `--alignment-output none` 可跳过可视化，更快。**

## 分析流程 | Pipeline

~~~text
查询序列(-i/-s) + 目标库(-r)
    │
    ▼
步骤1: 建库 makeblastdb → 01_database/{ref}.db
步骤2: 逐样本 BLAST 比对 → 02_blast/{sample}_{type}_results.tsv (16列, blastn带-task与-perc_identity)
步骤3: 合并全部样本 → blast_summary_results.tsv (18列, 补算查询覆盖度, 过滤<min_identity)
步骤4: 统计 + 排序 → blast_statistics.txt / blast_summary_results_sorted.tsv
步骤5: 高质量筛选 → blast_summary_results_sorted_high_quality.tsv
步骤6: 导出 Excel(软依赖) + 比对可视化(可选)
~~~

## 输出 | Output

~~~text
blast_output/
├── 00_pipeline_info/
│   └── software_versions.yml                 # 软件版本与参数记录
├── 01_database/
│   └── {ref}.db(.nhr/.phr 等)                # makeblastdb 建的目标库
├── 02_blast/
│   ├── {sample}_{type}_results.tsv           # 每样本原始比对(16列, 无表头, 末列qlen)
│   ├── blast_summary_results.tsv             # 合并汇总(18列, 含表头+查询覆盖度) —— 核心
│   ├── blast_summary_results_sorted.tsv      # 按覆盖度/相似度/E-value 排序
│   ├── blast_summary_results_sorted_high_quality.tsv  # 高质量命中
│   ├── blast_statistics.txt                  # 统计报告
│   └── blast_results.xlsx                    # 多sheet Excel(软依赖, 缺失则跳过)
├── 03_alignments/
│   ├── html/blast_alignments.html            # HTML 可视化(合并单文件, 默认; 概览+每样品tab)
│   └── text/all_samples_alignments.txt       # 文本比对摘要
└── 99_logs/                                  # 运行日志
~~~

> 旧模式（`--no-merge-html`）下 `html/` 目录为 `index.html` + 各样品 `{sample}_alignments.html`。合并文件是自包含的（样式/脚本/序列全部内联），分享给他人只发这一个文件即可。

- `blast_summary_results.tsv` 列（18 列）：样品名称 / 查询序列ID / 目标序列ID / 序列相似度(%) / 比对长度 / 错配数 / Gap数 / 查询起止 / 目标起止 / E-value / Bit_Score / 目标序列长度 / **查询覆盖度(%)** / 查询序列 / 目标序列 / 查询序列长度；低于 `--min-identity` 的行在合并阶段即被过滤

## 结果解读 | Interpreting Results

**通俗理解|In plain words:** 最常用 `blast_summary_results_sorted.tsv`（已按好坏排序）或 `..._high_quality.tsv`（已过滤）。看"序列相似度(%)"和"查询覆盖度(%)"两列：都高 = 几乎整条序列高度同源；相似度高但覆盖度低 = 只有一小段像（可能是结构域级别的保守）。

- 相似度 ≥ 90% 且查询覆盖度 ≥ 80%：基本可判为同源/同一基因
- 相似度 60–80%：分歧较大的同源（跨种/亚种常见），`--task blastn`（默认）才能报出来
- E-value 极小（如 ≤ 1e-20）：比对极显著，几乎不可能随机
- `blast_statistics.txt`：含总比对数、样品数、唯一查询/目标数、以及 E-value/相似度/覆盖度的分布区间计数，快速看整体
- `blast_results.xlsx`：`raw_results / summary / sorted / high_quality` 四个 sheet，适合在 Excel 里筛选

## 参数选择建议 | Parameter Guidance

- **找同源基因（常规）**：默认自动检测即可，`-i sequences/ -r genes.fa`
- **跨物种/品种间分歧比较**：默认 `--task blastn` 已够敏感能报 60%+ 的同源；追求速度且只关心高相似命中时换 `--task megablast`
- **只想要顶级命中**：收紧 `--high-quality-evalue 1e-30 --min-identity 90 --min-coverage 80`
- **想看全部原始比对不过滤**：`--min-identity 0`（注意：可视化默认每 query 只画 5 条，调 `--alignment-max-per-query`）
- **核酸 vs 蛋白不确定**：交给自动检测；确认时用 `--blast-type` 显式指定
- **只要表格、不要图**：`--alignment-output none`
- **结果要分享给同事**：默认产出的 `03_alignments/html/blast_alignments.html` 就是单个自包含文件，直接发送即可；需要旧的多文件结构时加 `--no-merge-html`
- **重跑并覆盖旧结果**：加 `--force`

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input, -i` | — | Path | 输入文件或目录路径｜Input file or directory path |
| `--sample-map-file, -s` | — | Path | 样品映射文件｜Sample mapping file |
| `--reference, -r` | 必填 | Path | 目标基因序列文件｜Target gene sequence file |
| `--output, -o` | `./blast_output` | Path | 输出目录路径｜Output directory |
| `--blast-type` | — | blastn/blastp/blastx/tblastn/tblastx | BLAST程序类型，默认根据输入文件自动检测｜BLAST program type (auto-detect from input files if not specified) |
| `--evalue, -e` | `1e-05` | float | E-value阈值｜E-value threshold |
| `--max-target-seqs` | `10` | int | 最大目标序列数｜Maximum target sequences |
| `--word-size` | — | int | 词大小，默认根据blast-type自动设置｜Word size (auto-set by blast-type if not specified) |
| `--threads, -t` | `12` | int | 线程数｜Number of threads |
| `--input-suffix` | `*.fa` | str | 输入文件后缀模式｜Input file suffix pattern |
| `--target-db-type` | — | nucl/prot | 目标数据库类型，默认根据blast-type自动设置｜Target database type (auto-set by blast-type if not specified) |
| `--task` | `blastn` | blastn/megablast/dc-megablast | blastn搜索任务，默认blastn(最敏感，可报出<70%分歧比对)｜blastn task, default blastn (most sensitive, reports <70% diverged hits) |
| `--min-identity` | `30.0` | float | 最小序列相似度｜Minimum sequence identity |
| `--min-coverage` | `50.0` | float | 最小覆盖度｜Minimum coverage |
| `--high-quality-evalue` | `1e-10` | float | 高质量比对E-value阈值｜High quality alignment E-value threshold |
| `--auto-detect-samples/--no-auto-detect-samples` | `True` |  | 自动检测样品名称｜Auto-detect sample names |
| `--sample-name-pattern` | `([^/]+?)(?:\.fa|\.fasta|\.fna)?$` | str | 样品名提取正则表达式｜Sample name extraction regex |
| `--sample-name` | — | str | 单文件输入时的样品名称｜Sample name for single-file input |
| `--makeblastdb-path` | `makeblastdb` | str | makeblastdb程序路径｜makeblastdb program path |
| `--blastn-path` | `blastn` | str | blastn程序路径｜blastn program path |
| `--blastp-path` | `blastp` | str | blastp程序路径｜blastp program path |
| `--blastx-path` | `blastx` | str | blastx程序路径｜blastx program path |
| `--tblastn-path` | `tblastn` | str | tblastn程序路径｜tblastn program path |
| `--tblastx-path` | `tblastx` | str | tblastx程序路径｜tblastx program path |
| `--alignment-output` | `both` | none/text/html/both | 比对可视化输出格式｜Alignment visualization output format |
| `--alignment-width` | `80` | int | 比对每行显示的字符数｜Characters per line in alignment display |
| `--alignment-min-identity` | `0.0` | float | 比对可视化最小相似度｜Minimum identity for alignment visualization |
| `--alignment-min-coverage` | `0.0` | float | 比对可视化最小覆盖度｜Minimum coverage for alignment visualization |
| `--alignment-max-per-sample` | `2000` | int | 每个样品最多显示的比对数｜Maximum alignments to display per sample |
| `--alignment-max-per-query` | `5` | int | 每条查询序列最多显示的比对数(按相似度取前N，保证所有query都有展示)｜Maximum alignments per query (top-N by identity, so every query is shown) |
| `--html-theme` | `modern` | modern/classic/dark | HTML主题样式｜HTML theme style |
| `--merge-html/--no-merge-html` | `True` |  | HTML输出合并为单个文件(默认)｜Merge HTML output into a single file (default) |
| `--verbose, -v` | — |  | 详细输出模式｜Verbose output mode |
| `--quiet` | — |  | 静默模式｜Quiet mode |
| `--log-level` | `INFO` | DEBUG/INFO/WARNING/ERROR/CRITICAL | 日志级别｜Log level |
| `--log-file` | — | Path | 日志文件路径｜Log file path |
| `--force, -f` | — |  | 强制覆盖已存在的结果｜Force overwrite existing results |
| `--dry-run` | — |  | 测试模式｜Test mode |
| `--version, -V` | — |  | 显示版本信息｜Show version information |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | — |  | 输入文件或目录路径｜Input file or directory path |
| `-s, --sample-map-file` | — |  | 样品映射文件(与-i二选一)｜Sample mapping file (alternative to -i) |
| `-r, --reference` | 必填 |  | 目标基因序列文件｜Target gene sequence file |
| `-o, --output-dir` | `./blast_output` |  | 输出目录｜Output directory |
| `--blast-type` | — | blastn/blastp/blastx/tblastn/tblastx | BLAST程序类型,默认自动检测｜BLAST program type (auto-detect if not specified) |
| `-e, --evalue` | `1e-05` | float | E-value阈值｜E-value threshold |
| `--max-target-seqs` | `10` | int | 最大目标序列数｜Maximum target sequences |
| `--word-size` | — | int | 词大小,默认按blast-type设置｜Word size (auto-set by blast-type) |
| `-t, --threads` | `12` | int | 线程数｜Number of threads |
| `--input-suffix` | `*.fa` |  | 输入文件后缀模式｜Input file suffix pattern |
| `--target-db-type` | — | nucl/prot | 目标数据库类型,默认按blast-type设置｜Target database type (auto-set by blast-type) |
| `--task` | `blastn` | blastn/megablast/dc-megablast | blastn搜索任务,默认blastn(最敏感,可报出<70%%分歧比对;megablast最快但仅适合高相似序列)｜blastn task, default blastn (most sensitive, reports <70%% diverged hits; megablast is fast but only for high-identity searches) |
| `--min-identity` | `30.0` | float | 最小序列相似度(%%),传给blastn -perc_identity并过滤汇总｜Minimum sequence identity (%%), passed to blastn -perc_identity and applied to summary |
| `--min-coverage` | `50.0` | float | 最小覆盖度(%%)｜Minimum coverage (%%) |
| `--high-quality-evalue` | `1e-10` | float | 高质量比对E-value阈值｜High quality alignment E-value threshold |
| `--sample-name` | — |  | 单文件输入时的样品名称｜Sample name for single-file input |
| `--no-auto-detect-samples` | — | store_false | 关闭自动检测样品名称｜Disable auto sample name detection |
| `--sample-name-pattern` | `([^/]+?)(?:\.fa|\.fasta|\.fna)?$` |  | 样品名提取正则表达式｜Sample name extraction regex |
| `--makeblastdb-path` | — |  | makeblastdb程序路径｜makeblastdb program path |
| `--blastn-path` | — |  | blastn程序路径｜blastn program path |
| `--blastp-path` | — |  | blastp程序路径｜blastp program path |
| `--blastx-path` | — |  | blastx程序路径｜blastx program path |
| `--tblastn-path` | — |  | tblastn程序路径｜tblastn program path |
| `--tblastx-path` | — |  | tblastx程序路径｜tblastx program path |
| `--alignment-output` | `both` | none/text/html/both | 比对可视化输出格式｜Alignment visualization output format |
| `--alignment-width` | `80` | int | 比对每行显示的字符数｜Characters per line in alignment display |
| `--alignment-min-identity` | `0.0` | float | 比对可视化最小相似度过滤｜Minimum identity for alignment visualization |
| `--alignment-min-coverage` | `0.0` | float | 比对可视化最小覆盖度过滤｜Minimum coverage for alignment visualization |
| `--alignment-max-per-sample` | `2000` | int | 每个样品最多显示的比对数｜Maximum alignments to display per sample |
| `--alignment-max-per-query` | `5` | int | 每条查询序列最多显示的比对数(按相似度取前N,保证所有query都有展示)｜Maximum alignments per query (top-N by identity, so every query is shown) |
| `--html-theme` | `modern` | modern/classic/dark | HTML主题样式｜HTML theme style |
| `--merge-html` | `True` | store_true | HTML输出合并为单个文件(默认)｜Merge HTML output into a single file (default) |
| `--no-merge-html` | — | store_false | 关闭HTML合并,输出index和分样品多文件｜Disable merging (legacy multi-file output) |
| `-v, --verbose` | `0` | count | 详细输出(-vv更详细)｜Verbose output (-vv for more) |
| `--quiet` | — | store_true | 静默模式｜Quiet mode |
| `--log-level` | `INFO` | DEBUG/INFO/WARNING/ERROR/CRITICAL | 日志级别｜Log level |
| `--log-file` | — |  | 日志文件路径｜Log file path |
| `-f, --force` | — | store_true | 强制覆盖已存在文件｜Force overwrite existing files |
| `--dry-run` | — | store_true | 模拟运行不执行｜Dry run without execution |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- BLAST+ 套件：`makeblastdb` + `blastn/blastp/blastx/tblastn/tblastx`（conda 环境默认 `Blast_v.2.16.0`，路径可经 `--*-path` 参数或环境变量覆盖）
- `Biopython`（`Bio.SeqIO`，用于自动检测序列类型）
- 软依赖：`pandas` + `openpyxl`（仅 Excel 输出，缺失会自动跳过，TSV 照常输出）

## 常见问题 | FAQ

**Q1：会不会断点续传？**
会。已建好的目标库（`.db.nhr`/`.db.phr` 存在）和已比对且非空的样品结果默认跳过；加 `--force` 可强制全部重跑。换参数重跑前建议先 `--force` 或删旧产物，否则会复用旧结果。

**Q2：为什么结果文件是空的？**
样品与目标库"零命中"是合法结果（rc=0 但无命中），会保留空文件并继续，不算失败。可放宽 `-e` 或检查序列类型是否匹配（核酸库用 blastp 当然比不到）。

**Q3：blast-type 是怎么自动判断的？**
读查询与参考的前几条序列，按碱基/氨基酸字符占比判断核酸还是蛋白，再按组合选程序（如 查询核酸+参考蛋白 → blastx）。只有同时给了 `-i` 和 `-r` 才自动检测；只给 `-s` 映射时默认 `blastn`，请用 `--blast-type` 显式指定。

**Q4：覆盖度列是怎么算的？**
是"查询序列被覆盖"的比例（`abs(qend-qstart+1)/qlen ×100`，上限 100），衡量整条查询序列对上了多少，和相似度是两个维度。v2.3.0 之前按目标序列长度算，对染色体级目标库恒 ≈0%，已废弃该口径。

**Q7：为什么结果里最小相似度总是卡在 70% 附近？怎么降到更低？**
v2.3.0 之前默认用 megablast 任务搜索，它只适合高相似序列，实测报不出 <70% 的比对（调 E-value、`--min-identity` 都没用，因为比对在 BLAST 内部就没被报告）。现在默认 `--task blastn`（最敏感，可报出 60% 左右的分歧同源）；旧结果如需补齐低相似命中，删除 `02_blast/` 或加 `--force` 重跑。

**Q8：改了参数重跑，结果怎么没变化？**
断点续传会复用旧的建库和比对结果（换参数不影响"输出文件已存在"的判断）。参数有变时加 `--force`，或先删除 `01_database/ 02_blast/`；v2.2.0 及更早的旧 raw 结果是 15 列（无 qlen 列），新版本无法直接复用，必须 `--force` 重跑。

**Q5：Excel 没生成？**
Excel 是软依赖，需要 `pandas` + `openpyxl`；没装会自动跳过（TSV 结果不受影响）。装了仍失败可能是某 sheet 超 Excel 行数上限，会被跳过并提示。

**Q6：之前跑的结果是 index.html + 分样品两个文件，怎么重新合并成单文件？**
不用重跑 BLAST：用相同参数（同一个 `-o`）重新运行一次即可。断点续传会跳过建库和比对步骤，只重新生成可视化，默认即产出合并单文件 `blast_alignments.html`；若想保留旧的多文件格式，加 `--no-merge-html`。
