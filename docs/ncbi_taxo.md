# NCBI 分类学注释 | NCBI Taxonomy Annotation

一句话理解：**给你的一堆序列比对结果（BLAST）或 accession 编号，查出它们各属于哪个物种、门、纲、目，并统计「谁家占了多少」**。

输入一个 BLAST 结果（制表符）或纯 accession 列表，输出每个 accession 的 TaxID 和完整分类谱系（界/门/纲/目/科/属/种），以及按层级的统计表。

## 功能概述 | Overview { #overview }

- 从 BLAST 结果（或 accession 列表）提取唯一 accession，批量查 NCBI TaxID 与分类谱系
- 用 taxonkit + accession2taxid 数据库本地查询，无需联网即可拿到界/门/纲/目/科/属/种
- 对 BLAST 输入可按「比对长度」过滤（默认只保留 >=1000 bp 的命中）
- 统计口径可选：按唯一 accession 数、按 BLAST hit 次数、或两者都算；层级可自选
- 可选功能：联网获取 accession 的序列描述并归类（线粒体/叶绿体/rRNA 等）
- 断点续传：无（每次从头跑）

## 快速开始 | Quick Start { #quick-start }

`@bash
biopytools ncbi-taxo -i blast_results.txt -o taxo_out
`@

最小输入：一个 BLAST 结果文件（制表符分隔）或 accession 列表，`-o` 指定输出文件前缀。

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语 | 通俗理解 |
|------|----------|
| accession | 序列在 NCBI 的「身份证号」，形如 NP_123456 或 XM_001234 |
| TaxID | 物种在 NCBI 分类系统里的「编号」，一个数字 |
| lineage（谱系） | 一个物种从界到种的「完整户口链」，如 界;门;纲;目;科;属;种 |
| BLAST 结果 | 序列比对软件的输出表，每行一条「命中」记录 |
| hit 次数 | 一个 accession 在 BLAST 结果里出现了几次（出现多=比对到它多次） |

## 输入 | Input { #input }

两种输入，程序默认自动识别（也可用 `--input-type` 手动指定）：

### BLAST 结果（默认识别）

制表符分隔的表，程序默认取**第 2 列**当 accession（可用 `--blast-column` 改），取**第 4 列**当比对长度做过滤（长度 <1000 bp 的行被丢弃，可用 `-l` 调阈值）。

### accession 列表

每行一个 accession 编号的纯文本。

`@text
NP_001234
XM_002345
XP_003456
`@

## 参数说明 | Parameters { #parameters }

### 必需参数 | Required

**通俗理解|In plain words:** 输入文件 + 输出前缀。输出前缀决定所有结果文件名（见「输出」）。

### 输入解析 | Input parsing

**通俗理解|In plain words:** `--input-type` 一般保持 auto 自动识别即可；`--blast-column` 管「accession 在第几列」（默认第 2 列，1 起算），只有你的 BLAST 表列顺序不同时才改；`-l/--min-length` 管「多短的命中算噪声丢掉」（默认 1000 bp，调小=保留更多短命中，调大=更严格）。

### 分类学数据库 | Taxonomy database

**通俗理解|In plain words:** `--taxid-db` 指向本地 accession2taxid 数据库文件（默认 `~/database/ncbi_taxonomy/nucl_gb.accession2taxid.gz`）。**第一次用之前必须先下载好这个库**（见 FAQ），否则程序会直接报错退出。`--lineage-format` 与 `--no-full-lineage` 控制谱系串的展示格式，一般不用动。

### 统计配置 | Statistics

**通俗理解|In plain words:** `--stats-by` 选「按哪些层级统计」（默认属 + 种，可选界到种任意组合）；`--stats-target` 选「按什么口径计数」——blast_hits 按命中次数、unique_accessions 按去重后的编号数、both 两者都给（默认）；`--stats-output` 选统计表存成 txt 还是 csv。

### 序列描述与性能 | Titles & performance

**通俗理解|In plain words:** `--fetch-titles` 会联网去 NCBI 抓每个 accession 的描述并归类（需网络 + edirect/Biopython），默认关闭；`-t` 线程数默认 4，一般不用动。

## 分析流程 | Pipeline { #pipeline }

`@text
步骤1: 提取 accession（BLAST 按列+长度过滤 / 列表直接读）
   |
   v
步骤2: zgrep 在 accession2taxid 库里查 TaxID
   |
   v
步骤3: taxonkit 查询每个 TaxID 的完整 lineage
   |
   v
步骤4(可选): --fetch-titles 联网抓序列描述并归类
   |
   v
步骤5: 按指定层级/口径统计 -> 写统计文件 + 汇总报告
`@

## 输出 | Output { #output }

`@text
<前缀>.accessions.txt       # 提取到的唯一 accession 列表
<前缀>.acc2taxid.txt        # accession -> TaxID 映射（查不到的 TaxID 为空）
<前缀>.taxonomy.txt         # accession / TaxID / 完整 lineage（三列）
<前缀>.statistics.txt       # 统计表（txt；--stats-output csv 时为 .statistics.csv）
<前缀>.accession2title.txt  # accession -> 序列描述（仅 --fetch-titles）
<前缀>.sequence_types.txt   # 序列类型归类统计（仅 --fetch-titles）
<前缀>.ncbi_taxo.log        # 运行日志
`@

## 结果解读 | Interpreting Results { #interpreting }

- **`<前缀>.taxonomy.txt`**：最核心结果，三列——accession、TaxID、完整 lineage。用「科/属/种」列对应回你关心的物种即可。
- **`<前缀>.statistics.txt`**：按层级的统计表，含每个分类名的数量与百分比。看「哪个属/种占大头」——占比高的就是你的数据里最主要的类群。
- **`<前缀>.acc2taxid.txt`**：查 TaxID 的中间表。**若某行 TaxID 为空，说明该 accession 没在数据库里查到**（可能是新提交序列、或库太旧），需要更新 `--taxid-db`。
- **好坏判据**：统计里「unclassified」或空 TaxID 比例越低，注释越完整；若大片查不到，先检查数据库是否最新、accession 列是否取错（比如取成了 query 列）。

## 参数选择建议 | Parameter Guidance { #guidance }

- **常规注释 BLAST 结果**：只给 `-i -o`，自动识别 + 默认过滤即可。
- **BLAST 表 accession 不在第 2 列**：加 `--blast-column N`（N 从 1 起算）。
- **想保留短比对（如 ITS/16S 片段）**：调小 `-l`（如 `-l 100`）。
- **想按科/目统计而非属/种**：`--stats-by order family`。
- **想顺便知道序列是什么类型（线粒体/叶绿体/rRNA）**：加 `--fetch-titles`（需联网）。

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input, -i` | 必填 |  | 输入文件（BLAST结果或accession列表）｜Input file (BLAST results or accession list) |
| `--output-prefix, -o` | 必填 |  | 输出文件前缀｜Output file prefix |
| `--input-type` | `auto` | auto/blast/accession | 输入文件类型｜Input file type |
| `--blast-column` | `2` | int | BLAST结果中accession所在的列（从1开始）｜Column containing accession in BLAST (1-based) |
| `--min-length, -l` | `1000` | int | 最小比对长度过滤（bp）｜Minimum alignment length filter (bp) |
| `--fetch-titles` | — |  | 获取accession的序列描述（需要edirect）｜Fetch accession titles (requires edirect) |
| `--taxid-db` | `~/database/ncbi_taxonomy/nucl_gb.accession2taxid.gz` |  | TaxID数据库路径｜TaxID database path |
| `--lineage-format` | `{k};{p};{c};{o};{f};{g};{s}` |  | 分类层级格式｜Lineage format string |
| `--no-full-lineage` | — |  | 不保留完整lineage｜Do not keep full lineage |
| `--stats-by` | `['genus', 'species']` | kingdom/phylum/class/order/family/genus/species | 统计层级（可多选）｜Statistics levels (can be specified multiple times) |
| `--stats-target` | `both` | blast_hits/unique_accessions/both | 统计对象｜Statistics target |
| `--stats-output` | `txt` | txt/csv | 统计输出格式｜Statistics output format |
| `--threads, -t` | `4` | int | 线程数｜Number of threads |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 输入文件路径（BLAST结果或accession列表）｜Input file path (BLAST results or accession list) |
| `-o, --output-prefix` | 必填 |  | 输出文件前缀｜Output file prefix |
| `--input-type` | `auto` | auto/blast/accession | 输入文件类型｜Input file type (auto/blast/accession) |
| `--blast-column` | `2` | int | BLAST结果中accession所在的列（从1开始）｜Column containing accession in BLAST results (1-based) |
| `-l, --min-length` | `1000` | int | 最小比对长度过滤（bp）｜Minimum alignment length filter (bp) |
| `--fetch-titles` | — | store_true | 获取accession的序列描述｜Fetch accession titles (requires edirect) |
| `--taxid-db` | `~/database/ncbi_taxonomy/nucl_gb.accession2taxid.gz` |  | TaxID数据库路径｜TaxID database path |
| `--lineage-format` | `{k};{p};{c};{o};{f};{g};{s}` |  | 分类层级格式｜Lineage format string |
| `--no-full-lineage` | — | store_true | 不保留完整lineage｜Do not keep full lineage |
| `--stats-by` | `['genus', 'species']` | kingdom/phylum/class/order/family/genus/species | 统计层级｜Statistics levels |
| `--stats-target` | `both` | blast_hits/unique_accessions/both | 统计对象｜Statistics target |
| `--stats-output` | `txt` | txt/csv | 统计输出格式｜Statistics output format |
| `-t, --threads` | `4` | int | 线程数｜Number of threads |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies { #dependencies }

- taxonkit（须在 PATH 中，conda 可装：`conda install -c bioconda taxonkit`）
- zgrep、zcat（系统自带，用于读压缩数据库）
- 本地 TaxID 数据库（默认 `~/database/ncbi_taxonomy/nucl_gb.accession2taxid.gz`，需提前下载）
- 可选：edirect + BioPython（仅 `--fetch-titles` 联网抓描述时需要）

## 常见问题 | FAQ { #faq }

**Q1：报「TaxID database not found」怎么办？**
`--taxid-db` 指向的数据库文件不存在。需要先下载 NCBI 的 `nucl_gb.accession2taxid.gz`（或 prot 版）放到对应路径，或用 `--taxid-db` 指向你已有的库。

**Q2：为什么很多 accession 查不到 TaxID？**
常见原因：数据库太旧（新提交的序列不在里面）；或 accession 列取错了（`--blast-column` 没指对，比如把 query 名当成了 accession）。先核对列，再考虑更新数据库。

**Q3：BLAST 结果里哪一列被当成了比对长度？**
程序固定取**第 4 列**当比对长度（标准 BLAST 出表格格式：query、subject、identity、alignment length……）。如果你的表不是这个格式，过滤逻辑可能失效，建议把 accession 和长度列整理成标准顺序。

**Q4：`--fetch-titles` 需要什么条件？**
需要能访问 NCBI（联网）、装了 Biopython（`pip install biopython`）。它会按 NCBI 速率限制分批查询，大批量时较慢；没网就跳过这个选项。

**Q5：统计里 blast_hits 和 unique_accessions 有什么区别？**
unique_accessions 是「去重后有多少个不同 accession 属于某类」；blast_hits 是「这些 accession 总共被命中多少次」。后者能反映「测序深度/丰度」倾向，默认 both 两个都出。
