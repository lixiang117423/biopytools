# eggnog-mapper - eggNOG 功能注释 | eggNOG functional annotation

一句话理解：**给一批蛋白(或基因)序列，逐个查「它在别的物种里叫什么、属于哪个功能家族、参与什么通路」，输出 GO、KEGG、COG、CAZy、Pfam 等一大堆功能标签**。
它解决的问题：基因预测完了，但还不知道每个基因「是干什么的」，eggnog-mapper 用海量直系同源数据库一次性给你标上通用的功能注释。

## 功能概述 | Overview { #overview }

- 基于 eggNOG 直系同源数据库，给蛋白/CDS/基因组序列做功能注释
- 一次输出多类注释：GO、KEGG(通路/模块/反应)、COG 分类、CAZy(碳水化合物酶)、Pfam 结构域、EC 酶号等
- 支持多种搜索引擎：mmseqs(默认，快)/ diamond / hmmer，以及 no_search、cache 模式
- 原生结果自动重排版为**中英双语表头**的 TSV 和 Excel(装了 openpyxl 才有 Excel)
- 支持断点续传(`--resume` + 产物存在性判断)与覆盖重跑(`--override`)

## 快速开始 | Quick Start { #quick-start }

```bash
biopytools eggnog-mapper -i proteins.faa -o out/
```

最小输入：一个蛋白 FASTA。输出目录会自动生成 `00_pipeline_info` / `01_emapper` / `99_logs` 三个子目录。

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语 | 通俗理解 |
|------|----------|
| 功能注释 | 给每个基因/蛋白贴「它是干什么的」标签 |
| 直系同源 | 不同物种里「同一个祖先基因的后代」，功能通常保守，可用已知的推未知 |
| GO | 基因本体论，用「分子功能/生物学过程/细胞定位」三套标准词条描述基因 |
| KEGG | 代谢与信号通路数据库，告诉你基因在哪些通路上 |
| COG | 按功能分类的直系同源群(如「能量代谢」「转录」) |
| CAZy | 碳水化合物活性酶数据库(降解糖/纤维素的酶) |
| Pfam | 蛋白结构域数据库，按保守的「功能模块」划分 |
| EC 酶号 | 酶的编号系统，一个号对应一种催化反应 |
| 搜索引擎 | 把你的序列拿去数据库「找最像的」，mmseqs 快、diamond 通用、hmmer 更灵敏 |

## 输入 | Input { #input }

### 输入 FASTA(-i，必填)

- 蛋白序列(`.faa` / `.fa` / `.fasta`)：默认输入类型
- CDS(`--itype CDS`)：可加 `--translate` 先翻译成蛋白再注释
- 基因组(`--itype genome`)或宏基因组(`--itype metagenome`)：同样配 `--translate`

### 数据库(需预下载)

eggnog-mapper 依赖本地 eggNOG 数据库(默认 `~/database/eggnog`，可用 `--data-dir` 改)，需包含 `eggnog.db`、`eggnog.taxa.db`，以及对应模式的搜索库(mmseqs 或 diamond)。未下载时程序会报错并提示下载命令。

## 参数说明 | Parameters { #parameters }

### 必需参数 | Required

**通俗理解|In plain words:** 输入 FASTA + 输出目录。输出前缀默认取输入文件名(可用 `--prefix` 改)。

### 输入类型 | Input type

**通俗理解|In plain words:** `--itype` 告诉程序你给的是什么。默认 `proteins` 直接注释；给 CDS/基因组/宏基因组时选对应类型，并加 `--translate` 让它先翻译成蛋白。**给蛋白时不要加 `--translate`(会被忽略并警告)**。

### 搜索模式与引擎 | Search mode & engine

**通俗理解|In plain words:** `-m` 选搜索引擎。`mmseqs`(默认)最快，适合绝大多数情况；`diamond` 是 BLAST 的加速替代；`hmmer` 用结构域模型更灵敏但慢得多；`no_search`/`cache` 用于已有比对结果的特殊场景。**一般用默认 mmseqs 即可**。`--sensmode`(默认 sensitive)是搜索灵敏度，`--seed-ortholog-evalue`(默认 0.001)是判定「像不像」的阈值，**一般不用动**。

### 运行与续传 | Runtime & resume

**通俗理解|In plain words:** `--cpu` 线程数(默认 12)。`--resume` 断点续传、`--override` 覆盖已有结果重跑、`--no-format` 跳过重排版只留原生文件。默认已有 `.emapper.annotations` 产物时自动跳过搜索(等价于隐式续传)。

### 路径覆盖 | Path overrides

**通俗理解|In plain words:** `--data-dir`(数据库目录，默认 `~/database/eggnog`)、`--emapper-path`(emapper.py 路径，默认 `~/miniforge3/envs/annot/bin/emapper.py`)。**一般不用动**，数据库装在别处时才需要。

## 分析流程 | Pipeline { #pipeline }

```text
输入: FASTA(蛋白/CDS/基因组) + eggNOG 数据库
    |
    v
步骤1: 搜索(断点续传, 已有 .emapper.annotations 则跳过)
  - 用 mmseqs/diamond/hmmer 把序列比对到 eggNOG 数据库
  - [--translate] 先把 CDS 翻译成蛋白
  - 生成 .emapper.annotations(原生注释) + seed_orthologs + hits
    |
    v
步骤2: 重排版(可 --no-format 跳过)
  - 解析原生注释 -> 中英双语表头 TSV(.cn.tsv)
  - 生成 Excel(.xlsx, 需 openpyxl)
    |
    v
步骤3: 写版本信息 -> 完成
```

## 输出 | Output { #output }

```text
out/
├── 00_pipeline_info/
│   └── software_versions.yml             # 软件版本与参数
├── 01_emapper/
│   ├── <prefix>.emapper.annotations      # 原生注释结果(制表符分隔)
│   ├── <prefix>.emapper.annotations.cn.tsv  # 中英双语表头 TSV(推荐直接用这个)
│   ├── <prefix>.emapper.annotations.xlsx    # Excel 版(装了 openpyxl 才有)
│   ├── <prefix>.emapper.seed_orthologs   # 种子直系同源命中
│   └── <prefix>.emapper.hits             # 比对命中详情
└── 99_logs/
    └── emapper.log                        # 运行日志
```

## 结果解读 | Interpreting Results { #results }

### 1. 注释结果(`01_emapper/<prefix>.emapper.annotations.cn.tsv`)

**通俗理解|In plain words:** 这是核心结果——每行一个输入序列，各列是它被标注上的功能信息。`.cn.tsv` 版把列名换成了中英双语，更适合人读。

关键列含义：

- `query`：你的序列 ID
- `seed_ortholog`：最相似的那个已知直系同源蛋白(证据来源)
- `evalue` / `score`：比对的可靠程度(evalue 越小越可信)
- `Preferred_name`：推荐基因名(如 `rpoB`)，最直观
- `Description`：功能一句话描述
- `GOs`：GO 词条；`KEGG_ko` / `KEGG_Pathway`：KEGG 通路；`EC_number`：酶号
- `CAZy`：碳水化合物酶；`PFAMs`：Pfam 结构域；`COG_category`：COG 功能类

- 有注释的序列比例高(如 >70%)说明注释成功、物种不算太偏门
- evalue 高(如 >1e-3)、注释稀疏的条目可信度低

### 2. 日志统计

**通俗理解|In plain words:** 日志末尾会打印「总序列数 / 有注释条数」，快速看注释覆盖率。

## 参数选择建议 | Parameter Guidance { #guidance }

- **常规蛋白注释**：`-i proteins.faa -o out/`，其余全默认
- **注释的是 CDS**：`-i genes.cds.fa --itype CDS --translate`
- **速度不够**：默认 mmseqs 已最快；确认用多线程 `--cpu 24`
- **要 Excel 给同事看**：确认装了 openpyxl(否则只有 TSV)，或直接用 `.cn.tsv` 在 Excel 里打开
- **数据库不在默认位置**：`--data-dir /path/to/eggnog`

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 输入FASTA(蛋白/CDS/基因组)｜Input FASTA |
| `-o, --output-dir` | 必填 | Path | 输出目录｜Output directory |
| `--itype` | `proteins` | proteins/CDS/genome/metagenome | 输入类型｜Input type |
| `--translate` | — |  | CDS翻译为蛋白｜Translate CDS |
| `-m, --mode` | `mmseqs` | mmseqs/diamond/hmmer/no_search/cache | 搜索模式｜Search mode |
| `--cpu` | `12` |  | 线程数｜Threads |
| `--sensmode` | `sensitive` |  | 灵敏度｜Sensitivity |
| `--seed-ortholog-evalue` | `0.001` |  | seed ortholog E值｜evalue |
| `--data-dir` | — |  | DB目录｜DB directory (default: ~/database/eggnog) |
| `--prefix` | — |  | 输出前缀｜Output prefix |
| `--emapper-path` | — |  | emapper.py路径｜emapper.py path override |
| `--resume` | — |  | 续传｜Resume |
| `--override` | — |  | 覆盖｜Override existing output |
| `--no-format` | — |  | 跳过重排版｜Skip reformat |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 输入FASTA(蛋白/CDS/基因组)｜Input FASTA |
| `-o, --output-dir` | 必填 |  | 输出目录｜Output directory |
| `--itype` | `proteins` | proteins/CDS/genome/metagenome | 输入类型｜Input type (default: proteins) |
| `--translate` | — | store_true | CDS翻译为蛋白(itype=CDS/genome/metagenome)｜Translate CDS |
| `-m, --mode` | `mmseqs` | mmseqs/diamond/hmmer/no_search/cache | 搜索模式｜Search mode (default: mmseqs) |
| `--cpu` | `12` | int | 线程数｜Threads (default: 12) |
| `--sensmode` | `sensitive` |  | 灵敏度｜Sensitivity (default: sensitive) |
| `--seed-ortholog-evalue` | `0.001` | float | seed ortholog E值｜evalue (default: 0.001) |
| `--data-dir` | — |  | DB目录｜DB directory (default: ~/database/eggnog) |
| `--prefix` | — |  | 输出前缀｜Output prefix (default: input stem) |
| `--emapper-path` | — |  | emapper.py路径｜emapper.py path override |
| `--resume` | — | store_true | 续传｜Resume |
| `--override` | — | store_true | 覆盖｜Override existing output |
| `--no-format` | — | store_true | 跳过重排版,只留原生产物｜Skip reformat |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies { #dependencies }

- **eggnog-mapper(emapper.py)**：conda 环境 `annot`(默认 `~/miniforge3/envs/annot/bin/emapper.py`，可用 `--emapper-path` 或环境变量 `EMAPPER_PATH` 覆盖)
- **diamond**：仅 `-m diamond` 模式需要
- **eggNOG 数据库**：默认 `~/database/eggnog`，需含 `eggnog.db` + `eggnog.taxa.db` + 对应搜索库(mmseqs/diamond)
- **openpyxl**：可选，生成 Excel 版注释时用

## 常见问题 | FAQ { #faq }

**Q1：中断后重跑会从头再来吗？**
不会。已有 `.emapper.annotations` 产物时会自动跳过搜索(除非加 `--override`)。也可显式加 `--resume` 让 emapper 自己续传。

**Q2：报「缺少 mmseqs 库 / eggnog.db」怎么办？**
数据库没下全。默认数据库目录 `~/database/eggnog` 需含 `eggnog.db`、`eggnog.taxa.db` 和 mmseqs/diamond 搜索库。按报错提示用 `download_eggnog_data.py` 下载(程序会打印具体命令)。

**Q3：`--translate` 对蛋白输入有用吗？**
没用，会被忽略并警告。只有输入类型是 CDS/genome/metagenome 时才需要加 `--translate` 先翻译。

**Q4：mmseqs 和 diamond 模式怎么选？**
默认 mmseqs 最快，绝大多数情况够用。diamond 是 BLAST 类比对，也很快。hmmer 最灵敏但慢。无特殊需求用默认 mmseqs。

**Q5：`--no-format` 是干什么的？**
跳过重排版，只留 emapper 的原生 `.emapper.annotations`(列头是英文)。想要中英双语 TSV/Excel 就别加这个参数。

**Q6：注释覆盖率低正常吗？**
取决于物种。模式物种(人、拟南芥、酵母)覆盖率可达 90%+；偏门物种或垃圾序列覆盖率低属正常。重点看有注释的条目质量(evalue 小、有 GO/KEGG 词条)，而不是单纯追求数量。
