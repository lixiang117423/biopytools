# FAPROTAX - 微生物群落功能注释 | FAPROTAX Functional Annotation

一句话理解：**把「每种细菌/古菌是谁」的清单，换算成「这些微生物整体上会干什么」的功能清单**，用于回答「这个环境里谁在产甲烷、谁在固氮、谁在还原硫酸盐」这类生态问题。

## 功能概述 | Overview { #overview }

- 输入一张 OTU/ASV 丰度表（BIOM 或 TSV），按分类学注释查表归类为功能组
- 核心依赖 FAPROTAX 官方脚本 collapse_table.py，本模块只是封装并规范化输出
- 支持 7 种标准化方式和 7 种组内聚合方式（参数透传给 collapse_table.py）
- 未匹配到任何功能组的记录可归入自定义分组（--group-leftovers-as）
- 断点续传：检测到已完成的功能表即跳过重跑（换参数需先删旧产物）
- 输出统一归入编号目录，并自动生成 software_versions.yml

## 快速开始 | Quick Start { #quick-start }

```bash
biopytools faprotaxtax -i otu_table.biom -o faprotaxtax_output
```

最小输入：一张含分类学注释的 OTU/ASV 表。默认数据库与脚本路径为 `~/software/FAPROTAX/FAPROTAX_1.2.12/` 下的 FAPROTAX.txt 和 collapse_table.py。

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语 | 通俗理解 |
|------|----------|
| OTU / ASV | 把序列几乎相同的一批微生物归成的一个「操作单元」，像把相似商品归为同一个「货号」 |
| 分类学注释 | 每个 OTU/ASV 在生命树上的「住址」（界门纲目科属种），决定它「姓什么」 |
| 功能组 | 某类微生物「会干什么」的标签（产甲烷、固氮、硫酸盐还原等），是 FAPROTAX 的输出维度 |
| FAPROTAX | 一张「细菌职业对照表」：按属/种名查它通常具备哪些生态功能 |
| BIOM | 微生物群落表的一种通用机器格式，除数值外还带元数据（可存分类学信息） |
| 标准化(normalize) | 把不同样本的总量拉到同一水平再比，类似「看百分比」而不是「看绝对人数」 |
| 聚合(average) | 多个 OTU 归到同一功能组后，怎么合并它们的数值（求和/平均/取最大等） |
| 元数据字段 | BIOM 表里额外附着的信息列，本工具用它来定位「分类学」那一列 |

## 输入 | Input { #input }

支持两种格式，程序按扩展名自动识别（.biom/.hbiom 等按 BIOM 处理，其余按 TSV 处理）。

- **BIOM**（推荐）：需包含分类学元数据字段，默认配合 `--collapse-by-metadata taxonomy` 指定字段名。
- **TSV（经典表格）**：第一列或某一列是分类学信息，配 `--collapse-by-metadata`（指定列名）或 `--row-names-are-in-column`（指定行名列）。

TSV 示例：

```text
OTU_ID    taxonomy                              S1    S2
OTU1      k__Bacteria;p__Proteobacteria;...     12    0
OTU2      k__Bacteria;p__Firmicutes;...         8     15
```

## 参数说明 | Parameters { #parameters }

### 必需与输出 | Required & output

**通俗理解|In plain words:** `-i` 是输入的丰度表，`-o` 是输出目录（默认 ./faprotaxtax_output）。这两个就够了，其余都有默认值。

### 数据库与脚本路径 | Database & script paths

**通俗理解|In plain words:** 指定 FAPROTAX 的两件「家当」——功能组数据库（FAPROTAX.txt）和核心脚本（collapse_table.py）。**默认路径通常已配好，一般不用动**；只有当你的 FAPROTAX 装在别处时才需要 `-g` 和 `--collapse-table-path` 指过去。

### 注释匹配 | Annotation matching

**通俗理解|In plain words:** 告诉程序「分类学信息在表里的哪一列」。BIOM 表用 `--collapse-by-metadata taxonomy`（字段名），经典 TSV 用 `--collapse-by-metadata` 或 `--row-names-are-in-column`。没匹配到任何功能组的记录，可用 `--group-leftovers-as` 给个名字（如 other）收拢起来，否则会被丢弃。

### 标准化与聚合 | Normalization & aggregation

**通俗理解|In plain words:** `-n/--normalize` 决定「先按列还是按行、折叠前还是折叠后」把数值归一化；`--average` 决定多个 OTU 归到同一功能组后怎么合并数值。**绝大多数情况用默认 none 即可**（先看原始相对丰度），只有当你想做严格的样本间可比或要控制稀有物种影响时才启用。

### 输出格式与执行 | Output format & execution

**通俗理解|In plain words:** `--output-format` 默认 auto（BIOM 输入输出 .biom，TSV 输入输出 .tsv）；`-t` 线程数默认 1（collapse_table.py 是单线程脚本，多线程意义不大）；`-f` 覆盖已有输出，`-v` 打印详细过程。

## 分析流程 | Pipeline { #pipeline }

```text
输入 OTU/ASV 表(BIOM/TSV)
    │
    ▼
断点续传检查(01_collapsed/functional_table 已存在则跳过)
    │
    ▼
自动识别输入格式 + 校验数据库/脚本路径
    │
    ▼
构建并执行 collapse_table.py(含 --group_members_defined_as words)
    │
    ▼
校验输出(functional_table + 可选 report)
    │
    ▼
重组文件到编号目录 + 生成 software_versions.yml + 清理 work/
```

## 输出 | Output { #output }

```text
faprotaxtax_output/
├── 00_pipeline_info/
│   └── software_versions.yml          # 脚本版本/数据库信息/运行参数
├── 01_collapsed/
│   └── functional_table.tsv           # 功能丰度表(BIOM输入时为 .biom)
├── 02_report/
│   └── faprotaxtax_report.txt         # 折叠报告(仅指定 --group-leftovers-as 时)
└── 99_logs/
    └── faprotaxtax_pipeline.log       # 运行日志
```

`work/` 为临时目录，运行结束自动清理。

## 结果解读 | Interpreting Results { #interpreting-results }

### 1. functional_table（核心结果）

**通俗理解|In plain words:** 这是最终交付表——行是功能组（如 methane oxidation、nitrogen fixation），列是样本，数值是各功能组的丰度。直接拿去做下游统计/绘图。

### 2. faprotaxtax_report.txt

**通俗理解|In plain words:** collapse_table.py 的折叠过程报告，记录了「哪些分类学单元被合并进了哪个功能组」。用 `--group-leftovers-as` 时会生成，用于核对未匹配记录的归并情况。

### 3. software_versions.yml 与日志

**通俗理解|In plain words:** 版本与参数存档，便于论文 Methods 和可复现；日志里能看到完整执行命令和耗时。

## 参数选择建议 | Parameter Guidance { #parameter-guidance }

- **常规分析**：全部默认即可，只需 `-i` 和 `-o`。
- **TSV 输入**：务必加 `--collapse-by-metadata taxonomy`（列名按你的表头改）。
- **想保留未匹配记录**：加 `--group-leftovers-as other`。
- **样本测序量差异大、想横向可比**：用 `-n columns_before_collapsing` 之类先按列标准化。
- **不想动默认路径**：FAPROTAX 装在别处时再用 `-g` / `--collapse-table-path` 覆盖。

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input-table` | 必填 |  | 输入OTU/ASV表（BIOM或TSV格式）｜Input OTU/ASV table (BIOM or TSV format) |
| `-o, --output-dir` | `./faprotaxtax_output` |  | 输出目录｜Output directory |
| `-g, --groups-file` | — |  | FAPROTAX功能组数据库文件路径｜Path to FAPROTAX groups database file |
| `--collapse-table-path` | — |  | collapse_table.py脚本路径｜Path to collapse_table.py script |
| `--collapse-by-metadata` | — |  | 用于功能注释的BIOM元数据字段名｜BIOM metadata field for functional annotation |
| `--group-leftovers-as` | — |  | 未匹配到功能组的记录归为此组名｜Group name for unmatched records |
| `-n, --normalize` | `none` | none/columns_before_collapsing/rows_before_collapsing/columns_after_collapsing/rows_after_collapsing/columns_before_collapsing_excluding_unassigned/rows_before_collapsing_excluding_unassigned | 标准化方式｜Normalization method |
| `--average` | `none` | none/across_records/across_group_members/across_used_group_members/maximum/minimum/minimum_across_records | 组内聚合方式｜Aggregation method within groups |
| `--row-names-are-in-column` | — |  | 包含行名的列名或索引｜Column name or index containing row names |
| `--output-format` | `auto` | auto/BIOM/classical | 输出格式｜Output format |
| `-t, --threads` | `1` | int | 线程数｜Number of threads |
| `-f, --force` | — |  | 覆盖已存在的输出文件｜Overwrite existing output files |
| `-v, --verbose` | — |  | 详细输出｜Verbose output |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input-table` | 必填 |  | 输入OTU/ASV表（BIOM或TSV格式）｜Input OTU/ASV table (BIOM or TSV format) |
| `-o, --output-dir` | `./faprotaxtax_output` |  | 输出目录｜Output directory (default: ./faprotaxtax_output) |
| `-g, --groups-file` | — |  | FAPROTAX功能组数据库文件路径｜Path to FAPROTAX functional groups database file |
| `--collapse-table-path` | — |  | collapse_table.py脚本路径｜Path to collapse_table.py script |
| `--collapse-by-metadata` | — |  | 用于功能注释的BIOM元数据字段名（如: taxonomy）｜BIOM metadata field for functional annotation (e.g., taxonomy) |
| `--group-leftovers-as` | — |  | 未匹配到功能组的记录归为此组名｜Group name for records not matching any functional group |
| `-n, --normalize` | `none` | none/columns_before_collapsing/rows_before_collapsing/columns_after_collapsing/rows_after_collapsing/columns_before_collapsing_excluding_unassigned/rows_before_collapsing_excluding_unassigned | 标准化方式｜Normalization method (default: none) |
| `--average` | `none` | none/across_records/across_group_members/across_used_group_members/maximum/minimum/minimum_across_records | 组内聚合方式｜Aggregation method within groups (default: none) |
| `--row-names-are-in-column` | — |  | 包含行名的列名或索引（仅经典TSV表格）｜Column name or index containing row names (for TSV tables only) |
| `--output-format` | `auto` | auto/BIOM/classical | 输出格式｜Output format (default: auto) |
| `-t, --threads` | `1` | int | 线程数｜Number of threads (default: 1) |
| `-f, --force` | — | store_true | 覆盖已存在的输出文件｜Overwrite existing output files |
| `-v, --verbose` | — | store_true | 详细输出｜Verbose output |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies { #dependencies }

- FAPROTAX 软件包（collapse_table.py + FAPROTAX.txt，默认 `~/software/FAPROTAX/FAPROTAX_1.2.12/`）
- Python 3（用当前解释器 `sys.executable` 直接运行 collapse_table.py，**不走 conda run**）

## 常见问题 | FAQ { #faq }

**Q1：支持断点续传吗？**
支持。检测到 `01_collapsed/functional_table.tsv`（或 .biom）已存在即跳过整个流程并只补写版本信息。换标准化/聚合参数重跑前，先删掉 `01_collapsed/` 里的旧产物，否则会复用旧结果。

**Q2：报「功能组数据库不存在」或「collapse_table.py 脚本不存在」？**
这是 config.validate 的检查。默认路径是 `~/software/FAPROTAX/FAPROTAX_1.2.12/`，若你的 FAPROTAX 装在别处，用 `-g` 和 `--collapse-table-path` 显式指定。

**Q3：TSV 输入为什么报错/结果为空？**
经典 TSV 表必须让程序知道哪一列是分类学信息。用 `--collapse-by-metadata`（列名）或 `--row-names-are-in-column`（行名所在列）指定。

**Q4：线程数设了 12 会更快吗？**
不会明显变快。collapse_table.py 是单线程 Python 脚本，`-t` 默认 1 即可。

**Q5：输出是 .tsv 还是 .biom？**
默认 auto：BIOM 输入且输出格式非 classical 时输出 .biom，否则 .tsv；可用 `--output-format` 强制指定 BIOM 或 classical。