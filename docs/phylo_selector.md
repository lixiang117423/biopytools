# 系统发育树样品选择 | Phylogenetic Tree Sample Selection

一句话理解：**从成百上千个样品里，自动挑出一批「代表性」样品，既覆盖得均匀、又不重复，方便你缩小后续分析规模**。

输入一份样品层级表（Excel）+ 一份 PCA 坐标，输出按「均匀间隔 + PCA 去重」挑出的样品名单。

## 功能概述 | Overview { #overview }

- 从样品全集里按「均匀间隔」抽出一批代表性样品（默认 150 个）
- 用 PCA 坐标做去重：在 PCA 空间里太接近（近乎重复）的样品会被替换
- 自动避免选中「顺序上相邻」的样品，保证覆盖更均匀
- 输出四种文件：样品名单（txt/csv）、详细报告、可视化索引
- 断点续传：无（每次重新选择）

## 快速开始 | Quick Start { #quick-start }

`@bash
biopytools phylo-selector -x hierarchy.xlsx -p pca.txt -o selected_samples -n 150
`@

最小输入：一个层级 Excel 文件（含 `label` 列）+ 一个 PCA 坐标文件。

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语 | 通俗理解 |
|------|----------|
| 样品 | 你要分析的一条序列/一个个体 |
| PCA 坐标 | 把高维差异压成几个「主成分」后的坐标，坐标近=整体相似 |
| PCA 距离 | 两个样品在主成分空间里的欧氏距离，越小越像 |
| 去重阈值 | 「多近算重复」的界限，距离小于它就认为两个样品是重复的 |
| 均匀间隔 | 从排好序的样品里每隔固定个数取一个，保证覆盖全貌 |

## 输入 | Input { #input }

### 层级文件（Excel，必需）

一个 `.xlsx` 文件，必须含一列名为 `label`（样品名），其余列是 `parent_1`、`parent_2` ……（层级父节点）。程序实际只用 `label` 列来确定样品名单和顺序。

### PCA 坐标文件（必需）

制表符或空格分隔，常见 PLINK 的 `.eigenvec` 格式：第1列 FID、第2列 IID（样品名，须与层级文件的 `label` 一致）、第3列起是 PC1、PC2 ……（程序取前 10 个主成分）。可带表头行（以 FID 或 # 开头会自动跳过）。

`@text
FID    IID       PC1       PC2
S1     S1        0.0123    -0.0045
S2     S2        -0.0087   0.0210
`@

### 分组文件（可选）

两列（样品名、分组名），当前版本会读入但选择算法暂未使用它。

## 参数说明 | Parameters { #parameters }

### 必需参数 | Required

**通俗理解|In plain words:** 层级表、PCA 表、输出前缀，三个必给。输出前缀决定所有结果文件名。

### 选择数量与分组 | Selection size & grouping

**通俗理解|In plain words:** `-n/--n-samples` 管「总共选几个」（默认 150）。**注意**：如果层级表里能在 PCA 里找到的有效样品不足 `-n`，会退而求其次把能选的都选上。`--hierarchy-level` 与 `--min-samples-per-group` 在当前实现里基本不参与计算（保留参数），**一般不用动**；`-g/--group-file` 目前读入后也暂未用于选择，可先不提供。

### PCA 去重 | PCA dedup

**通俗理解|In plain words:** `--pca-threshold` 管「多近算重复」——两个样品 PCA 距离小于它就认为重复、触发替换。默认 0.0001 已经很严（几乎只在坐标完全重合时才去重）；调大（如 0.001、0.01）会去掉更多「看着相似」的样品，保留更分散的代表。**一般不用动**。

### 输出控制 | Output control

**通俗理解|In plain words:** 三个开关分别关掉详细报告、CSV、可视化文件。默认全开，一般不用动。

## 分析流程 | Pipeline { #pipeline }

`@text
步骤1: 读层级表（取 label 列） -> 读 PCA -> 读分组(可选)
   |
   v
步骤2: 取「层级表 label 列里且 PCA 里也有坐标」的样品，作为候选全集
   |
   v
步骤3: 均匀间隔选择 -> 避免选中顺序相邻的样品
   |
   v
步骤4: PCA 去重（距离 < 阈值的对，保留靠前、替换靠后）
   |
   v
步骤5: 写结果（txt/csv/报告/可视化）
`@

## 输出 | Output { #output }

`@text
<前缀>.txt                 # 选中的样品名单（每行一个样品名，带 # 注释头）
<前缀>_report.txt          # 详细报告（策略、参数、样品列表，带序号）
<前缀>.csv                 # CSV 名单（Index,Sample 两列）
<前缀>_visualization.txt   # 可视化索引（每行一个序号，供画图标注）
`@

## 结果解读 | Interpreting Results { #interpreting }

- **`<前缀>.txt`**：最终要用的样品名单，直接把它当「取哪些样品」的清单喂给下游（如 `biopytools seq-extract` 抽取序列）。
- **`<前缀>.csv`**：和 txt 同内容，便于在 Excel/脚本里读。
- **`<前缀>_report.txt`**：人看的报告，含选择策略、参数、带序号的完整名单，核对用。
- **好坏判据**：日志末尾会打印「最小 PCA 距离 / 平均 PCA 距离」。若仍提示有样品对距离小于阈值，说明候选集里重复度较高；最终名单里的样品在 PCA 图上应大致分散、不扎堆。

## 参数选择建议 | Parameter Guidance { #guidance }

- **常规抽代表样品**：只给 `-x -p -o -n`，其余默认。
- **样品很多（成千上万）**：直接设 `-n` 为你要保留的数量。
- **发现选出的样品在 PCA 图上还扎堆**：把 `--pca-threshold` 调大一点（如 0.001、0.01）重新跑。
- **不需要 CSV/报告**：`--no-csv`、`--no-report` 关掉。

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-x, --hierarchy-file` | 必填 |  | 层级关系文件（必需）｜Hierarchy file (required) |
| `-p, --pca-file` | 必填 |  | PCA坐标文件（必需）｜PCA coordinates file (required) |
| `-o, --output-prefix` | 必填 |  | 输出文件前缀｜Output file prefix |
| `-n, --n-samples` | `150` | int | 选择样品总数｜Total number of samples to select |
| `-g, --group-file` | — |  | 样品分组表文件｜Sample group table file |
| `--hierarchy-level` | `10` | int | 层级深度（用于分组）｜Hierarchy level for grouping |
| `--pca-threshold` | `0.0001` | float | PCA去重阈值｜PCA dedup threshold |
| `--min-samples-per-group` | `1` | int | 每组最小样品数｜Minimum samples per group |
| `--no-report` | — |  | 不生成详细报告｜Do not generate detailed report |
| `--no-csv` | — |  | 不生成CSV文件｜Do not generate CSV file |
| `--no-visualization` | — |  | 不生成可视化｜Do not generate visualization |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-x, --hierarchy-file` | 必填 |  | 层级关系文件（必需）｜Hierarchy file (required) |
| `-p, --pca-file` | 必填 |  | PCA坐标文件（必需）｜PCA coordinates file (required) |
| `-o, --output-prefix` | 必填 |  | 输出文件前缀｜Output file prefix |
| `-n, --n-samples` | `150` | int | 选择样品总数｜Total number of samples to select (default: 150) |
| `-g, --group-file` | — |  | 样品分组表文件｜Sample group table file |
| `--hierarchy-level` | `10` | int | 层级深度（用于分组）｜Hierarchy level for grouping (default: 10) |
| `--pca-threshold` | `0.0001` | float | PCA去重阈值｜PCA dedup threshold (default: 0.0001) |
| `--min-samples-per-group` | `1` | int | 每组最小样品数｜Minimum samples per group (default: 1) |
| `--no-report` | — | store_true | 不生成详细报告｜Do not generate detailed report |
| `--no-csv` | — | store_true | 不生成CSV文件｜Do not generate CSV file |
| `--no-visualization` | — | store_true | 不生成可视化｜Do not generate visualization |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies { #dependencies }

- Python 3 + pandas、openpyxl（读 Excel）、numpy、scipy（PCA 距离计算）
- 无外部命令行软件依赖

## 常见问题 | FAQ { #faq }

**Q1：报错说「hierarchy file not found」或读不了 Excel？**
层级文件必须是 `.xlsx`（Excel 格式），并且要有一列名字精确叫 `label`。旧版 `.xls` 或改了列名都会失败。

**Q2：为什么选出来的样品数比 `-n` 少？**
因为「层级表 label 里、且 PCA 里也有坐标」的有效样品数不足 `-n`，程序会把所有有效样品都选上（日志会提示）。请核对 PCA 文件的 IID（第2列）是否和层级表 `label` 一致。

**Q3：`--hierarchy-level` 和 `--min-samples-per-group` 有用吗？**
当前版本里它们是保留参数，基本不参与计算，设了也不影响结果，放心忽略。

**Q4：PCA 文件要不要表头？**
可有可无。有表头（以 FID 或 # 开头）会自动跳过第一行；无表头直接从第一行读。关键是要有 FID、IID 两列 + 后面的 PC 值。

**Q5：选出来的样品能保证一定不重复吗？**
算法会用 PCA 去重 + 避免相邻，尽力而为，但不做硬性保证。日志末尾若提示仍有距离过近的对，可调大 `--pca-threshold` 重跑。
