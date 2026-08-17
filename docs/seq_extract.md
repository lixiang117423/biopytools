# seq_extract - 序列提取(seqkit 封装) | Sequence Extraction (seqkit wrapper)

一句话理解：**从一个大 FASTA 里按「名字或区间」把想要的序列捞出来——可以给一个 ID、一个 ID 清单、或一个 BED 区间文件，它自动判断你给的是哪种，用 seqkit 提取成新的 FASTA。**

## 功能概述 | Overview { #overview }

- seqkit 封装，自动识别三种查询：单个 ID、ID 文件（一列）、BED 文件（至少两列）
- 单 ID 用 seqkit grep -p，ID 文件用 seqkit grep -f，BED 用 seqkit subseq --bed
- 输出文件名可省略，自动推导为 {查询}.{目标}.fa
- 可用 --bed 强制按 BED 模式处理

## 快速开始 | Quick Start { #quick-start }

```bash
biopytools seq-extract -i gene.id.txt -s gene.fa -o gene.genomic.fa
```

最小输入：一个查询（ID / ID 文件 / BED 文件）+ 一个目标序列 FASTA。

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语 | 通俗理解 |
|------|----------|
| FASTA | 纯文本存序列的格式，> 开头是名字，下面是序列 |
| 序列 ID | FASTA 里 > 后面的名字，用来唯一指代一条序列 |
| ID 文件 | 一列 ID 的清单，一行一个，批量「点名」要哪些序列 |
| BED | 用「染色体 起始 终止」描述一段区间，按位置捞序列 |
| seqkit | 常用的序列处理命令行工具，这里负责实际提取 |

## 输入 | Input { #input }

- 查询(-i)：单个 ID（非文件路径字符串）、ID 文件（一列）、或 BED 文件（至少两列）。
- 目标序列 FASTA(-s)：要从中提取的序列库。

自动识别规则：-i 不是文件路径则当单个 ID；是文件则读首行，按制表符分列若列数 >= 2 判为 BED，否则判为 ID 文件。

示例（ID 文件）：

```text
geneA
geneB
geneC
```

示例（BED 文件，两列即可）：

```text
chr1	100	250
chr2	10	60
```

## 参数说明 | Parameters { #parameters }

### 必需与输出 | Required & output

**通俗理解|In plain words:** -i 是查询，-s 是目标 FASTA，两者必填；-o 输出文件可省略，省略时自动命名为 {查询文件名去扩展名}.{目标文件名去扩展名}.fa。

### 模式 | Mode

**通俗理解|In plain words:** --bed 强制按 BED 模式处理（跳过自动检测）。只有当自动检测判错、或你的 BED 文件首行恰好只有一列时才需要用，一般不用加。

## 分析流程 | Pipeline { #pipeline }

```text
自动检测查询类型（single_id / id_file / bed_file）
    |
    v
推导输出文件名（若未指定 -o）
    |
    v
按类型调用 seqkit：
  single_id -> seqkit grep -p <id>
  id_file   -> seqkit grep -f <id文件>
  bed_file  -> seqkit subseq --bed <bed文件>
    |
    v
写输出 FASTA
```

## 输出 | Output { #output }

```text
gene.genomic.fa      # 提取出的序列（文件名由 -o 指定或自动推导）
```

单文件 FASTA 输出，序列名与目标 FASTA 中的原始名字一致。

## 结果解读 | Interpreting Results { #interpreting-results }

**通俗理解|In plain words:** 输出就是你要的序列集合，条数应与查询命中数一致。

- ID 模式：输出里是 ID 命中（精确按完整 ID 列表匹配）的那些序列
- 单 ID 模式：用正则/子串匹配，可能把「名字里包含该字符串」的其它序列也带出来（见 FAQ）
- BED 模式：按区间坐标提取，返回每段区间对应的序列

好坏判据：输出条数大致等于期望条数即成功；明显偏多（单 ID 模式）或偏少（有 ID/BED 没命中）都要检查查询内容。

## 参数选择建议 | Parameter Guidance { #parameter-guidance }

- 按名单批量提取：用 ID 文件（每行一个完整 ID），最可靠
- 只取一个已知 ID：直接 -i geneA
- 按位置取：用 BED 文件（自动识别，必要时加 --bed）
- 不确定自动检测结果：先跑一次看日志里的「查询类型」提示

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 查询:单个ID、ID文件(一列)或BED文件(>=2列)｜Query: single ID, ID file (1 column), or BED file (>=2 columns) |
| `-s, --sequence` | 必填 |  | 目标序列FASTA文件｜Target sequence FASTA file |
| `-o, --output` | — |  | 输出文件(默认自动推导:{query}.{subject}.fa)｜Output file (default: auto-derived) |
| `--bed` | — |  | 强制BED模式(跳过自动检测)｜Force BED mode (skip auto-detection) |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 查询:单个ID、ID文件(一列)或BED文件(>=2列)｜Query: single ID, ID file (1 column), or BED file (>=2 columns) |
| `-s, --sequence` | 必填 |  | 目标序列FASTA文件｜Target sequence FASTA file |
| `-o, --output` | — |  | 输出文件(默认自动推导:{query}.{subject}.fa)｜Output file (default: auto-derived {query}.{subject}.fa) |
| `--bed` | — | store_true | 强制BED模式(跳过自动检测)｜Force BED mode (skip auto-detection) |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies { #dependencies }

- seqkit（默认命令 seqkit，可通过 SEQKIT_PATH 环境变量或配置文件指定路径，自动检测 conda 环境用 conda run 调用）
- Python 3

## 常见问题 | FAQ { #faq }

Q1：支持断点续传吗？
不支持。每次运行都重新提取并覆盖输出文件。

Q2：单 ID 模式为什么多提了一些序列？
单 ID 模式用 seqkit grep -p（正则/子串匹配），会把 ID 里「包含该字符串」的序列也带出来。想要精确匹配，改用 ID 文件（grep -f 按完整 ID 列表匹配）。

Q3：怎么确认它判断成了哪种模式？
运行日志里会打印「查询类型|Query type: single_id / id_file / bed_file」，据此判断；判断错了用 --bed 强制指定或调整输入文件内容。

Q4：输出文件名是什么？
指定了 -o 就用它；否则自动为 {查询文件名去扩展名}.{目标文件名去扩展名}.fa。当 -i 是单个 ID（非文件）时，查询部分就是该 ID 字符串本身。

Q5：seqkit 找不到怎么办？
默认调用 seqkit（可在 PATH 或 conda 环境里）；也可用 SEQKIT_PATH 环境变量或 ~/.config/biopytools/config.yml 指定路径。