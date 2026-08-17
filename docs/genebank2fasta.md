# GenBank 序列提取 | GenBank Sequence Extraction

一句话理解：**把一堆 GenBank 格式文件里的基因批量「拆解」出来**——抽每个基因的 CDS（编码序列）和蛋白序列，按样本/按基因分别整理成 FASTA，还能顺手统计核心基因、生成系统发育矩阵。

## 功能概述 | Overview { #overview }

- 从 GenBank 文件（.gb/.gbk/.genbank）提取所有 CDS 特征并翻译成蛋白序列
- 按样本（by_sample）和按基因（by_gene）两种方式分别输出，默认都开
- 默认跳过标记为 unknown 的基因，可按蛋白长度过滤
- 自动识别核心基因（在 ≥90% 样本中出现的基因），生成列表
- 可选 `--phylo` 生成系统发育超矩阵（supermatrix + 分区文件）
- 纯 Python（Biopython）实现，多线程并行处理

## 快速开始 | Quick Start { #quick-start }

```bash
biopytools genebank2fasta -i /path/to/genbank -o ./output
```

最小输入：一个装着若干 GenBank 文件的目录。

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语<br>Term | 通俗理解 |
|---|---|
| GenBank 文件 | NCBI 的标准注释格式，一个文件里既有序列又有基因注释，像「序列+说明书」合订本 |
| CDS | 编码序列，真正翻译成蛋白的那段 DNA |
| 翻译表 | 密码子→氨基酸的对照规则，细菌常用第 11 号表，工具默认按它翻译 |
| by_sample / by_gene | 两种整理方式：按「谁家样本」分，或按「什么基因」跨样本合并 |
| 核心基因 | 在绝大多数样本里都存在的基因，像物种的「必备零件清单」 |
| 超矩阵 | 把多个基因的序列头尾拼成一条，用于建系统发育树 |

## 输入 | Input { #input }

一个目录，目录下放 GenBank 文件（扩展名 `.gb`/`.gbk`/`.genbank`）。工具自动扫描目录找文件，样本名取文件名（去掉扩展名）。

## 参数说明 | Parameters { #parameters }

### 输出组织 | Output organization

**通俗理解|In plain words:** 默认按样本和按基因各出一套，方便不同用途；`--no-sample-sep` 关掉按样本，`--no-gene-sep` 关掉按基因。`--keep-unknown` 打开后，那些名字就叫 unknown 的基因也会被保留（默认丢弃，因为它们通常没命名价值）。

相关参数：`--no-sample-sep`、`--no-gene-sep`、`--keep-unknown`。

### 过滤与可选分析 | Filtering & optional analysis

**通俗理解|In plain words:** `--min-length` 是蛋白最短长度（氨基酸数），设大一点能滤掉碎片化的短 ORF；`--phylo` 打开后额外生成系统发育矩阵，供后续建树。`-t` 线程数影响并行解析速度，默认 12 够用。

相关参数：`--min-length`（默认 10）、`--phylo`、`-t/--threads`（默认 12）。

## 分析流程 | Pipeline { #pipeline }

```text
GenBank 目录
    │
    ▼
扫描并解析每个文件(Biopython, 并行)
    │
    ▼
提取 CDS + 翻译蛋白 + 长度过滤
    │
    ▼
按样本 / 按基因 写 FASTA
    │
    ▼
统计报告 + 核心基因列表 (+ 可选系统发育矩阵)
```

## 输出 | Output { #output }

```text
output/
├── cds/
│   ├── by_sample/{sample}_cds.fasta     # 每样本一个 CDS 文件
│   └── by_gene/{gene}.fasta             # 每基因一个 CDS 文件(跨样本合并)
├── pep/
│   ├── by_sample/{sample}_pep.fasta     # 每样本一个蛋白文件
│   └── by_gene/{gene}.fasta             # 每基因一个蛋白文件
├── extraction_summary.txt               # 提取统计报告
├── core_genes_list.txt                  # 核心基因列表(≥90% 样本)
└── phylogenetic_matrix/                 # 仅 --phylo 时
    ├── supermatrix.fasta                # 拼接超矩阵
    └── partitions.txt                   # 分区文件(各基因在矩阵中的区间)
```

FASTA 头格式：按样本为 `>{gene} {organism}|{product}|{长度}{bp/aa}`；按基因为 `>{sample} {organism} {gene}|{product}|{长度}{bp/aa}`。

## 结果解读 | Interpreting Results { #interpreting-results }

- `by_sample/`：看单个样本提取了多少基因、序列是否完整
- `by_gene/`：把同名基因跨样本放在一个文件里，直接用于多序列比对
- `extraction_summary.txt`：含总体统计、基因统计（按样本数排序）、样本统计（Top 20）、核心基因清单
- `core_genes_list.txt`：在 ≥90% 样本中出现的基因，可用于 `--phylo` 或单拷贝基因建树
- CDS 长度不是 3 的倍数、蛋白过短等都会记 warning，但不阻断，看日志可排查

## 参数选择建议 | Parameter Guidance { #parameter-guidance }

- 做比较基因组/系统发育：保持默认（双输出）+ `--phylo`，用核心基因建超矩阵
- 只要每个样本的蛋白序列：`--no-gene-sep`，省掉按基因的冗余输出
- 输入含大量碎片 ORF：调大 `--min-length`（如 30、50）过滤噪音
- 需要保留未命名基因：加 `--keep-unknown`

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input, -i` | 必填 |  | GenBank文件目录｜Input GenBank files directory |
| `--output-dir, -o` | `./genbank_output` | Path | 输出目录｜Output directory |
| `--threads, -t` | `12` | int | 并行线程数｜Number of parallel threads |
| `--min-length` | `10` | int | 最小蛋白长度(氨基酸)｜Minimum protein length (amino acids) |
| `--phylo, --create-phylogenetic-matrix` | — |  | 创建系统发育分析矩阵｜Create phylogenetic analysis matrix |
| `--no-sample-sep` | — |  | 不按样本分离输出｜Do not separate output by sample |
| `--no-gene-sep` | — |  | 不按基因分离输出｜Do not separate output by gene |
| `--keep-unknown` | — |  | 保留未知基因｜Keep unknown genes |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 输入GenBank文件目录｜Input GenBank files directory |
| `-o, --output-dir` | `./genbank_output` |  | 输出目录｜Output directory |
| `-t, --threads` | `88` | int | 并行线程数｜Number of parallel threads |
| `--min-length` | `10` | int | 最小蛋白质长度(氨基酸)｜Minimum protein length (amino acids) |
| `--phylo, --create-phylogenetic-matrix` | — | store_true | 创建系统发育分析矩阵｜Create phylogenetic analysis matrix |
| `--no-sample-sep` | — | store_true | 不按样品分离输出｜Do not separate output by sample |
| `--no-gene-sep` | — | store_true | 不按基因分离输出｜Do not separate output by gene |
| `--keep-unknown` | — | store_true | 保留unknown基因｜Keep unknown genes |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies { #dependencies }

- Python 3 + Biopython（`Bio.SeqIO`）
- 无 conda 环境、无其他外部生信软件依赖

## 常见问题 | FAQ { #faq }

**Q1：会断点续传吗？**
不会。每次运行都重新解析并覆盖输出，多线程并行处理速度快，无需续传。

**Q2：为什么有些基因没被提取？**
三处会被丢弃：基因名是 unknown（除非 `--keep-unknown`）、蛋白长度小于 `--min-length`、CDS 提取异常。逐一对应检查日志 warning。

**Q3：翻译表用的是什么？**
优先读 GenBank 里每个 CDS 的 `transl_table` 限定符，缺省回退到第 11 号表（细菌标准）。若你的物种用别的表且文件里没标，翻译结果可能不对。

**Q4：`--phylo` 依赖什么？**
它读取 `core_genes_list.txt`（核心基因）来拼超矩阵，该列表在同一轮运行的统计阶段先生成，正常一次带 `--phylo` 即可；核心基因为空时不会产出矩阵。
