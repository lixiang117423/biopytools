# VCF 变异基因注释 | VCF Variant Gene Annotation

一句话理解：**把 VCF 里的每个变异位点，按 GFF 注释标注它落在哪个基因、哪个部位（外显子/内含子/UTR/基因间区）**，快速把「一堆变异」变成「哪些变异影响了哪些基因」。

## 功能概述 | Overview { #overview }

- 基于 GFF 注释，为每个 VCF 变异位点标注所属基因与基因部位
- 部位按优先级取最底层：外显子 > CDS > UTR > 内含子 > 转录本 > 基因
- 支持 .vcf/.vcf.gz 和 .gff/.gff.gz，纯 Python 流式处理，内存占用低
- 多等位位点（ALT 含逗号）自动拆成多行
- 未落在任何基因的位点标注为 intergenic

## 快速开始 | Quick Start { #quick-start }

```bash
biopytools vcf2gene -i variants.vcf -g annotation.gff -o annotated.txt
```

最小输入：VCF + GFF 注释，指定输出文件。

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语<br>Term | 通俗理解 |
|---|---|
| VCF | 记录变异位点的标准格式，每个位点一行 |
| 外显子 | 编码蛋白的「正文」片段 |
| 内含子 | 外显子之间会被切掉的「间隔」片段 |
| UTR | 非翻译区，位于编码区前后、不翻译成蛋白 |
| 基因间区 | 两个基因之间的「空白地带」 |
| 注释优先级 | 一个位点可能同时落在外显子和基因上，取更细粒度的那个（外显子优先） |

## 输入 | Input { #input }

### VCF 文件

标准 VCF（`.vcf` 或 `.vcf.gz`），读取 CHROM/POS/REF/ALT 四列做注释。

### GFF 注释文件

GFF3 格式（`.gff` 或 `.gff.gz`），属性须用 `=` 分隔（`ID=...;Parent=...`）。需要基因有 ID，且子特征通过 Parent 关联到基因，才能正确追溯 gene_id：

```text
Chr01	source	gene	1000	5000	.	+	.	ID=Gene001
Chr01	source	mRNA	1000	5000	.	+	.	ID=Gene001.1;Parent=Gene001
```

## 参数说明 | Parameters { #parameters }

### 必需参数 | Required

**通俗理解|In plain words:** 三个必填：VCF 提供位点、GFF 提供注释、输出文件写结果。

相关参数：`-i/--vcf`、`-g/--gff`、`-o/--output`。

### 线程 | Threads

**通俗理解|In plain words:** `-t/--threads` 当前为预留参数，尚未实际参与计算，用默认值即可。

相关参数：`-t/--threads`（默认 12）。

## 分析流程 | Pipeline { #pipeline }

```text
GFF 注释文件
    |
    v
解析基因/特征区间,建立索引
    |
    v
流式逐行读 VCF
    |
    v
二分查找位点所在特征(按优先级取最底层)
    |
    v
写输出: Chr Pos Ref Alt Gene Type
```

## 输出 | Output { #output }

单个制表符分隔文本文件，表头 `Chr\tPos\tRef\tAlt\tGene\tType`：

| 列名 | 含义 |
|---|---|
| Chr / Pos | 变异所在染色体与位置 |
| Ref / Alt | 参考与替代等位基因（多 ALT 拆成多行） |
| Gene | 所属基因 ID，基因间区为 `intergenic` |
| Type | 命中部位类型（exon/cds/five_prime_utr/three_prime_utr/utr/intron/mrna/transcript/gene），基因间区为 `intergenic` |

## 结果解读 | Interpreting Results { #interpreting-results }

- `Type=exon`/`cds`：变异落在编码区，最可能影响蛋白
- `Type=intron`：落在内含子，可能影响剪接
- `Type=*_utr`/`utr`：落在非翻译区
- `Gene=intergenic`：不落在任何基因，属基因间区变异
- 若结果里大量 `intergenic` 且本应有基因命中，多半是 GFF 格式问题（见 FAQ）

## 参数选择建议 | Parameter Guidance { #parameter-guidance }

- 默认参数即可满足绝大多数需求
- 只关心编码区变异：对输出按 `Type` 列过滤 exon/cds
- 输入巨大时无需特殊设置，工具本身流式处理

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--vcf, -i` | 必填 |  | 输入VCF文件｜Input VCF file path |
| `--gff, -g` | 必填 |  | 输入GFF注释文件｜Input GFF annotation file path |
| `--output, -o` | 必填 |  | 输出结果文件｜Output result file path |
| `--threads, -t` | `12` | int | 线程数｜Number of threads (reserved for future use) |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --vcf` | 必填 |  | 输入VCF文件路径｜Input VCF file path |
| `-g, --gff` | 必填 |  | 输入GFF注释文件路径｜Input GFF annotation file path |
| `-o, --output` | 必填 |  | 输出结果文件路径｜Output result file path |
| `-t, --threads` | `12` | int | 线程数｜Number of threads (reserved for future use) |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies { #dependencies }

- Python 3 标准库（gzip、bisect、collections）
- 无外部生信软件、无 conda 环境依赖

## 常见问题 | FAQ { #faq }

**Q1：会断点续传吗？**
不会。每次运行重新解析并覆盖输出文件。

**Q2：为什么很多位点标成了 intergenic？**
最可能是 GFF 用了 GTF 风格（`gene_id "xxx"`，空格+引号），本工具只认 GFF3 的 `=` 风格，导致基因/特征关系解析失败。换成 GFF3 格式即可。

**Q3：多等位位点怎么处理？**
ALT 列含多个等位基因时，每个 ALT 拆成独立一行输出，Gene/Type 相同。

**Q4：`-t` 线程数有用吗？**
当前未实际使用，是预留参数，保持默认即可。
