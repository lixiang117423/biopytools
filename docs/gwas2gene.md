# GWAS 候选基因筛选 | GWAS Candidate Gene Finder (gwas2gene)

一句话理解：**把 GWAS 找出的显著位点，就近对应到基因组上的基因**。输入 GWAS 结果和基因注释(GFF3)，自动挑出显著位点并抓取其上下游窗口内的候选基因，可选附带基因功能注释。

## 功能概述 | Overview

- 按 P 值阈值筛选显著 SNP，自动识别染色体/位置/P 值列（列名或列号均可）
- 在指定上下游窗口内从 GFF3 注释提取候选基因（含转录本 ID）
- 自动统一 GWAS 与 GFF 的染色体命名（Chr1 vs 1），无需手动对齐
- 计算 SNP 到基因的链方向感知距离（0=基因内，负=上游，正=下游）
- 可选整合功能注释文件（基因 ID → 功能描述）

## 快速开始 | Quick Start

```bash
biopytools gwas2gene -g gwas.txt -p P -t 1e-5 --gff annotation.gff3 -o candidates.tsv
```

最小输入：一个 GWAS 结果文件 + 一个 GFF3 注释文件 + 输出文件路径；-p 指定 P 值列（列名或 1-based 列号）。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| GWAS 结果 | 一张「每个位点与性状相关程度」的打分表，关键是 P 值那一列 |
| GFF3 | 基因组注释文件，记录「哪些位置有哪些基因、起止坐标、方向」 |
| 上下游窗口 | 以显著位点为中心，左右各延伸一段距离（如 200kb）划定的搜索范围 |
| 距离(distance) | SNP 到基因的距离：0=落在基因里，负=基因上游，正=基因下游 |

## 输入 | Input

### GWAS 结果文件

制表符分隔且**带表头**，需含染色体列、位置列、P 值列。染色体/位置列名含 chr/chrom、pos/position/bp 等关键词即自动识别，P 值列用 -p 指定：

```text
CHR    POS         SNP          P
1      14370       rs6054257    2.1e-08
1      17330       .            1.5e-06
```

### GFF3 注释文件

标准 GFF3，解析 gene 与 mRNA/transcript 行。

### 功能注释文件（可选）

两列：基因 ID（或转录本 ID）+ 功能描述。

## 输出 | Output

输出单个 TSV 文件，每行一个「显著 SNP × 窗口内基因」的组合：

```text
snp_id  snp_chrom  snp_pos  snp_pval  gene_id  gene_name  transcript_id  gene_chrom  gene_start  gene_end  gene_strand  distance  [function]
```

| 列 | 含义 |
|------|------|
| snp_id / snp_chrom / snp_pos / snp_pval | 显著 SNP 的名称、染色体、位置、P 值 |
| gene_id / gene_name | 基因 ID 与名称 |
| transcript_id | 该基因的转录本 ID（多个用逗号分隔，无则 NA） |
| gene_chrom / gene_start / gene_end / gene_strand | 基因坐标与链方向 |
| distance | SNP 到基因的距离（0=基因内，负=上游，正=下游） |
| function | 功能描述（仅当提供 -f 功能文件） |

## 结果解读 | Interpreting Results

**优先看 distance 为 0 的记录**（SNP 直接落在基因内部），其次看 |distance| 较小、且 function 有意义的记录。同一 SNP 可能对应多个基因，距离越近、功能越相关的基因越值得优先验证。

## 参数选择建议 | Parameter Guidance

- -t / --threshold：默认 1e-5；显著位点太多就调小（如 1e-6），太少就调大
- -w / --window：默认 200000（200kb）；基因稀疏的物种可加大到 500kb，基因密集则减小
- -p / --pval-col：给列名（大小写不敏感）或 1-based 列号均可，一般给列名最稳妥
- -f / --func：有基因功能注释表时加上，输出多一列 function 方便判断

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--gwas, -g` | 必填 | Path | GWAS结果文件路径｜GWAS result file path |
| `--pval-col, -p` | 必填 | str | P值所在列名或列号（1-based）｜P-value column name or index (1-based) |
| `--threshold, -t` | `1e-05` | float | P值阈值｜P-value threshold |
| `--window, -w` | `200000` | int | 上下游窗口大小｜Window size upstream/downstream in bp |
| `--gff` | 必填 | Path | GFF3注释文件路径｜GFF3 annotation file path |
| `--output, -o` | 必填 | str | 输出文件路径｜Output file path |
| `--func, -f` | — | Path | 功能注释文件路径｜Function annotation file path (optional) |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--gwas` | 必填 |  | GWAS结果文件路径｜GWAS result file path |
| `--pval-col` | 必填 |  | P值所在列名或列号（1-based）｜P-value column name or index (1-based) |
| `--threshold` | `1e-05` | float | P值阈值｜P-value threshold (default: 1e-5) |
| `--window` | `200000` | int | 上下游窗口大小｜Window size upstream/downstream in bp (default: 200000) |
| `--gff` | 必填 |  | GFF3注释文件路径｜GFF3 annotation file path |
| `--output` | 必填 |  | 输出文件路径｜Output file path |
| `--func` | — |  | 功能注释文件路径｜Function annotation file path (optional) |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- Python 3（标准库即可，无第三方依赖）

## 常见问题 | FAQ

**Q1：报「无法自动识别染色体列或位置列」？**
程序按列名关键词猜染色体(chr/chrom)和位置(pos/position/bp)列。请确保 GWAS 结果有表头、且这些列名包含相应关键词；或把列名改成 CHR、POS 这样的标准写法。

**Q2：-p 给列名还是列号？**
两者都支持：列号是 1-based（第 1 列=1）；列名大小写不敏感。给列名更不易出错。

**Q3：GWAS 和 GFF 的染色体一个写 1、一个写 Chr1，能对上吗？**
能。程序会自动检测双方格式并统一（Chr1↔1），无需手动改文件。

**Q4：输出里 distance 是负数什么意思？**
distance 按链方向计算：0=SNP 落在基因内，负值=基因上游，正值=基因下游。配合基因链方向看，负值表示在转录起始方向的上游。
