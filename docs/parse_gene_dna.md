# parse_gene_dna — 基因序列提取 | Gene DNA Sequence Extraction

一句话理解：**按注释把基因的 DNA 序列从基因组里「切」出来**——输入基因组 FASTA + GFF 注释，提取指定特征类型（默认 gene）的序列区间，输出成 FASTA，供下游比对、建库、序列分析使用。

## 功能概述 | Overview

- 从基因组 FASTA 按 GFF 注释提取基因序列
- 支持指定特征类型（默认 `gene`，也可提取 mRNA、exon 等）
- 支持最小长度过滤（`--min-length`）
- 处理负链基因（自动取反向互补），跳过坐标越界或不存在的染色体
- FASTA 头带基因名与坐标信息，便于溯源

## 快速开始 | Quick Start

```bash
biopytools parse-gene-dna -g genome.fasta -f annotation.gff -o genes.fasta
```

最小输入：一个基因组 FASTA + 一个 GFF 注释文件 + 一个输出 FASTA 路径。

## 零基础概念速览 | Concepts in plain words

不熟悉生信术语的话，先花两分钟看这张表：

| 术语 | 通俗理解 |
|------|----------|
| 特征类型 (feature type) | GFF 里「这一行描述什么」，如 gene（基因）、mRNA（转录本）、exon（外显子） |
| 坐标区间 | 基因在染色体上的「门牌号范围」，如 1000–5000 |
| 正链 / 负链 | 基因的「朝向」；负链基因要取反向互补才是正确的序列 |
| 反向互补 | 把 DNA 序列翻转并替换互补碱基（A-T、C-G），得到反方向那条链 |

## 输入 | Input

### 基因组 FASTA（-g, --genome）

标准 FASTA，序列名须与 GFF 第一列的 seqid 一致（不一致的会跳过并警告）。

### GFF 注释（-f, --gff）

标准 GFF 格式（9 列制表符分隔），提取第 3 列与 `--feature-type` 匹配的行。

```text
chr1    source    gene    1000    5000    .    +    .    ID=gene1;Name=example
chr1    source    gene    6000    9000    .    -    .    ID=gene2
```

## 参数说明 | Parameters

### 必需参数 | Required

**通俗理解|In plain words:** 三个必填项：基因组、GFF、输出文件。缺一不可。

### 提取控制 | Extraction controls

**通俗理解|In plain words:** `--feature-type` 决定「切什么」——默认 `gene` 提取整个基因区（含内含子）；想只取编码区就换成 `CDS` 或 `exon`。`--min-length` 是「短于多少就丢掉」的过滤线，0 表示不过滤，想排除过短片段就设一个阈值。`--line-width` 是 FASTA 每行的碱基数（纯排版，不影响序列本身），一般不用动。

### 运行选项 | Runtime options

**通俗理解|In plain words:** `--threads` 目前提取为纯 Python 单线程实现，该参数实际未用于并行，可忽略；`--verbose` 打开后会在日志里详细打印每个被跳过的基因及原因，排查问题时有用。

## 分析流程 | Pipeline

**通俗理解|In plain words:** 四步：读基因组 → 读 GFF 挑出目标特征 → 逐个切序列并过滤 → 写 FASTA。

```text
输入 基因组 FASTA + GFF
    │
    ▼
1. 加载基因组序列
    │
    ▼
2. 解析 GFF，筛选指定 feature type
    │
    ▼
3. 逐个提取序列区间（负链取反向互补），按最小长度过滤
    │
    ▼
4. 写出 FASTA（按指定行宽折行）
```

## 输出 | Output

```text
genes.fasta    # 提取出的基因序列（核心）
tmp/           # 日志相关（ExtractionLogger 输出到输出文件同目录）
```

### 关键文件说明 | Key files

- `genes.fasta`：提取出的基因 DNA 序列。FASTA 头格式为 `{gene_id} {gene_name} {seqid}:{start}-{end}({strand}) length={len}`，既含基因 ID 又有坐标，便于溯源

## 结果解读 | Interpreting Results

**通俗理解|In plain words:** 跑完看「总共提取了多少条」。每条序列的头里带着坐标和方向，负链基因已经是反向互补后的正确序列，可直接用于下游。

- 输出条数应等于 GFF 中目标 feature type 的数量减去被跳过的（染色体不存在、坐标越界、长度不足）
- 负链（strand 为 `-`）的基因，程序已自动取反向互补，序列方向与转录方向一致
- 用 `-v` 可看到每个被跳过基因的具体原因，便于核对是坐标问题还是长度过滤

## 参数选择建议 | Parameter Guidance

- **只要基因全长 DNA**：默认参数即可（`--feature-type gene`）
- **只要编码区**：`--feature-type CDS`
- **排除过短片段**：设 `--min-length`，如 `--min-length 100`
- **排查提取数量异常**：加 `-v` 看逐条跳过原因

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--genome, -g` | 必填 |  | 基因组FASTA文件路径｜Genome FASTA file path |
| `--gff, -f` | 必填 |  | GFF注释文件路径｜GFF annotation file path |
| `--output, -o` | 必填 | Path | 输出FASTA文件路径｜Output FASTA file path |
| `--feature-type` | `gene` |  | 要提取的特征类型｜Feature type to extract |
| `--min-length` | `0` | int | 最小基因长度过滤｜Minimum gene length filter |
| `--threads, -t` | `12` | int | 线程数｜Number of threads |
| `--line-width` | `60` | int | FASTA序列行宽度｜FASTA sequence line width |
| `--verbose, -v` | — |  | 显示详细信息｜Show verbose output |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-g, --genome` | 必填 |  | 基因组FASTA文件路径｜Genome FASTA file path |
| `-f, --gff` | 必填 |  | GFF注释文件路径｜GFF annotation file path |
| `-o, --output` | 必填 |  | 输出FASTA文件路径｜Output FASTA file path |
| `--feature-type` | `gene` |  | 要提取的特征类型｜Feature type to extract |
| `--min-length` | `0` | int | 最小基因长度过滤｜Minimum gene length filter |
| `-t, --threads` | `88` | int | 线程数｜Number of threads |
| `--line-width` | `60` | int | FASTA序列行宽度｜FASTA sequence line width |
| `-v, --verbose` | — | store_true | 显示详细信息｜Show verbose output |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- Python 3（标准库实现，无外部生信软件依赖）

## 常见问题 | FAQ

**Q1：提取出来的条数比 GFF 里的基因少，为什么？**
可能原因：某些基因的染色体 seqid 在基因组 FASTA 里不存在、坐标超出染色体长度、或长度低于 `--min-length`。加 `-v` 可逐条看到跳过原因。

**Q2：负链基因的序列方向对吗？**
对。负链基因程序会自动取反向互补，输出即为转录方向的正确序列。

**Q3：`--feature-type` 能填什么？**
任意 GFF 第三列出现的特征名（不区分大小写），如 `gene`、`mRNA`、`exon`、`CDS`。默认 `gene`。

**Q4：`--threads` 有用吗？**
当前提取逻辑为纯 Python 单线程实现，`--threads` 参数被接收但未实际并行，可忽略，不影响结果。

**Q5：GFF 里坐标是 1-based 吗？**
是的，GFF 坐标约定为 1-based（闭区间），程序按此直接切片，无需手动换算。
