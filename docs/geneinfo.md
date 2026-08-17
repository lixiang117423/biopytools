# geneinfo - GFF3 基因转录本信息提取 | GFF3 Gene/Transcript Information Extraction

一句话理解：**把一份 GFF3 注释文件里「每个基因有哪些转录本、各自在染色体什么位置」的零散信息，整理成一张一行一个转录本的表格（TSV）**，方便直接用 Excel 或脚本做统计和下游分析。

## 功能概述 | Overview

- 从 GFF3 文件（或一个目录里的多个 GFF3 文件）中提取基因与转录本的整合信息
- 输出一张制表符分隔的 TSV，一行对应一个转录本，同时带上它所属基因的坐标与链方向
- 支持单文件与批量（目录）两种输入方式，目录下自动扫描 .gff3/.gff 及 .gz 压缩版本
- 转录本的父基因找不到时（孤儿转录本）不报错，基因坐标记 NA，照常输出
- 附带生成一份文本总结报告，统计样品数、转录本数、涉及基因数、孤儿转录本比例等
- 纯 Python 实现，无需任何外部软件或 conda 环境

## 快速开始 | Quick Start

```bash
biopytools geneinfo -i genome.gff3 -o gene_info.tsv
```

## 零基础概念速览 | Concepts in plain words

不熟悉生信术语的话，先看这张表：

| 术语 | 通俗理解 |
|------|----------|
| GFF3 | 基因组注释的「标准账本」格式，每一行用制表符分隔 9 列，记录一条特征（基因、转录本、外显子等） |
| 特征类型 | 第 3 列，表示这一行是什么东西：gene（基因）、mRNA/transcript（转录本）、exon（外显子）、CDS（编码区）等 |
| 基因 gene | 基因组上一段能发挥功能的区域，像一本书 |
| 转录本 transcript | 同一个基因可能「读出」多个版本（可变剪接），每个版本是一个转录本，像一本书的不同修订版 |
| ID / Parent | 特征之间的「父子关系」：转录本用 Parent=基因ID 挂在基因下面；本工具靠这两个属性把基因和转录本对上号 |
| 孤儿转录本 | 转录本的 Parent 指向的基因在文件里找不到，像「认领不到家长的孩子」 |
| 链方向 Strand | 基因在 DNA 哪条链上读（+ 或 -），像单行道的行车方向 |
| TSV | 制表符分隔的表格文件，Excel 可直接打开 |

## 输入 | Input

### GFF3 文件

标准 GFF3 格式，支持普通与 gzip 压缩版本，扩展名识别 .gff3.gz、.gff.gz、.gff3、.gff。-i 可传单个文件，也可传一个目录（自动扫描该目录下所有 GFF 文件并批量提取）。

关键要求：转录本行（如 mRNA/transcript）必须有 ID 和 Parent 两个属性，Parent 的值要与某个基因行的 ID 对应，这样才能正确关联基因信息。

```text
##gff-version 3
chr1    source  gene    1000    5000    .    +    .    ID=gene-001
chr1    source  mRNA    1000    5000    .    +    .    ID=mRNA-001;Parent=gene-001
chr1    source  CDS     1200    4800    .    +    0    ID=CDS-001;Parent=mRNA-001
```

## 参数说明 | Parameters

### 输入输出 | Input & output

**通俗理解|In plain words:** 这是两个必填项。-i 告诉程序「注释文件在哪」（单个文件或一个目录都行），-o 告诉程序「结果表格写到哪」。一般不用动其他参数，这两个填对就能跑。

### 特征类型 | Feature types

**通俗理解|In plain words:** 决定「什么样的行算基因、什么样的行算转录本」。默认基因类型是 gene（第 3 列等于 gene 的行当成基因）、转录本类型是 mRNA 和 transcript 两种。绝大多数 GFF3 用默认值即可；只有当你的注释文件用了别的类型名（例如基因叫 protein_coding_gene，或转录本只有 lnc_RNA）时才需要改。

### 日志与调试 | Logging & debugging

**通俗理解|In plain words:** 控制程序「说多少话」。-v 说得更多（调试用），--quiet 只说错误，--log-file 把日志写到指定文件，--log-level 设日志级别。一般不用动。

### 高级选项 | Advanced

**通俗理解|In plain words:** --dry-run 是「预演模式」，只做参数校验不真正执行；-t/--threads 是线程数。注意：本工具是单遍顺序扫描，线程数参数目前不影响实际处理速度，保持默认即可。

## 分析流程 | Pipeline

**通俗理解|In plain words:** 整个过程就是「扫两遍同一份文件」——第一遍先记下所有基因的位置，第二遍再给每个转录本补上它父基因的位置信息。

```text
输入 GFF3 文件或目录
    |
    v
第一遍: 收集所有基因 (type == gene-type) 的 ID/染色体/坐标/链
    |
    v
第二遍: 收集所有转录本 (type in transcript-types)，按 Parent 关联父基因
    |   (找不到父基因的记为孤儿转录本，基因坐标填 NA)
    v
写入 TSV 结果表 + 生成文本总结报告
```

## 输出 | Output

输出目录（-o 所在目录）下会生成三个文件：

```text
<输出目录>/
├── gene_info.tsv                    # 主结果表（文件名由 -o 指定）
├── gff_extraction_summary.txt       # 提取统计总结报告
└── gff_processing_<时间戳>.log      # 运行日志
```

主结果表 TSV 的列（制表符分隔，9 列）：

```text
Sample    Gene_ID    Transcript_ID    Chromosome    Strand    Gene_Start    Gene_End    Transcript_Start    Transcript_End
```

- Sample：样品名（单文件取文件名去掉扩展名；批量模式取各文件名）
- Gene_ID / Transcript_ID：基因与转录本的原始 ID
- Chromosome / Strand：染色体与链方向
- Gene_Start / Gene_End：父基因坐标；孤儿转录本这两列为 NA
- Transcript_Start / Transcript_End：转录本自身坐标

## 结果解读 | Interpreting Results

**通俗理解|In plain words:** 主表就是「每个转录本的户口信息」，一行一个转录本，看它属于哪个基因、在染色体哪一段。

- 看总转录本数与涉及基因数：两者越接近，说明平均每个基因只有一个转录本；差距大说明可变剪接普遍
- 看孤儿转录本比例（总结报告里）：比例很高（比如超过 10%）通常意味着注释文件的 Parent 关系不全或类型名不匹配，需要回头检查输入 GFF3 或调整转录本类型参数
- Gene_Start/Gene_End 为 NA 的行就是孤儿转录本，想排查时可单独筛出这些 Transcript_ID
- 平均每个基因转录本数（报告中）可快速判断注释复杂度

## 参数选择建议 | Parameter Guidance

- 默认跑法：biopytools geneinfo -i input.gff3 -o out.tsv 即可，不需要动任何其他参数
- 批量：-i 直接传目录，程序自动把所有 GFF 文件合并成一张大表，Sample 列区分来源
- 类型名不标准：先 head 看一下 GFF3 第 3 列实际用了哪些值，再用 --gene-type / --transcript-types 指定
- 只想快速看能提取多少：可先 --dry-run 校验参数，再用 -v 观察过程

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input, -i` | 必填 | Path | 输入GFF3文件或目录路径｜Input GFF3 file or directory path |
| `--output, -o` | 必填 | Path | 输出TSV文件路径｜Output TSV file path |
| `--gene-type` | `gene` |  | 基因特征类型｜Gene feature type |
| `--transcript-types` | `['mRNA', 'transcript']` |  | 转录本特征类型｜Transcript feature types |
| `--verbose, -v` | — |  | 增加输出详细程度｜Increase output verbosity |
| `--quiet` | — |  | 静默模式，仅输出错误信息｜Quiet mode, only output errors |
| `--log-file` | — | str | 日志文件路径｜Log file path |
| `--log-level` | `INFO` | DEBUG/INFO/WARNING/ERROR/CRITICAL | 日志级别｜Log level |
| `--dry-run` | — |  | 试运行模式，不实际执行｜Dry run mode, no actual execution |
| `--threads, -t` | `12` | int | 线程数｜Number of threads |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input, -i` | 必填 |  | 输入GFF3文件或目录路径｜Input GFF3 file or directory path |
| `--output, -o` | 必填 |  | 输出的TSV文件路径｜Output TSV file path |
| `--gene-type` | `gene` |  | 基因特征类型｜Gene feature type |
| `--transcript-types` | `['mRNA', 'transcript']` |  | 转录本特征类型列表｜Transcript feature types list |
| `-v, --verbose` | `0` | count | 增加输出详细程度｜Increase output verbosity |
| `--quiet` | — | store_true | 静默模式，仅输出错误信息｜Quiet mode, only output errors |
| `--log-file` | — | str | 日志文件路径｜Log file path |
| `--log-level` | `INFO` | DEBUG/INFO/WARNING/ERROR/CRITICAL | 日志级别｜Log level |
| `--dry-run` | — | store_true | 试运行模式，不实际执行｜Dry run mode, no actual execution |
| `-t, --threads` | `1` | int | 线程数｜Number of threads |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- Python 3（纯标准库实现）
- 无外部软件、无 conda 环境依赖

## 常见问题 | FAQ

**Q1：为什么结果里有些 Gene_Start/Gene_End 是 NA？**
那是「孤儿转录本」——它的 Parent 指向的基因在本文件里找不到。不是程序出错，是注释文件本身的父子关系不全。若比例很高，检查 --gene-type/--transcript-types 是否与你文件的类型名一致。

**Q2：-i 传目录时，会处理哪些文件？**
只处理扩展名为 .gff3、.gff、.gff3.gz、.gff.gz 的文件（大小写不敏感），按文件名排序处理。

**Q3：-t/--threads 能加速吗？**
不能。本工具是单遍顺序扫描，线程数参数不影响处理速度，保持默认即可。

**Q4：输入是 gzip 压缩的 GFF3 能直接读吗？**
能。程序按扩展名 .gz 自动以 gzip 方式解压读取。