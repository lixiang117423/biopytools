# 基因信息+序列合并表 | Gene Info + Sequence Merged Table

一句话理解：**把每个基因的「户口信息 + 全长 DNA + 上下游区间 + CDS + 蛋白」打包成一张表**，同时导出 gene.fa / cds.fa / pep.fa，方便下游一次性取用。

## 功能概述 | Overview { #overview }

- 输出一张「每行一个转录本」的 TSV：坐标 + Gene_DNA + Region + CDS + Protein
- 同时写出 gene.fa / region.fa（工具切片）与 cds.fa / pep.fa（gffread 产出）
- 支持 GFF3（含 .gz），自动识别 gene / mRNA / transcript / CDS 特征
- 可选每基因只保留最长转录本、按基因 DNA 长度过滤
- 上游/下游侧翼长度可调，区间按链方向取 5'→3'
- gffread 路径自动检测（conda 环境自动识别），也可 `--gffread` 覆盖

## 快速开始 | Quick Start { #quick-start }

```bash
biopytools gene-table -g genome.fa -f input.gff -o out.tsv
```

最小输入：基因组 FASTA + GFF3 注释 + 输出路径。

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语<br>Term | 通俗理解 |
|---|---|
| 合并表 | 把分散在基因组/注释里的信息拼到一张表，每行一个转录本 |
| Gene_DNA | 基因全长 DNA（含内含子、UTR） |
| Region | 上游 + 基因 + 下游 拼起来的一段，用来找启动子/调控区 |
| CDS / Protein | 编码序列 / 翻译出的蛋白序列 |
| 最长转录本 | 一个基因可能剪出多个转录本，取最长的那个做代表 |
| gffread | 一个从 GFF+基因组抽取 CDS/蛋白的经典小工具 |

## 输入 | Input { #input }

- **基因组** `-g`：FASTA（序列 ID 须与 GFF 第 1 列一致）
- **GFF3** `-f`：含 gene / mRNA / transcript / CDS 特征，支持 `.gz` 压缩
- 输出路径 `-o` 可以是文件（`out.tsv`）或目录（自动命名 `{prefix}.gene_table.tsv`）

## 参数说明 | Parameters { #parameters }

### 输出与样本列 | Output & sample

**通俗理解|In plain words:** `-o` 给文件就用它当表名、目录就自动生成 `{prefix}.gene_table.tsv`；`--prefix` 同时决定各 .fa 文件名和表里的 Sample 列，默认从 GFF 文件名推断。**一般不用动。**

相关参数：`-o/--output`（必需）、`--prefix`。

### 转录本与过滤 | Transcript & filtering

**通俗理解|In plain words:** `--longest-only` 打开后每个基因只留最长转录本（表行数变少、每基因一行）；`--min-length` 过滤基因 DNA 过短的（0=不过滤），用来甩掉碎片注释。

相关参数：`--longest-only`、`--min-length`（默认 0）。

### 侧翼区间 | Flanking region

**通俗理解|In plain words:** `--upstream`/`--downstream` 控制 Region 列在基因两侧各取多长（bp）。取多少取决于下游用途：找启动子一般上游 2000–3000、下游 500–1000。**默认 3000/1000 是启动子分析的常用值。**

相关参数：`--upstream`（默认 3000）、`--downstream`（默认 1000）。

### 工具路径 | Tool path

**通俗理解|In plain words:** `--gffread` 显式指定 gffread 路径（默认自动检测，含 conda 环境识别）。**一般不用动，除非 gffread 不在默认位置。**

相关参数：`--gffread`。

## 分析流程 | Pipeline { #pipeline }

```text
基因组 FASTA + GFF3
    │
    ▼
解析 GFF(基因/转录本/CDS 长度)
    │
    ▼
切片:基因全长 DNA + 上下游区间(链定向)
    │
    ▼
gffread 出 CDS / 蛋白
    │
    ▼
按转录本合并写 TSV + gene.fa/region.fa
```

## 输出 | Output { #output }

```text
输出目录/
├── {prefix}.gene_table.tsv    # 合并表(主结果,每行一个转录本)
├── {prefix}.gene.fa           # 基因全长 DNA
├── {prefix}.region.fa         # 上游+基因+下游区间
├── {prefix}.cds.fa            # CDS 序列(gffread)
└── {prefix}.pep.fa            # 蛋白序列(gffread)
```

TSV 列（13 列）：`Sample`、`Gene_ID`、`Transcript_ID`、`Chromosome`、`Strand`、`Gene_Start`、`Gene_End`、`Transcript_Start`、`Transcript_End`、`Gene_DNA`、`Region`、`CDS`、`Protein`。缺失的序列填 `NA`。

## 结果解读 | Interpreting Results { #interpreting-results }

- 表里一行 = 一个转录本；同一个基因有多个转录本时会占多行（除非 `--longest-only`）
- `Gene_DNA`/`Region`/`CDS`/`Protein` 为 `NA` 表示该转录本没拿到对应序列（如无 CDS 的非编码转录本，或染色体 ID 对不上）
- 日志末尾会打印各列 `NA` 计数：`CDS`/`Protein` 的 NA 数 = 非蛋白编码转录本数，属正常
- `Sample` 列固定为 prefix，方便多物种/多样本拼表时区分来源

## 参数选择建议 | Parameter Guidance { #parameter-guidance }

- 只要蛋白编码基因的代表序列：`--longest-only`
- 做启动子/调控区分析：保持默认 `--upstream 3000 --downstream 1000`（或按需调整）
- 清洗碎片注释：设 `--min-length`（如 300、500）
- 输出目录而非文件：`-o out_dir/`，自动命名多份 .fa 与 TSV

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--genome, -g` | 必填 |  | 基因组 FASTA｜Genome FASTA |
| `--gff, -f` | 必填 |  | GFF3 注释(支持 .gz)｜GFF3 annotation (gz supported) |
| `--output, -o` | 必填 |  | 输出表路径(或目录)｜Output table path (or directory) |
| `--prefix` | — |  | 输出前缀 + Sample 列(默认取 GFF 文件名)｜Output prefix + Sample column |
| `--longest-only` | — |  | 每基因仅保留最长转录本(默认全部)｜Keep only longest transcript per gene |
| `--min-length` | `0` | int | 基因 DNA 最小长度过滤｜Min gene-DNA length filter |
| `--upstream` | `3000` | int | 上游侧翼长度｜Upstream flank length |
| `--downstream` | `1000` | int | 下游侧翼长度｜Downstream flank length |
| `--gffread` | — |  | gffread 路径(默认自动检测)｜gffread path (auto-detected) |
| `--log-file` | — |  | 日志文件｜Log file |
| `--log-level` | `INFO` | DEBUG/INFO/WARNING/ERROR/CRITICAL | 日志级别｜Log level |
| `--verbose, -v` | — |  | 详细日志｜Verbose |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-g, --genome` | 必填 |  | 基因组 FASTA｜Genome FASTA |
| `-f, --gff` | 必填 |  | GFF3 注释(支持 .gz)｜GFF3 annotation (gz supported) |
| `-o, --output` | 必填 |  | 输出表路径(或目录)｜Output table path (or directory) |
| `--prefix` | — |  | 输出文件前缀 + Sample 列(默认取 GFF 文件名)｜Output prefix + Sample column |
| `--longest-only` | — | store_true | 每基因仅保留最长转录本(默认全部)｜Keep only longest transcript per gene |
| `--transcript-types` | `['mRNA', 'transcript']` |  | 视为转录本的 feature 类型｜Feature types treated as transcripts |
| `--gene-type` | `gene` |  | 基因 feature 类型｜Gene feature type |
| `--min-length` | `0` | int | 基因 DNA 最小长度过滤(0=不过滤)｜Min gene-DNA length filter |
| `--upstream` | `3000` | int | 上游侧翼长度(默认3000)｜Upstream flank length |
| `--downstream` | `1000` | int | 下游侧翼长度(默认1000)｜Downstream flank length |
| `--gffread` | — |  | gffread 路径(默认自动检测)｜gffread path (auto-detected by default) |
| `--log-file` | — |  | 日志文件｜Log file |
| `--log-level` | `INFO` | DEBUG/INFO/WARNING/ERROR/CRITICAL | 日志级别｜Log level |
| `-v, --verbose` | — | store_true | 详细日志｜Verbose |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies { #dependencies }

- gffread（默认 `~/.local/bin/gffread`，环境变量 `GFFREAD_PATH` 可覆盖；conda 环境自动识别）
- Python 3（无其他必需 Python 依赖）

## 常见问题 | FAQ { #faq }

**Q1：会断点续传吗？**
不会。单遍解析即出结果，重跑直接覆盖输出。

**Q2：CDS/Protein 列出现很多 NA？**
通常是这些转录本没有 CDS 特征（非编码 RNA），或基因组染色体 ID 与 GFF 第 1 列不一致。检查 GFF 里是否有 CDS 行、ID 是否匹配。

**Q3：`-o` 到底当文件还是目录？**
以 `/` 结尾或已存在的目录会当目录（自动命名）；带扩展名的路径当表文件路径。要精确控制就用 `--prefix` + `-o` 目录。

**Q4：gffread 找不到？**
默认找 `~/.local/bin/gffread`（可设 `GFFREAD_PATH`），或直接用 `--gffread /path/to/gffread` 指定。
