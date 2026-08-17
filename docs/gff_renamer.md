# gff-renamer - GFF 文件 ID 规范化 | GFF File ID Standardization

一句话理解：**把一份 GFF 注释文件里 gene、mRNA、exon、CDS 等所有特征的 ID 统一重命名成一套规范编号**（例如 CDRT_Ov12g000010、CDRT_Ov12g000010.mRNA1、CDRT_Ov12g000010.exon1），让不同来源的注释能用同一套命名规则对接下游分析。

## 功能概述 | Overview

- 把 gene、mRNA、exon、CDS、intron、UTR、codon 等特征的 ID 统一重命名并保持父子关系
- 基因编号按染色体独立从 1 开始，每个基因递增 10（000010、000020……），保证 ID 稳定可排序
- 三种命名格式可选（standard / simple / compact），并可加版本号前缀（如 v1-）
- 默认先跑 AGAT 清洗（agat_convert_sp_gxf2gxf.pl），把输入 GFF 规整后再重命名，可用 --skip-clean 跳过
- 默认开启 prefer_mrna 去冗余：对含 mRNA 的基因丢弃冗余的 transcript 变体及其子特征
- 可选输出 mRNA 新旧 ID 映射表，方便追溯；支持自定义染色体映射文件
- 并行处理（默认 12 线程），大文件自动分块多进程

## 快速开始 | Quick Start

```bash
biopytools gff-renamer -i input.gff -o output.gff -p CDRT -s Ov
```

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| GFF / GFF3 | 基因组注释的标准文本格式，每行 9 列记录一个特征 |
| 特征 feature | 注释里的一行：gene（基因）、mRNA/transcript（转录本）、exon（外显子）、CDS（编码区）、intron（内含子）、UTR（非翻译区）等 |
| ID / Parent | 特征的「身份证号」和「家长」：子特征用 Parent 指向父特征，重命名时父子要一起换 |
| ID 前缀 prefix | 你给这套注释起的「代号」，例如物种项目名 CDRT、AGIS |
| 物种缩写 species | 物种的短名（如 Ov 水稻、Os 另一个物种），会拼进每个基因 ID 里 |
| 命名格式 | 基因 ID 的排版风格（standard/simple/compact），只影响 ID 长相，不影响含义 |
| AGAT | 一个注释格式整理工具，能把杂乱的 GFF 规整成标准 GFF3 |

## 输入 | Input

### GFF 文件

标准 GFF 格式（.gff 或 .gff3）。要求特征行有 ID 属性，转录本和子特征有 Parent 属性。工具会自动识别多种来源的 ID 风格：

- EGAPx 风格：ID=gene-XXX、ID=mRNA-XXX-R1
- BRAKER 风格：ID=gXXX.t1
- EviAnn 风格：ID=LOC_00000001、XXX-mRNA-1

```text
##gff-version 3
chr1    source  gene    1000    5000    .    +    .    ID=gene-001
chr1    source  mRNA    1000    5000    .    +    .    ID=mRNA-001;Parent=gene-001
chr1    source  CDS     1200    4800    .    +    0    ID=CDS-001;Parent=mRNA-001
```

### 染色体映射文件（可选）

两列文本（制表符或空格分隔）：第一列原始序列名，第二列标准化后的名称。用于把非标准染色体名（如 scaffold）规范成统一编号：

```text
chr1    Chr1
scaffold_1    Scf1
```

## 参数说明 | Parameters

### 输入输出与必填命名参数 | Required

**通俗理解|In plain words:** -i 输入、-o 输出、-p 前缀、-s 物种缩写四个必填。-p 和 -s 决定了生成 ID 的「长相」，填错会导致下游识别不了，所以一定要按项目约定填对。

### 命名格式与版本 | Naming & version

**通俗理解|In plain words:** --naming-format 决定基因 ID 的排版（standard 带物种缩写最完整，simple/compact 更短）；-v/--version 在所有 ID 最前面加 v版本- 前缀，适合同一套注释迭代多个版本时区分。一般用默认 standard 即可，除非下游对 ID 长度有要求。

### AGAT 清洗与去冗余 | Cleaning & dedup

**通俗理解|In plain words:** --skip-clean 跳过默认的 AGAT 清洗（输入已经很标准、想省时间时才用）；--prefer-mrna（默认开）会自动丢弃那些「有 mRNA 却还带一条冗余 transcript」的重复注释，通常能减少重复。一般不用动。

### 染色体与 UTR | Chromosome & UTR

**通俗理解|In plain words:** --chr-mapping 传入染色体映射文件来规范染色体名（scaffold/contig 较多时才用）；--include-utr 控制 UTR 特征的处理细节，默认关闭一般不用开。

### mRNA 映射输出 | mRNA mapping

**通俗理解|In plain words:** --output-mrna-mapping 输出一张「新 ID 对旧 ID」的对照表，方便追溯重命名前后对应关系；需要时才开。

### 并行线程 | Threads

**通俗理解|In plain words:** -t/--threads 控制并行进程数，默认 12，大文件可调高提速，但别超过机器核数。

## 分析流程 | Pipeline

**通俗理解|In plain words:** 先（可选）用 AGAT 把输入洗干净，再扫一遍文件建立「旧 ID → 新 ID」对照表，最后按对照表逐行替换并写出。

```text
输入 GFF 文件
    |
    v
步骤0: AGAT 清洗 (默认开, --skip-clean 跳过)
    |
    v
步骤1: 读取 GFF 全部行
    |
    v
步骤1.5: prefer_mrna 去冗余 (默认开, 丢弃冗余 transcript 及其子特征)
    |
    v
步骤2: 解析特征, 建立 旧ID -> 新ID 映射 (基因按染色体排序编号)
    |
    v
步骤3: 应用映射 (多进程并行替换 ID 与 Parent)
    |
    v
步骤4: 写出重命名后的 GFF (清理子特征 Name, 按位置重排)
    |
    v
步骤5: (可选) 生成 mRNA 新旧 ID 映射表
```

## 输出 | Output

```text
<输出目录>/
├── output.gff                          # 重命名后的 GFF 文件 (文件名由 -o 指定)
└── output_mrna_mapping.tsv             # (可选) mRNA 新旧 ID 映射表
```

### 新 ID 命名规则

基因 ID（standard 格式）：前缀_物种缩写+染色体号g+编号，例如前缀 CDRT、物种 Ov、12 号染色体第 1 个基因：

```text
CDRT_Ov12g000010      # 第 1 个基因 (编号 1*10 = 10)
CDRT_Ov12g000020      # 第 2 个基因
```

转录本与子特征 ID 由基因 ID 派生：

```text
CDRT_Ov12g000010.mRNA1            # 第 1 个 mRNA
CDRT_Ov12g000010.exon1            # 第 1 个外显子
CDRT_Ov12g000010.cds1             # 第 1 个 CDS
CDRT_Ov12g000010.five_prime_UTR1  # 5' UTR
CDRT_Ov12g000010.start1           # 起始密码子
```

### mRNA 映射表（--output-mrna-mapping 时生成）

每行：新 ID、旧 ID、以及原始第 9 列属性（按分号拆分），方便追溯。

## 结果解读 | Interpreting Results

**通俗理解|In plain words:** 重命名后的 GFF 应该保持「结构不变、ID 变规范」——特征数量、坐标、链方向都不该变，只有第 9 列的 ID/Parent/Name 变了。

- 用 diff 或 grep 核对：基因数、转录本数应与输入一致（除被 prefer_mrna 丢弃的冗余 transcript 外）
- 每个转录本的子特征 ID 应连续编号（.exon1、.exon2……），且 Parent 指向正确的转录本新 ID
- 日志会打印各类特征计数（genes/transcripts/exons/introns/CDS/codons/UTRs），可与输入对照确认无丢失
- prefer_mrna 去冗余会减少部分 transcript 及子特征，属预期行为；若不想丢弃，用 --no-prefer-mrna

## 参数选择建议 | Parameter Guidance

- 常规项目：biopytools gff-renamer -i in.gff -o out.gff -p <前缀> -s <物种缩写>，其余全默认
- 输入是 EGAPx/NCBI 注释（常有冗余 transcript）：保持 --prefer-mrna 默认开
- 输入已经很规范、想跳过 AGAT 提速：加 --skip-clean
- 需要追溯或后续做 gene 名转换：加 --output-mrna-mapping
- 染色体名混乱（scaffold/contig 混排）：准备映射文件并加 --chr-mapping
- 同一套注释要发多个版本：加 -v 1、-v 2 区分

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 输入GFF文件｜Input GFF file |
| `-o, --output` | 必填 |  | 输出GFF文件｜Output GFF file |
| `-p, --prefix` | 必填 |  | ID前缀｜ID prefix (e.g., CDRT, AGIS) |
| `-s, --species` | 必填 |  | 物种缩写｜Species abbreviation (e.g., Ov, Os) |
| `-t, --threads` | `12` | int | 并行线程数｜Number of parallel threads |
| `--output-mrna-mapping` | `False` |  | 输出mRNA映射文件｜Output mRNA mapping file |
| `--mrna-mapping-file` | — |  | mRNA映射文件路径（可选）｜mRNA mapping file path (optional) |
| `--chr-mapping` | — |  | 染色体映射文件路径｜Chromosome mapping file path |
| `--naming-format` | `standard` | standard/simple/compact | 命名格式｜Naming format (standard/simple/compact) |
| `--include-utr` | `False` |  | 包含UTR特征｜Include UTR features (five_prime_UTR, three_prime_UTR) |
| `--skip-clean` | `False` |  | 跳过AGAT清洗步骤｜Skip AGAT GFF cleaning step |
| `--prefer-mrna/--no-prefer-mrna` | `True` |  | 对含mRNA的基因丢弃冗余transcript(misc_RNA)变体及其子特征(默认开)｜Drop redundant transcript (misc_RNA) variants and their children from genes that have mRNA (default on) |
| `-v, --version` | — |  | 版本号前缀｜Version tag (e.g., 1, 1.2); 指定后在所有ID最前面加 v{version}- (如 -v 1 → v1-CDRT_Ov12g000010)｜when set, prepends v{version}- to every ID (e.g., -v 1 → v1-CDRT_Ov12g000010) |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 输入GFF文件｜Input GFF file |
| `-o, --output` | 必填 |  | 输出GFF文件｜Output GFF file |
| `-p, --prefix` | 必填 |  | ID前缀｜ID prefix (e.g., CDRT, AGIS) |
| `-s, --species` | 必填 |  | 物种缩写｜Species abbreviation (e.g., Ov, Os) |
| `-t, --threads` | `12` | int | 并行线程数｜Number of parallel threads |
| `--output-mrna-mapping` | — | store_true | 输出mRNA映射文件｜Output mRNA mapping file |
| `--mrna-mapping-file` | — |  | mRNA映射文件路径（可选，默认为输出文件名_mrna_mapping.tsv）｜mRNA mapping file path (optional, defaults to output_filename_mrna_mapping.tsv) |
| `--chr-mapping` | — |  | 染色体映射文件路径｜Chromosome mapping file path |
| `--naming-format` | `standard` | standard/simple/compact | 命名格式｜Naming format (standard/simple/compact, default: standard) |
| `--include-utr` | — | store_true | 包含UTR特征｜Include UTR features (five_prime_UTR, three_prime_UTR) |
| `--skip-clean` | — | store_true | 跳过AGAT清洗步骤｜Skip AGAT GFF cleaning step |
| `--prefer-mrna` | `True` | store_true | 默认开启:对含mRNA的基因丢弃冗余transcript(misc_RNA)变体及其子特征;仅含transcript的基因保留｜Default on: drop redundant transcript (misc_RNA) variants and their children from genes that have mRNA; transcript-only genes are kept |
| `--no-prefer-mrna` | — | store_false | 禁用prefer_mrna去冗余｜Disable prefer_mrna deduplication |
| `-v, --version` | — |  | 版本号前缀｜Version tag (e.g., 1, 1.2); 指定后在所有ID最前面加 v{version}- (如 -v 1 → v1-CDRT_Ov12g000010)｜when set, prepends v{version}- to every ID (e.g., -v 1 → v1-CDRT_Ov12g000010) |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- Python 3
- AGAT（agat_convert_sp_gxf2gxf.pl），默认路径 ~/miniforge3/envs/annot/bin/，conda 环境 annot
- 若加 --skip-clean 可跳过 AGAT，纯 Python 也能完成重命名

## 常见问题 | FAQ

**Q1：为什么输出里少了一些 transcript？**
默认开启了 prefer_mrna 去冗余：对含 mRNA 的基因，丢弃其冗余的 transcript（misc_RNA 变体）及这些转录本的子特征。仅含 transcript（无 mRNA）的基因保留。若想全部保留，用 --no-prefer-mrna。

**Q2：AGAT 清洗失败怎么办？**
确认 annot 环境里装了 agat_convert_sp_gxf2gxf.pl，或先用 --skip-clean 跳过清洗直接重命名。

**Q3：-v 版本号有什么限制？**
版本号不能包含空白字符（会破坏 GFF 属性 ID），例如 -v 1、-v 1.2 都合法。

**Q4：能断点续传吗？**
不能。本工具是单次运行、一次性写出输出，重跑会覆盖旧输出文件。

**Q5：基因编号为什么是 10、20、30 而不是 1、2、3？**
编号固定按 10 递增（000010、000020……），给后续插入新基因留出空位，这是设计约定，无法通过参数改步长（改步长请用 renamegff 工具）。