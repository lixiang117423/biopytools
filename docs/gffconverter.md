# renamegff - GFF 基因 ID 整理 | GFF Gene ID Organization

一句话理解：**给一份 GFF 注释里的每个基因重新分配一套规范的编号 ID（如 OV53_Ov01G000010），并把它的转录本、外显子、CDS 等子特征的 ID 也一并改写成挂靠新基因 ID 的格式**，让注释 ID 整齐、连续、可预测。

## 功能概述 | Overview

- 按染色体给基因重新编号，每条染色体从起始编号开始、按步长递增（默认 10、20、30……）
- 新基因 ID 形如 物种名称_物种缩写+染色体号G+六位编号，例如 OV53_Ov01G000010
- 转录本、exon、CDS、UTR 等子特征的 ID 随基因新 ID 一并改写并保持父子关系
- 检测同一坐标的 gene 与非编码 RNA（lncRNA/tRNA/rRNA 等）重复，自动合并（给基因加 Note、跳过多余 ncRNA 行）
- 纯 Python 实现，无需外部软件或 conda 环境

## 快速开始 | Quick Start

```bash
biopytools renamegff -i input.gff -o output.gff -s OV53 -p Ov
```

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| 物种名称 species-name | 项目的全名或物种名（如 OV53、Ecoli_K12），会出现在每个基因 ID 开头 |
| 物种缩写 species-prefix | 物种的短名（如 Ov、Ec），拼在染色体号前面 |
| 起始编号 start-num | 每条染色体第一个基因的编号（默认 10） |
| 步长 step | 相邻基因编号的间隔（默认 10，即 10、20、30） |
| 染色体号 chr_num | 从染色体名里提取的数字，补成两位（chr1 变 01） |
| 子特征 child feature | 挂在转录本下面的 exon、CDS、UTR 等，ID 由父转录本 ID 派生 |

## 输入 | Input

### GFF 文件

标准 GFF 格式（.gff 或 .gff3）。要求基因行有 ID，转录本/子特征行有 Parent。工具按染色体名提取数字编号（支持 Chr1、chr1、OV12_、LG1、scaffold1 等格式，提取不到则默认 01）。

```text
##gff-version 3
chr1    source  gene    1000    5000    .    +    .    ID=gene-001;Name=hypothetical_protein
chr1    source  mRNA    1000    5000    .    +    .    ID=mRNA-001;Parent=gene-001
chr1    source  CDS     1200    4800    .    +    0    ID=CDS-001;Parent=mRNA-001
```

## 参数说明 | Parameters

### 输入输出与物种信息 | Required

**通俗理解|In plain words:** -i 输入、-o 输出、-s 物种名称、-p 物种缩写四个必填。-s 和 -p 一起决定新基因 ID 的「长相」（物种名称_物种缩写+染色体号G+编号），务必按项目约定填对。

### 编号参数 | Numbering

**通俗理解|In plain words:** --start-num 是每条染色体第一个基因的起始编号，--step 是相邻基因的编号间隔。调大 --step 给编号之间留更多空位（方便以后插入新基因），调小则编号更紧凑。默认都是 10，一般不用动；除非你有特定的编号起始要求。

### 运行与调试 | Runtime & debugging

**通俗理解|In plain words:** -t/--threads 线程数、-v/--verbose 详细模式、--keep-intermediate 保留中间文件这三个参数目前已被接收但转换逻辑中未实际使用（转换是单线程顺序处理）。--show-sample N 只打印前 N 个基因的原始 ID 预览后退出，不真正执行转换，可用来先看看输入长什么样。

## 分析流程 | Pipeline

**通俗理解|In plain words:** 解析 GFF → 按染色体给基因重新编号 → 建立新旧 ID 映射 → 按新 ID 重写每一行，途中顺便合并同坐标重复的 gene 与 ncRNA。

```text
输入 GFF 文件
    |
    v
步骤1: 解析 GFF (保留头部注释, 校验 9 列格式与坐标)
    |
    v
步骤2: 按染色体分组基因, 从 start-num 起、按 step 递增生成新基因 ID
    |
    v
步骤3: 建立完整 ID 映射 (基因 -> 转录本 -> 子特征), 检测并合并同坐标的 gene/ncRNA 重复
    |
    v
步骤4: 按新 ID 重写所有特征并输出 GFF 文件
```

## 输出 | Output

```text
<输出目录>/
├── output.gff                    # 整理后的 GFF 文件 (文件名由 -o 指定)
└── gff_conversion.log            # 运行日志
```

### 新 ID 命名规则

基因 ID = 物种名称_物种缩写+染色体号G+六位编号，例如物种名称 OV53、缩写 Ov、1 号染色体、起始编号 10：

```text
OV53_Ov01G000010      # 1 号染色体第 1 个基因 (编号 10)
OV53_Ov01G000020      # 第 2 个基因 (编号 20)
OV53_Ov02G000010      # 2 号染色体第 1 个基因
```

子特征 ID 由父 ID 派生：

```text
OV53_Ov01G000010.mRNA1        # 第 1 个 mRNA
OV53_Ov01G000010.exon1        # 第 1 个外显子
OV53_Ov01G000010.cds1         # 第 1 个 CDS
OV53_Ov01G000010.5utrp1       # 5' UTR
OV53_Ov01G000010.3utrp1       # 3' UTR
```

## 结果解读 | Interpreting Results

**通俗理解|In plain words:** 整理后的 GFF 结构不变、坐标不变，只有第 9 列的 ID/Parent/Name 换成了新编号。

- 检查基因 ID 是否连续且步长一致（10、20、30……），每个染色体独立从起始编号开始
- 检查子特征 Parent 是否指向新的转录本 ID，且每个转录本下 exon/CDS 编号连续
- 日志会打印输入输出的行数统计与基因/转录本/外显子计数，可对照确认无丢失
- 若同一坐标同时有 gene 和单个 ncRNA，ncRNA 行会被合并（基因加 Note），这是预期行为

## 参数选择建议 | Parameter Guidance

- 常规整理：biopytools renamegff -i in.gff -o out.gff -s <物种名> -p <缩写>，其余全默认
- 需要预留编号空位：调大 --step（如 20、50）；需要紧凑连续：--step 1、--start-num 1
- 想先看输入里有哪些基因 ID：--show-sample 10 预览前 10 个后退出
- 注意：本工具与 gff-renamer 的 ID 格式不同（本工具 ID 形如 OV53_Ov01G000010，带物种名称和 G），按你项目的规范选其一

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input, -i` | 必填 |  | GFF｜Input GFF file path |
| `--output, -o` | 必填 | Path | GFF｜Output GFF file path |
| `--species-name, -s` | 必填 | str | (: OV53)｜Species name (e.g., OV53) |
| `--species-prefix, -p` | 必填 | str | (: Ov)｜Species prefix (e.g., Ov) |
| `--start-num` | `10` | int | Starting number for gene numbering (default: 10)｜Starting number for gene numbering (default: 10) |
| `--step` | `10` | int | Step size for gene numbering (default: 10)｜Step size for gene numbering (default: 10) |
| `--threads, -t` | `12` | int | Number of threads (default: 88)｜Number of threads (default: 88) |
| `--verbose, -v` | — |  | Verbose output mode｜Verbose output mode |
| `--keep-intermediate` | — |  | Keep intermediate files｜Keep intermediate files |
| `--show-sample` | — | int | N｜Show N conversion samples and exit |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 输入GFF文件路径｜Input GFF file path |
| `-o, --output` | 必填 |  | 输出GFF文件路径｜Output GFF file path |
| `-s, --species-name` | 必填 |  | 物种名称 (如: OV53)｜Species name (e.g., OV53) |
| `-p, --species-prefix` | 必填 |  | 物种缩写 (如: Ov)｜Species prefix (e.g., Ov) |
| `--start-num` | `10` | int | 起始编号｜Starting number for gene numbering |
| `--step` | `10` | int | 编号步长｜Step size for gene numbering |
| `-t, --threads` | `88` | int | 线程数｜Number of threads |
| `-v, --verbose` | — | store_true | 详细输出模式｜Verbose output mode |
| `--keep-intermediate` | — | store_true | 保留中间文件｜Keep intermediate files |
| `--show-sample` | — | int | 显示N个转换示例然后退出｜Show N conversion samples and exit |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- Python 3（纯标准库实现）
- 无外部软件、无 conda 环境依赖

## 常见问题 | FAQ

**Q1：--threads / --keep-intermediate / --verbose 为什么不起作用？**
这三个参数目前被 CLI 接收并传入配置，但转换逻辑是单线程顺序处理、不写中间文件，所以暂未实际生效。属已知现象，按默认即可。

**Q2：CLI 帮助里提到的 parsed_features.tmp、id_mappings.txt、validation_report.txt 在哪？**
那些是早期帮助文档里写的中间文件，当前版本代码并不会真正生成它们（--keep-intermediate 未生效）。

**Q3：能断点续传吗？**
不能。单次运行、一次性写出，重跑覆盖输出文件。

**Q4：基因编号为什么从 10 开始？**
默认起始编号 10、步长 10，是为了给编号之间留空位便于后续插入。可用 --start-num 和 --step 自定义。

**Q5：和 gff-renamer 有什么区别？**
两者都是 GFF ID 规范化，但 ID 格式不同：本工具是 物种名称_物种缩写+染色体G+编号（OV53_Ov01G000010），gff-renamer 是 前缀_物种缩写+染色体g+编号（CDRT_Ov01g000010）且带 AGAT 清洗与去冗余。按项目约定选用。