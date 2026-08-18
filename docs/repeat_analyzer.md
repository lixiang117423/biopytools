# RepeatAnalyzer - 重复序列分析 | Repeat Sequence Analysis

一句话理解：**给一个基因组自动跑一套完整的重复序列分析**——从头建重复库(RepeatModeler)、补全长 LTR(LTR_FINDER/LTRharvest/LTR_retriever)、合并库、RepeatMasker 注释屏蔽、TEsorter 分类，一站式拿到「重复库 + 注释 + 分类」。

## 功能概述 | Overview { #overview }

- 6 工具串联：RepeatModeler → LTR_FINDER/LTRharvest/LTR_retriever → 库合并 → RepeatMasker → TEsorter
- 可跳过 RepeatModeler（`--skip-modeler`）或 LTR 分析（`--skip-ltr`），灵活适配不同需求
- 工具路径自动检测（PATH / 常见路径），也可用 `--xxx-path` 手动指定
- 生成分析汇总报告 `repeat_analysis_summary.txt`（基因组统计 + 各库 + 各输出文件清单）
- 注意：**无断点续传**，每次重跑从头再来（见 FAQ）

## 快速开始 | Quick Start { #quickstart }

```bash
biopytools repeat-analyzer -i genome.fasta -o repeat_results
```

最小输入：一个基因组 FASTA。默认跑完整流程（RepeatModeler + LTR + RepeatMasker + TEsorter）。

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语 | 通俗理解 |
|------|----------|
| 重复序列 | 基因组里反复出现的大量相似片段（主要是转座子），基因预测前要先把它们标出来 |
| 从头建库(de novo) | 不看参考、直接从这份基因组里「归纳」出重复序列的代表库 |
| RepeatModeler | 最常用的从头建库软件，自动识别并合并各类重复序列 |
| LTR 逆转座子 | 最大的一类转座子，两端有重复的「长末端重复」；需专门的 LTR_FINDER/LTRharvest 找全长拷贝 |
| RepeatMasker | 拿重复库去基因组里「比对着色」，标出每个重复在哪、屏蔽掉 |
| 屏蔽(masking) | 把重复区用小写字母标记(软屏蔽)，基因预测软件会跳过这些区 |
| TEsorter | 把重复库里的每条序列归类到具体的转座子家族（如 LTR/Copia、LTR/Gypsy） |

## 输入 | Input { #input }

### 基因组 FASTA

必需，标准 FASTA 格式（可含多序列）：

```text
>chr1
ATGCGTACGTACGTAGCTA...
>chr2
TTGCAAGCTAGCATCGATC...
```

## 参数说明 | Parameters { #parameters }

### 流程控制 | Pipeline control

**通俗理解|In plain words:** `--skip-modeler` 跳过最耗时的从头建库（你已有库、或只想跑 LTR 时用）；`--skip-ltr` 跳过 LTR 结构分析（不关心全长 LTR 时用）。两者都开则只跑 RepeatMasker + TEsorter（需已有库）。**默认都不跳过，跑完整流程。**

本组参数：`--skip-modeler`、`--skip-ltr`（均为开关，默认关）。

### 运行参数 | Runtime

**通俗理解|In plain words:** `-t` 线程数，RepeatModeler/RepeatMasker/TEsorter 都很吃核，**越大越快但内存占用也越大**；CLI 默认 12。

本组参数：`-t/--threads`（默认 12）。

### 工具路径 | Tool paths

**通俗理解|In plain words:** 六个软件的路径。程序会先从 PATH 自动检测，找不到时用这些参数手动指定。**软件都装好且能 `which` 到时不用动**；装在非标准位置时才逐个指定。注意 `--ltrharvest-path` 默认是复合命令 `gt ltrharvest`。

本组参数：`--repeatmodeler-path`（默认 RepeatModeler）、`--ltr-finder-path`（默认 ltr_finder）、`--ltrharvest-path`（默认 gt ltrharvest）、`--ltr-retriever-path`（默认 LTR_retriever）、`--repeatmasker-path`（默认 RepeatMasker）、`--tesorter-path`（默认 TEsorter）。

## 分析流程 | Pipeline { #pipeline }

**通俗理解|In plain words:** 先统计基因组，再从头建库、补 LTR，把两个库合并成最终库，用它做屏蔽和分类。

```text
基因组 FASTA
    |
    v
1. 依赖检查 + 基因组统计(序列数/长度/N50/GC)
    |
    v
2. RepeatModeler 从头建库(可 --skip-modeler) -> {base}_db-families.fa
    |
    v
3. LTR 分析(可 --skip-ltr):
     LTR_FINDER + LTRharvest -> LTR_retriever 筛选 -> LTR 库
    |
    v
4. 合并库(单库直接用, 多库合并) -> {base}_repeat_combined.fa
    |
    v
5. RepeatMasker 屏蔽 -> repeatmasker_output/
    |
    v
6. TEsorter 分类   -> tesorter_output/
    |
    v
7. 生成汇总 repeat_analysis_summary.txt
```

其中 `{base}` 为 `<基因组名>_repeat`（如 `genome.fasta` → `genome_repeat`）。

## 输出 | Output { #output }

```text
<输出目录>/
├── {base}_db*                        # BuildDatabase 数据库文件(二进制)
├── {base}_db-families.fa             # RepeatModeler 从头库
├── {base}_ltrfinder.scn              # LTR_FINDER 候选
├── {base}_ltrharvest.out             # LTRharvest 候选(含 .inner/.gff3)
├── {base}_repeat_combined.fa         # 合并后的最终重复库
├── repeatmasker_output/
│   ├── <基因组>.out                  # RepeatMasker 注释(核心)
│   ├── <基因组>.out.gff              # 注释 GFF
│   ├── <基因组>.masked               # 屏蔽后基因组(软屏蔽小写)
│   └── <基因组>.tbl                  # 统计表
├── tesorter_output/
│   ├── {base}_tesorter.cls.lib       # 分类后的库
│   ├── {base}_tesorter.cls.tsv       # 分类结果表(核心)
│   ├── {base}_tesorter.rexdb.cls.lib # rexdb-plant 分类库
│   └── {base}_tesorter.rexdb.cls.tsv # rexdb-plant 分类表
└── repeat_analysis_summary.txt       # 汇总报告 + repeat_analysis.log 日志
```

## 结果解读 | Interpreting Results { #interpreting }

### 1. 合并库（`{base}_repeat_combined.fa`）

**通俗理解|In plain words:** 最终「重复序列代表名单」，RepeatModeler 的从头库 + LTR 全长库合并而成，可直接复用于其它基因组屏蔽。

- 序列头里带 `#RepeatModeler` / `#LTR` 来源标识，方便追溯

### 2. RepeatMasker 注释（`repeatmasker_output/`）

- `<基因组>.out`：每个重复的位置、类型、分歧度；`<基因组>.tbl`：按类型汇总的占比表
- **总重复占比**是核心指标（多数植物 30%-80%）；异常低提示建库漏检，异常高提示假阳性
- `<基因组>.masked`：软屏蔽基因组，直接用于下游基因预测

### 3. TEsorter 分类表（`tesorter_output/*.cls.tsv`）

**通俗理解|In plain words:** 把库里的每条序列归到具体转座子家族（如 LTR/Copia、LTR/Gypsy、TIR/Mutator），TSV 表里一列是序列 ID、一列是家族名。

- 看各家族数量分布，判断该基因组转座子谱系构成；TSV 可导入 Excel/R 做饼图

### 4. 汇总报告（`repeat_analysis_summary.txt`）

基因组统计（序列数/总长/N50/GC）+ 各库路径 + RepeatMasker/TEsorter 输出清单，一次看全流程产物。

## 参数选择建议 | Parameter Guidance { #guidance }

- **`--skip-modeler`**：已有现成重复库、或基因组巨大只关心 LTR 时开启，可省下 RepeatModeler 的大量时间
- **`--skip-ltr`**：不关心全长 LTR 结构（如动物基因组 LTR 较少）时开启
- **两者都开**：直接用外部库跑 RepeatMasker + TEsorter（需确保有可用库文件）
- **`--threads`**：RepeatModeler 是大头，建议 24-64；内存不足时适当降低
- **工具路径**：默认从 PATH 检测；conda 环境未激活时，用 `--xxx-path` 指向绝对路径（如 `~/miniforge3/envs/<env>/bin/RepeatModeler`）

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input, -i` | 必填 |  | FASTA｜Input genome FASTA file path |
| `--output, -o` | 必填 |  | 输出目录｜Output directory |
| `--threads, -t` | `12` | int | (: 88)｜Number of threads (default: 88) |
| `--skip-modeler` | — |  | RepeatModeler｜Skip RepeatModeler step |
| `--skip-ltr` | — |  | LTR｜Skip LTR analysis step |
| `--repeatmodeler-path` | `RepeatModeler` |  | RepeatModeler (: RepeatModeler)｜RepeatModeler program path (default: RepeatModeler) |
| `--ltr-finder-path` | `ltr_finder` |  | LTR_FINDER (: ltr_finder)｜LTR_FINDER program path (default: ltr_finder) |
| `--ltrharvest-path` | `gt ltrharvest` |  | LTRharvest (: gt ltrharvest)｜LTRharvest program path (default: gt ltrharvest) |
| `--ltr-retriever-path` | `LTR_retriever` |  | LTR_retriever (: LTR_retriever)｜LTR_retriever program path (default: LTR_retriever) |
| `--repeatmasker-path` | `RepeatMasker` |  | RepeatMasker (: RepeatMasker)｜RepeatMasker program path (default: RepeatMasker) |
| `--tesorter-path` | `TEsorter` |  | TEsorter (: TEsorter)｜TEsorter program path (default: TEsorter) |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 输入基因组FASTA文件路径｜Input genome FASTA file path |
| `-o, --output` | `./repeat_output` |  | 输出目录｜Output directory |
| `-t, --threads` | `88` | int | 线程数｜Number of threads |
| `--skip-modeler` | — | store_true | 跳过RepeatModeler步骤｜Skip RepeatModeler step |
| `--skip-ltr` | — | store_true | 跳过LTR分析步骤｜Skip LTR analysis step |
| `--repeatmodeler-path` | — |  | RepeatModeler程序路径｜RepeatModeler program path |
| `--ltr-finder-path` | — |  | LTR_FINDER程序路径｜LTR_FINDER program path |
| `--ltrharvest-path` | — |  | GenomeTools gt二进制路径｜GenomeTools gt binary path |
| `--ltr-retriever-path` | — |  | LTR_retriever程序路径｜LTR_retriever program path |
| `--repeatmasker-path` | — |  | RepeatMasker程序路径｜RepeatMasker program path |
| `--tesorter-path` | — |  | TEsorter程序路径｜TEsorter program path |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies { #dependencies }

- RepeatModeler + BuildDatabase
- LTR_FINDER（`ltr_finder`）
- LTRharvest（GenomeTools `gt ltrharvest`）
- LTR_retriever
- RepeatMasker
- TEsorter（含 rexdb-plant 数据库）
- 以上工具需在 PATH 可找到，或用 `--xxx-path` 指定；模块本身不绑定特定 conda 环境名

## 常见问题 | FAQ { #faq }

**Q1：重跑会复用之前的结果吗？**
不会。本模块**无断点续传**，每次运行都从头执行所有步骤。重跑前请确认输出目录可覆盖，或换新目录。

**Q2：报「缺少依赖软件 Missing dependencies」？**
某工具不在 PATH 或未安装。用 `--xxx-path` 指定对应工具的完整路径（如 `--repeatmodeler-path ~/miniforge3/envs/<env>/bin/RepeatModeler`）。

**Q3：`--ltrharvest-path` 为什么默认是 `gt ltrharvest`？**
LTRharvest 是 GenomeTools 的一个子命令，需通过 `gt` 调用。指定路径时也写完整命令形式，如 `--ltrharvest-path '~/.../bin/gt ltrharvest'`。

**Q4：RepeatModeler 返回非零退出码但库文件生成了？**
RepeatModeler 偶发「文件已生成但返回非零」，本模块会以库文件是否存在为准，存在即视为成功，可放心继续。

**Q5：线程数默认到底是多少？**
通过 `biopytools repeat-analyzer` 调用时默认 12 线程；直接 `python -m biopytools.repeat_analyzer` 时默认 88。建议显式传 `-t` 避免歧义。
