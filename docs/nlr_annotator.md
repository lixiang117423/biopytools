# nlr_annotator — NLR 基因预测 | NLR Gene Prediction (NLR-Annotator)

一句话理解：**从 DNA 或 CDS 序列里把「抗病相关 NLR 基因」自动识别出来**——封装 NLR-Annotator（基于保守 motif 的 Java 工具），支持单个文件或整目录批量预测，并把各样本结果合并成一张汇总表。

## 功能概述 | Overview

- 封装 NLR-Annotator，从 DNA / CDS 序列预测 NLR（核苷酸结合亮氨酸重复）基因
- 支持单文件或目录批量处理（目录模式自动按样本名分目录）
- 输出清洗：自动加表头、去重并排序 motif 列表
- 可选输出 GFF / BED / motifs BED / motif 比对 FASTA
- 目录模式自动生成多样本汇总表 `nlr_annotator_summary.tsv`
- 断点续传：已完成样本自动跳过；`--merge-only` 可只合并已有结果补汇总

## 快速开始 | Quick Start

```bash
biopytools nlr-annotator -i genome.cds.fa -o output_dir/
```

最小输入：一个 DNA 或 CDS FASTA 文件（或目录）+ 一个输出目录。

## 零基础概念速览 | Concepts in plain words

不熟悉生信术语的话，先花两分钟看这张表：

| 术语 | 通俗理解 |
|------|----------|
| NLR 基因 | 植物免疫系统的「报警器」，识别病原并启动抗病反应 |
| motif | 一段保守的短序列「签名」，像身份证上的特征位点；NLR 基因靠多个特征 motif 来识别 |
| CDS | 基因的编码区核酸序列 |
| NBS-LRR | NLR 的经典结构：一个结合核苷酸的域 + 一串富亮氨酸重复 |
| motif 组合 | 几个特征 motif 以特定顺序、距离出现，构成判定 NLR 的依据 |

## 输入 | Input

### 输入文件或目录（-i, --input）

DNA 或 CDS 的 FASTA 文件；也可以是目录（目录模式按 `--sample-suffix` 匹配文件，默认 `*.fa`）。

```text
>gene1
ATGACGT...
>gene2
ATGGCTAG...
```

## 参数说明 | Parameters

### 必需与批处理 | Required & batch

**通俗理解|In plain words:** `-i` 是要预测的序列文件或目录；`-o` 是输出目录。目录模式下 `--sample-suffix` 控制匹配哪些文件（默认 `*.fa`），一般不用动。`--merge-only` 是「只合并已有结果、不再跑预测」，用于批量跑到一半被杀、只需补一张汇总表的情况。

### 工具路径 | Tool paths

**通俗理解|In plain words:** 这三件套是 NLR-Annotator 运行必需的外部资源——JAR 主程序、`mot.txt`（motif 定义）、`store.txt`（存储配置）。都有默认路径（`~/software/NLR-Annotator/` 下），装在不同位置才需指定。`--java-path` 是 Java 解释器，用 conda 环境的 Java 时填 `~/miniforge3/envs/xxx/bin/java`。

### 运行与距离参数 | Runtime & distance

**通俗理解|In plain words:** `--threads` 是并行线程数；`--num-seqs-per-thread` 是每线程处理的序列数（影响内存占用）。三个距离参数控制 motif 组合的判定边界：`--distance-within-motif-combination` 是「同一个组合内 motif 之间允许隔多远」，`--distance-for-elongating` 是「向外延伸多远」，`--distance-between-motif-combinations` 是「两个组合之间允许隔多远」。**这些都是 NLR-Annotator 的判定细则，默认值经原工具验证，一般不用动。**

### 可选输出 | Optional outputs

**通俗理解|In plain words:** 默认只出 TSV 结果表。需要坐标格式就开 `--output-gff` / `--output-bed`；想看每个预测用了哪些 motif 开 `--output-motifs`；想看 motif 的序列比对开 `--output-alignment`。**默认都不开，按需勾选。**

## 分析流程 | Pipeline

**通俗理解|In plain words:** 先收集输入文件，然后逐个样本调用 NLR-Annotator（Java），最后把各样本结果清洗并合并成汇总表。

```text
输入 DNA/CDS 文件或目录
    │
    ▼
1. 收集输入文件（单文件或目录批量）
    │
    ▼
2. 逐样本运行 NLR-Annotator（java -jar，断点续传跳过已完成）
    │
    ▼
3. 清洗输出：加表头、去重并排序 motif
    │
    ▼
4. 目录模式：合并所有样本结果 -> nlr_annotator_summary.tsv
```

## 输出 | Output

单文件模式：

```text
output/
├── {sample}.nlr_annotator.tsv      # 清洗后的 NLR 结果表（核心）
├── {sample}.gff                    # 仅 --output-gff
├── {sample}.bed                    # 仅 --output-bed
├── {sample}_motifs.bed             # 仅 --output-motifs
├── {sample}_alignment.fa           # 仅 --output-alignment
└── 99_logs/
    └── nlr_annotator.log
```

目录批处理模式：

```text
output/
├── {sample1}/
│   ├── {sample1}.nlr_annotator.tsv
│   └── 99_logs/{sample1}.nlr_annotator.log
├── {sample2}/
│   └── ...
└── nlr_annotator_summary.tsv       # 多样本汇总表（核心）
```

### 关键文件说明 | Key files

- `{sample}.nlr_annotator.tsv`：每个样本预测出的 NLR 基因表，表头为 `gene_id / nlr_id / type / start / end / strand / motifs`
- `nlr_annotator_summary.tsv`：目录模式的汇总表，比单样本表多一列 `sample`，一行对应一个样本的一条 NLR 记录

## 结果解读 | Interpreting Results

**通俗理解|In plain words:** 结果表里一行 = 预测出的一个 NLR 基因。看 `type` 判断类型，看 `motifs` 了解它带哪些特征 motif。

- `gene_id`：输入序列里的基因/序列 ID
- `nlr_id`：预测出的 NLR 编号
- `type`：NLR 类型（如 CNL、TNL、RNL 等）
- `start / end / strand`：该 NLR 在输入序列上的坐标区间与方向
- `motifs`：命中的 motif 列表（已去重、排序、去掉 `motif_` 前缀），motif 越多通常结构越完整

## 参数选择建议 | Parameter Guidance

- **单样本预测**：`-i 文件 -o 目录` 即可，默认参数
- **多基因组批量**：`-i 目录 -o 目录`，程序自动按样本分目录并出汇总表
- **批量中途被杀、只想补汇总**：`--merge-only -i 结果目录 -o 目录`，跳过 Java 预测直接合并
- **Java 用 conda 环境**：`--java-path ~/miniforge3/envs/xxx/bin/java`
- **距离参数**：除非你对 NLR-Annotator 的判定逻辑很熟悉，否则保持默认

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 输入DNA/CDS FASTA文件或目录｜Input DNA/CDS FASTA file or directory |
| `-o, --output-dir` | `./output` |  | 输出目录｜Output directory |
| `-t, --threads` | `12` | int | 线程数｜Number of threads |
| `--sample-suffix` | `*.fa` |  | 目录模式下文件匹配后缀｜File match suffix for directory mode |
| `--merge-only` | — |  | 只合并已有结果TSV(*.nlr_annotator.tsv),不运行NLR-Annotator｜Merge existing result TSVs only, skip NLR-Annotator |
| `--output-gff` | — |  | 输出GFF文件｜Output GFF file |
| `--output-bed` | — |  | 输出BED文件｜Output BED file |
| `--output-motifs` | — |  | 输出motifs BED文件｜Output motifs BED file |
| `--output-alignment` | — |  | 输出motif比对FASTA｜Output motif alignment FASTA |
| `--jar-path` | `` |  | NLR-Annotator JAR文件路径｜NLR-Annotator JAR file path |
| `--mot-file` | `` |  | mot.txt配置文件路径｜mot.txt config file path |
| `--store-file` | `` |  | store.txt配置文件路径｜store.txt config file path |
| `--java-path` | `java` |  | Java解释器路径(conda env用~/miniforge3/envs/xxx/bin/java)｜Java interpreter path |
| `--num-seqs-per-thread` | `1000` | int | 每线程处理序列数｜Sequences per thread |
| `--distance-within-motif-combination` | `500` | int | motif组合内距离｜Distance within motif combination |
| `--distance-for-elongating` | `2500` | int | 延伸距离｜Distance for elongating |
| `--distance-between-motif-combinations` | `50000` | int | motif组合间距离｜Distance between motif combinations |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 输入DNA/CDS FASTA文件或目录｜Input DNA/CDS FASTA file or directory |
| `-o, --output-dir` | `./output` |  | 输出目录｜Output directory (default: ./output) |
| `--sample-suffix` | `*.fa` |  | 目录模式下文件匹配后缀｜File match suffix for directory mode (default: *.fa) |
| `--merge-only` | — | store_true | 只合并已有结果TSV(*.nlr_annotator.tsv),不运行NLR-Annotator｜Merge existing result TSVs only, skip NLR-Annotator |
| `--jar-path` | `` |  | NLR-Annotator JAR文件路径｜NLR-Annotator JAR file path |
| `--mot-file` | `` |  | mot.txt配置文件路径｜mot.txt config file path |
| `--store-file` | `` |  | store.txt配置文件路径｜store.txt config file path |
| `--java-path` | `java` |  | Java解释器路径(默认系统java;conda env用~/miniforge3/envs/xxx/bin/java)｜Java interpreter path |
| `-t, --threads` | `12` | int | 线程数｜Number of threads (default: 12) |
| `--num-seqs-per-thread` | `1000` | int | 每线程处理序列数｜Sequences per thread (default: 1000) |
| `--output-gff` | — | store_true | 输出GFF文件｜Output GFF file |
| `--output-bed` | — | store_true | 输出BED文件｜Output BED file |
| `--output-motifs` | — | store_true | 输出motifs BED文件｜Output motifs BED file |
| `--output-alignment` | — | store_true | 输出motif比对FASTA｜Output motif alignment FASTA |
| `--distance-within-motif-combination` | `500` | int | motif组合内距离｜Distance within motif combination (default: 500) |
| `--distance-for-elongating` | `2500` | int | 延伸距离｜Distance for elongating (default: 2500) |
| `--distance-between-motif-combinations` | `50000` | int | motif组合间距离｜Distance between motif combinations (default: 50000) |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- Java（默认系统 `java`，可 `--java-path` 指定）
- NLR-Annotator JAR（默认 `~/software/NLR-Annotator/NLR-Annotator-v2.1b.jar`，可 `NLR_ANNOTATOR_PATH` 覆盖）
- `mot.txt`（默认 `~/software/NLR-Annotator/mot.txt`）
- `store.txt`（默认 `~/software/NLR-Annotator/store.txt`）

## 常见问题 | FAQ

**Q1：报错「JAR 文件不存在 / mot.txt 不存在 / store.txt 不存在」？**
NLR-Annotator 需要这三件套齐全。检查默认路径 `~/software/NLR-Annotator/` 下是否有这三个文件，装在不同位置就用 `--jar-path` / `--mot-file` / `--store-file` 指定。

**Q2：Java 报错或找不到命令？**
用 conda 环境的 Java 时，需显式 `--java-path ~/miniforge3/envs/xxx/bin/java`（默认的裸 `java` 可能指向系统 Java，版本或环境不对）。

**Q3：批量跑到一半被杀了，怎么接着跑？**
重跑同一命令即可，已完成样本会因断点续传自动跳过；若只是想补一张汇总表，用 `--merge-only` 直接合并已有 `*.nlr_annotator.tsv`。

**Q4：结果表里的 motifs 列为什么和原工具不一样？**
程序做了清洗：加了表头、把 motif 去重并排序、去掉 `motif_` 前缀，让结果更规整、便于下游处理。

**Q5：目录模式怎么确认匹配到了哪些文件？**
程序会按 `--sample-suffix`（默认 `*.fa`）在目录里 glob 匹配并打印「发现 N 个文件」。若你的文件后缀不是 `.fa`，用 `--sample-suffix` 调整，例如 `*.cds.fa`。
