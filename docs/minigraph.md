# Minigraph 泛基因组图构建与分析 | Minigraph Pangenome Graph Construction and Analysis

一句话理解：**用一张「路网图」表示多个基因组的共同结构，然后在这张图上做四件事——建图、找结构变异（SV）、提取变异气泡、把序列比对回图**，是轻量级泛基因组分析的通用工具。

## 功能概述 | Overview

- 封装 minigraph 与 gfatools，提供四个子命令：build（建图）、call（SV 调用）、bubble（SV 气泡提取）、map（序列到图映射）
- 建图以一条参考为骨架、逐个样本叠加差异，输出标准 GFA 图文件
- 自动检测 minigraph/gfatools 的 conda 环境并包装调用（无需手动指定环境名）
- 支持追加模式（--append-mode）在图已存在时继续叠加新样本
- 日志写入 99_logs/minigraph_pipeline.log，命令执行完整记录

## 快速开始 | Quick Start

```bash
biopytools minigraph build --ref ref.fa --samples s1.fa s2.fa -o graph.gfa
```

最小输入：一个参考基因组 FASTA + 一个或多个样本基因组 FASTA（--samples 可给多个），输出一张 GFA 图。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| 泛基因组图 | 把多个基因组「叠成一张路网」：共有路段只画一次，差异处画成分叉 |
| GFA | 描述这张路网的标准文本：S 是路段、L 是连接、P 是每条基因组的行走路径 |
| 参考基因组 | 图的「主干道」，先以它为骨架建图，再往上叠加样本差异 |
| 结构变异（SV） | 大片段层面的差异：插入、缺失、倒位、重复，比单个 SNP 大得多 |
| 气泡（bubble） | 图上的「分叉-再汇合」结构，通常对应一个结构变异位点 |
| GAF | 序列比对到图上的结果格式（图版本的 SAM/PAF） |
| 预设（preset） | 一组为特定场景调好的默认参数，建图/映射用不同的预设 |
| BED | 用染色体+起止坐标描述区间的一种列格式，SV 结果常用它表示 |

## 输入 | Input

各子命令的输入都是标准 FASTA / GFA：

- build：--ref 参考 FASTA + --samples 样本 FASTA（可多个）
- call：--graph-gfa 已建好的图 + --samples 样本 FASTA（可多个）
- bubble：--graph-gfa 已建好的图
- map：--graph-gfa 图 + --queries 查询 FASTA（可多个）

FASTA 要求序列名唯一、无非法字符；程序会校验文件存在与非空。

## 参数说明 | Parameters

### build 子命令 | Build graph

**通俗理解|In plain words:** 这是最常用的命令。--preset 选「图建多细」：g 只抓大结构、gs 中等、ggs 最细（默认，抓得最全但也更慢更大）。--min-identity / --min-aln-len / --max-gap 是三个相似度/长度/断口阈值，**一般用默认即可**，图太碎或太大时才调。

- --ref：参考基因组 FASTA（必填）
- --samples：样本基因组 FASTA 列表（必填，可多个）
- -o, --output-gfa：输出 GFA 路径，默认 ./pangenome.gfa
- --preset：g / gs / ggs，默认 ggs
- --min-identity：最小序列相似度，默认 0.9
- --min-aln-len：最小比对长度，默认 100000
- --max-gap：最大 gap 大小，默认 1000000
- -t, --threads：线程数，默认 16
- --batch-size：批处理大小（MB，minigraph 的 -K）
- --keep-intermediate：保留中间文件
- --append-mode：追加模式（在已有 GFA 上继续叠加）

### call 子命令 | Call structural variants

**通俗理解|In plain words:** 对每个样本，把它比对到图上、找出它与图不一致的地方（结构变异），每个样本出一个 BED 文件。--preset 目前只有 asm（针对组装到图的比对），不用改。

- --graph-gfa：图 GFA 文件（必填）
- --samples：样本 FASTA 列表（必填）
- -o, --output-dir：输出目录，默认 ./minigraph_call（每个样本出 {样本名}.bed）
- --preset：asm（固定）
- -t, --threads：线程数，默认 16

### bubble 子命令 | Extract SV bubbles

**通俗理解|In plain words:** 从图里把「分叉-再汇合」的气泡结构全部挖出来，得到一张 SV 位点清单（BED）。

- --graph-gfa：图 GFA 文件（必填）
- -o, --output-bed：输出 BED 路径，默认 ./sv_bubbles.bed

### map 子命令 | Map sequences to graph

**通俗理解|In plain words:** 把新的序列（reads 或组装）比对到已有的图上，得到 GAF 结果。--preset 按数据类型选：sr 短读、lr 长读（默认）、map-pb PacBio、map-ont 纳米孔、asm 组装。

- --graph-gfa：图 GFA 文件（必填）
- --queries：查询 FASTA 列表（必填）
- -o, --output-gaf：输出 GAF 路径，默认 ./mapping.gaf
- --preset：sr / lr / map-pb / map-ont / asm，默认 lr
- --max-intron-len：最大内含子长度（asm 预设时用）
- -t, --threads：线程数，默认 16
- --batch-size：批处理大小（MB）

### 工具路径 | Tool paths

**通俗理解|In plain words:** 程序会自动探测 minigraph/gfatools 所在的 conda 环境并包装调用。只有你的工具不在 PATH、装在自己知道的地方时，才用这两条参数指路径。

- --minigraph-path：minigraph 工具路径，默认 minigraph
- --gfatools-path：gfatools 工具路径，默认 gfatools

## 分析流程 | Pipeline

```text
build:  ref.fa + samples  ->  minigraph -ggs ...  ->  pangenome.gfa
    |
    v
call:   graph.gfa + samples  ->  minigraph -xasm --call ...  ->  {sample}.bed
    |
    v
bubble: graph.gfa  ->  gfatools bubble ...  ->  sv_bubbles.bed
    |
    v
map:    graph.gfa + queries  ->  minigraph -xlr ...  ->  mapping.gaf
```

四个子命令相互独立，可按需单跑；典型用法是先 build 建图，再 bubble 提 SV、map 做比对。

## 输出 | Output

```text
# build
pangenome.gfa                  # 泛基因组图（S/L/P 记录）
99_logs/minigraph_pipeline.log # 运行日志

# call
minigraph_call/
├── s1.bed                     # 每个样本的结构变异 BED
└── s2.bed

# bubble
sv_bubbles.bed                 # 图中所有 SV 气泡

# map
mapping.gaf                    # 序列到图比对结果
```

## 结果解读 | Interpreting Results

- **pangenome.gfa**：核心结果。用 grep 统计 S/L/P 开头的行数（分别对应片段、连接、路径）；P 行数量应等于输入基因组数。图能在 Bandage/vg 里正常展开即结构合理。
- **{样本名}.bed**：该样本相对图的结构变异清单，行数越多 = 该样本与参考差异越大。可用于下游 SV 注释。
- **sv_bubbles.bed**：图里所有 SV 位点，配合 build 的图一起看，气泡数代表样本间的结构变异规模。
- **mapping.gaf**：比对上图的结果，每行一条比对，可用于图上的变异检测或覆盖分析。
- **好坏判据**：命令返回 0、输出文件生成且非空、日志无 ERROR 即成功。GFA 中 P 行数不足 = 有样本没成功并进图。

## 参数选择建议 | Parameter Guidance

- **--preset（build）**：默认 ggs 最全；样本间差异小、只想看大结构用 g/gs（更快更小）；追求完整 SV 谱用 ggs。
- **--min-identity / --min-aln-len**：图里碎片化气泡太多可调高（更保守）；想抓远缘差异可调低。一般不动。
- **--append-mode**：分批加样本时用，避免每次从头重建；注意追加的图要与原图参考一致。
- **--keep-intermediate**：排查建图细节时保留中间产物，正常跑不必加。
- **--batch-size（-K）**：内存紧张时调小（如 500M），速度略有取舍，一般不用设。

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--ref` | 必填 |  | 参考基因组FASTA文件｜Reference genome FASTA file |
| `--samples` | 必填 |  | 样本基因组FASTA文件｜Sample genome FASTA files |
| `-o, --output-gfa` | `./pangenome.gfa` |  | 输出GFA文件路径｜Output GFA file path |
| `--preset` | `ggs` | g/gs/ggs | 图构建预设｜Graph building preset |
| `--min-identity` | `0.9` | float | 最小序列相似度｜Minimum sequence identity |
| `--min-aln-len` | `100000` | int | 最小比对长度｜Minimum alignment length |
| `--max-gap` | `1000000` | int | 最大gap大小｜Maximum gap size |
| `-t, --threads` | `16` | int | 线程数｜Number of threads |
| `--batch-size` | — | int | 批处理大小(MB)｜Batch size (MB) |
| `--minigraph-path` | `minigraph` |  | minigraph工具路径｜minigraph tool path |
| `--gfatools-path` | `gfatools` |  | gfatools工具路径｜gfatools tool path |
| `--keep-intermediate` | — |  | 保留中间文件｜Keep intermediate files |
| `--append-mode` | — |  | 追加模式｜Append mode |
| `--graph-gfa` | 必填 |  | 泛基因组图GFA文件｜Pangenome graph GFA file |
| `-o, --output-dir` | `./minigraph_call` |  | 输出目录｜Output directory |
| `-o, --output-bed` | `./sv_bubbles.bed` |  | 输出BED文件路径｜Output BED file path |
| `--queries` | 必填 |  | 查询序列FASTA文件｜Query sequence FASTA files |
| `-o, --output-gaf` | `./mapping.gaf` |  | 输出GAF文件路径｜Output GAF file path |
| `--max-intron-len` | — | int | 最大内含子长度｜Maximum intron length |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--ref` | 必填 |  | 参考基因组FASTA文件｜Reference genome FASTA file |
| `--samples` | 必填 |  | 样本基因组FASTA文件列表｜Sample genome FASTA file list |
| `-o, --output-gfa` | `./pangenome.gfa` |  | 输出GFA文件路径｜Output GFA file path |
| `--preset` | `ggs` | g/gs/ggs | 图构建预设｜Graph building preset |
| `--min-identity` | `0.9` | float | 最小序列相似度｜Minimum sequence identity |
| `--min-aln-len` | `100000` | int | 最小比对长度｜Minimum alignment length |
| `--max-gap` | `1000000` | int | 最大gap大小｜Maximum gap size |
| `-t, --threads` | `16` | int | 线程数｜Number of threads |
| `--batch-size` | — | int | 批处理大小(MB)｜Batch size (MB) |
| `--minigraph-path` | `minigraph` |  | minigraph工具路径｜minigraph tool path |
| `--gfatools-path` | `gfatools` |  | gfatools工具路径｜gfatools tool path |
| `--keep-intermediate` | — | store_true | 保留中间文件｜Keep intermediate files |
| `--append-mode` | — | store_true | 追加模式｜Append mode |
| `--graph-gfa` | 必填 |  | 泛基因组图GFA文件｜Pangenome graph GFA file |
| `-o, --output-dir` | `./minigraph_call` |  | 输出目录｜Output directory |
| `-o, --output-bed` | `./sv_bubbles.bed` |  | 输出BED文件路径｜Output BED file path |
| `--queries` | 必填 |  | 查询序列FASTA文件列表｜Query sequence FASTA file list |
| `-o, --output-gaf` | `./mapping.gaf` |  | 输出GAF文件路径｜Output GAF file path |
| `--max-intron-len` | — | int | 最大内含子长度｜Maximum intron length |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- **minigraph**（建图/SV 调用/映射核心，自动探测 conda 环境；若无则在 PATH 中查找）
- **gfatools**（bubble 提取用；可选，不影响其他子命令）
- Python 3（封装脚本自身）

## 常见问题 | FAQ

**Q1：有没有断点续传？**
本工具不设断点续传，每次都会重新执行对应子命令。若不想重算整张图，用 build 的 --append-mode 在已有 GFA 上追加新样本；否则换输出文件路径保存。

**Q2：minigraph 报「not available」？**
程序用 which 探测 minigraph 并尝试 conda 包装。确认 which minigraph 能找到它，或已装在某 conda 环境且该环境被当前 conda 发现；找不到时用 --minigraph-path 指到实际路径。

**Q3：GFA 里 P 行比样本数少？**
说明有样本没成功并入图（可能序列太短/太相似被合并，或参数过严）。检查日志与 min-identity/min-aln-len 设置。

**Q4：call 出的 BED 文件在哪、叫什么？**
在 --output-dir 下，每个样本一个，文件名 = 样本 FASTA 去掉扩展名 + .bed。

**Q5：map 的 preset 怎么选？**
短读选 sr，PacBio 长读选 map-pb，纳米孔选 map-ont，一般长读选 lr（默认），组装到图选 asm（可配合 --max-intron-len）。
