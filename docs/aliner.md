# a-liner 共线性可视化 | a-liner Synteny Visualization (FASTA→minimap2→plot)

一句话理解：**把两条基因组的染色体"并排"画成连线图，用彩色连线标出"同一段基因在谁家哪条染色体上"的对应关系（共线性）**，一眼看出两个物种之间的基因组重排。

## 功能概述 | Overview

- 输入参考（ref）+ 查询（query）两条基因组 FASTA，自动跑 minimap2 比对、再用 a-liner 画共线性图
- 支持指定染色体或区段（如 `chrZ` 或 `chrZ:1-30000000`），ref 与 query 按配对顺序交错排列
- 内置校验：序列名是否存在、区段是否越界、ref/query 序列数量必须一致
- 断点续传：minimap2 的 PAF 文件已存在且非空时自动跳过
- 输出 PDF 图 + 版本元数据（software_versions.yml + pipeline_params.yaml）

## 快速开始 | Quick Start

```bash
biopytools aliner --ref ref.fa --query query.fa --ref-seqs chr1,chr2 --query-seqs Chr1,Chr2 -o output
```

最小输入：两条基因组 FASTA + 两侧要画的序列名（逗号分隔，两侧数量必须一致）。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| 共线性(synteny) | 两个物种基因组里"基因排列顺序相似"的一段区域，像两本同源书按相同章节顺序排 |
| minimap2 | 快速 DNA 序列比对软件，像"查重引擎"，找出两段序列的相似片段 |
| PAF | 比对结果文件格式（Pairwise mApping Format），逐行记录每条比对片段的坐标与质量 |
| identity | 比对一致性；100=完全相同，70≈七成碱基一致 |
| preset | minimap2 预设参数组；asm5=近缘（>95% 一致），asm20=远缘 |
| 参考/查询(ref/query) | 两条待比较的基因组；ref 画上方，query 画下方 |

## 输入 | Input

- `--ref`：参考基因组 FASTA（会被 samtools faidx 建索引取长度）
- `--query`：查询基因组 FASTA
- `--ref-seqs` / `--query-seqs`：两侧要展示的序列，逗号分隔；支持整条 `chrZ` 或区段 `chrZ:1-30000000`（1-based 闭区间），**两侧数量必须相等**（配对交错排列）

## 参数说明 | Parameters

### 输入参数 | Input

**通俗理解|In plain words:** 四个参数里 `--ref` 和 `--query` 是两条待比较的基因组文件，`--ref-seqs` 和 `--query-seqs` 是"这次只画哪些染色体/区段"。ref 和 query 两侧按位置一一配对，所以**两侧列表长度必须相同**，否则直接报错。

### 输出参数 | Output

**通俗理解|In plain words:** `-o` 指定结果目录，`--out-prefix` 是最终 PDF 和中间 PAF 的文件名前缀。**一般用默认前缀 `synteny` 即可**；同一目录跑多组比较时才改前缀避免互相覆盖。

### 比对参数 | Alignment

**通俗理解|In plain words:** `--preset` 按两个基因组的亲缘关系选：越近选 asm5（更敏感、更严），越远选 asm20（更宽松）。`--min-identity` 和 `--min-alignment-len` 是 a-liner 画图时对 PAF 的过滤线——低于 identity 或短于长度的比对片段不画。**近缘物种用默认值即可**；画出来的连线太稀疏时再调低这两个阈值。

### 可视化参数 | Visualization

**通俗理解|In plain words:** `--colormap`（0-5）只是换配色，`--figure-size` 设图的宽高（高度填 0 表示按内容自适应）。**一般不用动**，图太挤时再放大宽度。

### 高级参数 | Advanced

**通俗理解|In plain words:** `--extra-args` 是把字符串原样透传给底层 a-liner 的兜底开关，用于传本文档未覆盖的 a-liner 高级参数。**普通用户用不到**。

## 分析流程 | Pipeline

```text
参考FASTA + 查询FASTA + 序列规格
    │
    ▼
samtools faidx 取各序列长度（无 .fai 先建索引）
    │
    ▼
校验：序列名存在？区段越界？ref/query 数量一致？
    │
    ▼
minimap2 比对 → PAF（断点续传：PAF 已存在则跳过）
    │
    ▼
生成 sequence_config.txt（ref/query 配对交错排布）
    │
    ▼
a-liner 出图 → PDF
    │
    ▼
写版本元数据（software_versions.yml / pipeline_params.yaml）
```

## 输出 | Output

```text
output/
├── 00_pipeline_info/
│   ├── software_versions.yml       # 各软件版本与路径
│   └── pipeline_params.yaml        # 本次运行参数快照
├── 01_alignment/
│   ├── synteny.paf                 # minimap2 比对结果（可复用，断点续传依据）
│   └── sequence_config.txt         # 传给 a-liner 的序列排布表（n/ID/start/end/strand/name）
├── 02_aliner/
│   └── synteny.pdf                 # 最终共线性图
└── 99_logs/
    └── aliner_pipeline.log         # 运行日志
```

## 结果解读 | Interpreting Results

### 1. 共线性图（`synteny.pdf`）

**通俗理解|In plain words:** 这是核心结果。上方一排是 ref 序列，下方一排是 query 序列，中间的彩色连线表示"ref 这一段和 query 那一段是同一个东西"。

- 连线整体平行、成束：两个基因组这段基本同序（共线性好）
- 连线交叉打结：发生了倒位（方向反了）或易位（跑到别的染色体）
- 连线稀疏/断裂：该区域分化较大，或被 identity/长度阈值过滤掉了

### 2. 比对文件（`synteny.paf`）

**通俗理解|In plain words:** 中间产物，PAF 第 12 列以后有 block 数等指标。普通用户无需细看；改比对参数重跑前删掉它，否则会复用旧比对。

### 3. 序列排布表（`sequence_config.txt`）

**通俗理解|In plain words:** a-liner 的"舞台布置表"，决定哪条序列画在上面、哪条画在下面、画全长还是画区段。想核对画的是不是想要的区段，看这个文件的 start/end 列即可。

## 参数选择建议 | Parameter Guidance

- `--preset`：同物种/近缘物种用默认 `asm5`；跨属或更远用 `asm10`/`asm20`
- `--min-identity`：近缘默认 70 即可；想多看远缘保守片段降到 50-60，想只看高置信片段升到 80-90
- `--min-alignment-len`：想看大片共线性保持默认 1000；要排除小片段杂音可调大
- `--figure-size`：高度填 0 自适应即可，图太挤把宽调到 8-12
- `--ref-seqs`/`--query-seqs`：只想看某条染色体就用 `chrZ`，只想看局部用 `chrZ:1-30000000`

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--ref` | 必填 |  | 参考基因组FASTA｜Reference genome FASTA |
| `--query` | 必填 |  | 查询基因组FASTA｜Query genome FASTA |
| `--ref-seqs` | 必填 |  | ref侧序列(逗号分隔,如 chrZ,chrW 或 chrZ:1-30000000)｜ref-side seqs |
| `--query-seqs` | 必填 |  | query侧序列(逗号分隔,与ref等长)｜query-side seqs |
| `-o, --output-dir` | `./aliner_output` |  | 输出目录｜Output directory |
| `--out-prefix` | `synteny` |  | 输出前缀｜Output prefix |
| `--preset` | `asm5` | asm5/asm10/asm20 | minimap2预设(近缘asm5)｜minimap2 preset |
| `--min-identity` | `70` | int | identity阈值%%｜identity threshold |
| `--min-alignment-len` | `1000` | int | 最小比对长度｜min alignment length |
| `--colormap` | `5` | 0/1/2/3/4/5 | 配色｜colormap |
| `--figure-size` | `[6, 0]` | float | 图尺寸[宽 高]｜figure size |
| `-t, --threads` | `12` | int | 线程数｜Threads |
| `--extra-args` | `` |  | 透传a-liner参数｜pass-through args |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--ref` | 必填 |  | 参考基因组FASTA｜Reference genome FASTA |
| `--query` | 必填 |  | 查询基因组FASTA｜Query genome FASTA |
| `--ref-seqs` | 必填 |  | ref侧序列（逗号分隔，如 chrZ,chrW 或 chrZ:1-30000000）｜ref-side seqs (comma-separated) |
| `--query-seqs` | 必填 |  | query侧序列（逗号分隔，与ref等长）｜query-side seqs (comma-separated, equal length to ref) |
| `-o, --output-dir` | `./aliner_output` |  | 输出目录｜Output directory |
| `--out-prefix` | `synteny` |  | 输出文件前缀｜Output prefix |
| `--preset` | `asm5` | asm5/asm10/asm20 | minimap2预设（近缘asm5/远缘asm10,asm20）｜minimap2 preset |
| `--min-identity` | `70` | int | a-liner identity阈值(%%)｜identity threshold |
| `--min-alignment-len` | `1000` | int | 最小比对长度｜min alignment length |
| `--colormap` | `5` | 0/1/2/3/4/5 | 配色(0-5)｜colormap |
| `--figure-size` | `[6, 0]` | float | 图尺寸[宽 高]，高0自适应｜figure size [w h], h=0 auto |
| `--threads, -t` | `12` | int | 线程数｜Threads |
| `--extra-args` | `` |  | 透传给a-liner的额外参数｜pass-through args to a-liner |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- Python 3
- minimap2（conda 环境 `align`，默认 `~/miniforge3/envs/align/bin/minimap2`）
- samtools（conda 环境 `align`，默认 `~/miniforge3/envs/align/bin/samtools`，用于取序列长度）
- a-liner（固定 conda 环境 `a-liner`）

## 常见问题 | FAQ

**Q1：换比对参数重跑，为什么结果没变？**
断点续传按 `01_alignment/synteny.paf` 是否存在判断。改了 `--preset`/`--threads` 等比对参数后，先删掉旧 PAF 再重跑，否则会复用旧比对。

**Q2：报"ref-seqs 与 query-seqs 数量不一致"？**
两侧是配对关系，逗号分隔后数量必须相等。例如 ref 画 2 条、query 也要画 2 条。

**Q3：报"seq not in ref FASTA"？**
序列名大小写要完全一致（FASTA 头里 `>chr1` 和参数里 `Chr1` 不同）。日志会打印可用序列名列表，照抄即可。

**Q4：a-liner 报错或找不到？**
a-liner 走固定 conda 环境 `a-liner`，确认该环境已安装且可 `conda run -n a-liner a-liner --version`。

**Q5：想传 a-liner 的高级参数怎么办？**
用 `--extra-args` 原样透传，例如 `--extra-args "--font_size 8"`（空格分隔的参数列表）。
