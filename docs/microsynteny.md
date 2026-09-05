# Microsynteny 微观共线性分析 | Microsynteny Analysis

一句话理解：**给定几个物种里你关心的「目标基因」，自动找到它们附近的小范围（微观）共线性区域，并画成一张圈图，帮你看清这些基因在演化中是否「排排坐」地保守下来**。

## 功能概述 | Overview

- 基于 **JCVI** 的经典微观共线性流程（GFF→BED→CDS→比对→anchors→blocks），再交给 **pyCirclize** 画圈图
- 四步流水线，每步都有完成标记，支持断点续传（`--step` 可单独跑某一步）
- 输入极简：一个基因组文件夹（每个物种一套 .fa + .gff）+ 一张目标基因清单
- 自动做物种两两组合的共线性分析，最后把含目标基因的染色体段整合进一张多物种圈图
- 目标基因周围可指定「左右各延伸多少个基因」来控制显示范围

## 快速开始 | Quick Start

```bash
biopytools microsynteny -i ./genome_data -g genes.txt
```

最小输入：一个基因组文件夹（含 A.fa、A.gff、B.fa、B.gff 等）+ 一张两列基因清单（物种 ID + 基因 ID）。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| 微观共线性 | 不看整条染色体，只看一小段（一个基因家族附近）里基因的排列是否保守 |
| GFF | 基因注释文件：每条记录标注某段序列上有什么（基因、外显子等）及其位置 |
| BED | 更精简的位置表：染色体 + 起止 + 名字 + 方向，比对程序爱吃这种格式 |
| CDS | 编码序列，即真正翻译成蛋白质的那部分 DNA，比对用的是它 |
| anchors | 两个物种间「一对一匹配」的基因对（锚点），是共线性的骨架 |
| blocks | 由相邻锚点串成的「共线性块」，一段同源且顺序一致的基因串 |
| cscore | 共线性打分阈值：越高越严格，只保留「几乎确定」的锚点对 |
| 延伸基因数 | 以目标基因为中心，左右各再带多少个基因一起画，决定圈图窗口大小 |

## 输入 | Input

### 基因组文件夹（-i/--genome-folder）

每个物种放一对同名文件：序列（.fa 或 .fasta）+ 注释（.gff 或 .gff3）。**文件名前缀（去扩展名）就是物种 ID**，必须与基因清单里的物种 ID 一致。

```text
genome_data/
├── Ccu.fa
├── Ccu.gff
├── Ceq.fa
└── Ceq.gff
```

GFF 中基因的 ID 属性（ID=...）会被用作基因 ID，务必与基因清单里的基因 ID 一致。

### 基因清单（-g/--gene-list）

两列，制表符或空格分隔：**物种 ID + 基因 ID**。

```text
Ccu    gene001
Ccu    gene002
Ceq    gene003
```

- `#` 开头与空行跳过
- 一个物种可列多个基因；至少要有一个物种有目标基因

## 参数说明 | Parameters

### 必需参数 | Required

**通俗理解|In plain words:** 告诉程序「基因组放在哪、关心哪些基因」。两者缺一不可。

- -i, --genome-folder：基因组文件夹（含 {物种}.fa 与 {物种}.gff）
- -g, --gene-list：目标基因清单（两列：物种ID<空格/制表符>基因ID）

### 环境路径 | Environment paths

**通俗理解|In plain words:** 本工具拆成两段跑：前 3 步用 JCVI 环境，第 4 步绘图用 pyCirclize 环境。这两个路径指向对应 conda 环境，**默认位置装好了就一般不用动**。

- -j, --jcvi-path：JCVI 环境路径，默认 ~/miniforge3/envs/jcvi_v.1.5.7
- --pycirclize-path：pyCirclize 环境路径，默认 ~/miniforge3/envs/pycirclize_v.1.10.1

### 分析参数 | Analysis

**通俗理解|In plain words:** --cscore 是「锚点要多可靠才算数」，调大 = 只保留最确定的共线性、图更干净但可能丢信号；--extend-genes 是「以目标基因为中心向两边看多少个基因」，调大 = 窗口更大、图更全，调小 = 更聚焦。**默认值适合多数场景，一般不用动**。

- -t, --threads：线程数，默认 12
- --extend-genes：目标基因左右各延伸的基因数，默认 30
- --cscore：共线性分数阈值（0-1），默认 0.99

### 步骤控制与日志 | Step control & logging

**通俗理解|In plain words:** --step 让你只跑某一步（1 预处理 / 2 共线性 / 3 区块 / 4 绘图），适合断点续传或单独重画图；--log-level 控制日志详细程度。

- --step：1 / 2 / 3 / 4，只跑指定步骤
- --log-level：DEBUG / INFO / WARNING / ERROR / CRITICAL，默认 INFO
- -o, --output-dir：输出目录，默认 ./microsynteny_output

## 分析流程 | Pipeline

```text
基因组文件夹 + 基因清单
    |
    v
步骤1 预处理: GFF -> BED -> 提取 CDS -> 合并所有物种 BED (JCVI)
    |
    v
步骤2 共线性: 物种两两 JCVI ortholog 比对 -> .anchors/.lifted.anchors
    |
    v
步骤3 区块: JCVI mcscan 提取共线性 blocks
    |
    v
步骤4 绘图: pyCirclize 筛选含目标基因的染色体段 -> 多物种圈图
```

每一步结束都会写一个 .step_N.done 标记文件；重跑时已完成的步骤自动跳过。

## 输出 | Output

```text
microsynteny_output/
├── 1_preprocess/
│   ├── Ccu.bed            # 各物种 GFF 转出的 BED
│   ├── Ccu.cds            # 各物种提取的 CDS 序列
│   └── all_species.bed    # 合并后的 BED
├── 2_synteny/
│   ├── Ccu_Ceq.last       # LAST 比对原始结果
│   ├── Ccu_Ceq.anchors    # 过滤后的锚点
│   └── Ccu_Ceq_lifted.anchors  # 抬升后的锚点
├── 3_blocks/
│   ├── Ccu_Ceq.blocks     # 共线性块（基因串）
│   └── blocks_index.txt   # 块文件索引
├── 4_plot/
│   ├── multi_species_synteny.png   # 最终圈图（300 dpi）
│   ├── multi_species_synteny.svg   # 矢量版
│   └── multi_species_synteny.pdf   # PDF 版
└── logs/
    └── microsynteny.log   # 运行日志
```

物种两两组合会各自生成一套 anchors/blocks 文件（文件名 = 物种1.物种2）。

## 结果解读 | Interpreting Results

- **4_plot/multi_species_synteny.png**：最终圈图。每个扇区是一条含目标基因的染色体，目标基因高亮、旁边的基因半透明，橙色连线串起不同物种间的共线性锚点。**连线密集且方向一致 = 该区域微观共线性保守**；连线稀疏/杂乱 = 重排较多。
- **3_blocks/*.blocks**：每行一个共线性块（一串有序的基因 ID），基因 ID 数量越多 = 该块越大。目标基因出现在哪些块里，说明它在哪些物种间有保守邻域。
- **2_synteny/*.anchors**：锚点对数量。锚点太少（接近 0）可能 cscore 太严或物种太远；可下调 --cscore 重跑步骤 2。
- **好坏判据**：四步全部完成（无 ERROR、日志末尾「所有步骤成功完成」），且 4_plot 下 png/svg/pdf 均非空，即成功。

## 参数选择建议 | Parameter Guidance

- **--cscore**：默认 0.99 很严格，适合近缘物种。若锚点/块少得可怜，降到 0.9 或 0.7 试试；近缘且想收紧用默认。
- **--extend-genes**：默认 30 平衡「上下文充分」与「图不太挤」。目标基因密集时可调到 10–20；想看更大邻域背景调到 50–100。
- **--step**：只改绘图或只调 cscore 时，用它跳过无需重跑的步骤（注意下游步骤依赖上游产物）。
- **--jcvi-path / --pycirclize-path**：默认装好不动；只在自己装新版本时覆盖。

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --genome-folder` | 必填 |  | 基因组文件夹路径｜Genome folder path (包含A.fa, A.gff等｜containing A.fa, A.gff, etc.) |
| `-g, --gene-list` | 必填 |  | 目标基因列表文件｜Target gene list file (两列｜two cols: species_id gene_id) |
| `-o, --output-dir` | `./microsynteny_output` |  | 输出目录路径｜Output directory path |
| `-j, --jcvi-path` | `~/miniforge3/envs/jcvi_v.1.5.7` |  | JCVI环境路径｜JCVI environment path |
| `-t, --threads` | `12` | int | 线程数｜Number of threads |
| `--extend-genes` | `30` | int | 延伸基因数｜Number of genes to extend on each side |
| `--cscore` | `0.99` | float | 共线性分数阈值｜Synteny score threshold (0-1) |
| `--step` | — | 1/2/3/4 | 运行指定步骤｜Run specific step only: 1: 数据预处理｜Data preprocessing 2: 共线性分析｜Synteny analysis 3: 区块提取｜Block extraction 4: 绘图｜Plotting |
| `--log-level` | `INFO` | DEBUG/INFO/WARNING/ERROR/CRITICAL | 日志级别｜Log level |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --genome-folder` | 必填 |  | 基因组文件夹路径｜Genome folder path (containing A.fa, A.gff, B.fa, B.gff, ...) |
| `-g, --gene-list` | 必填 |  | 目标基因列表文件｜Target gene list file (two columns: species_id<tab>gene_id) |
| `-o, --output-dir` | `./microsynteny_output` |  | 输出目录路径｜Output directory path (default: ./microsynteny_output) |
| `-j, --jcvi-path` | `~/miniforge3/envs/jcvi_v.1.5.7` |  | JCVI环境路径｜JCVI environment path (default: ~/miniforge3/envs/jcvi_v.1.5.7) |
| `--pycirclize-path` | `~/miniforge3/envs/pycirclize_v.1.10.1` |  | pyCirclize环境路径｜pyCirclize environment path (default: ~/miniforge3/envs/pycirclize_v.1.10.1) |
| `-t, --threads` | `12` | int | 线程数｜Number of threads (default: 12) |
| `--extend-genes` | `30` | int | 延伸基因数｜Number of genes to extend on each side (default: 30) |
| `--cscore` | `0.99` | float | 共线性分数阈值(0-1)｜Synteny score threshold 0-1 (default: 0.99) |
| `--step` | — | 1/2/3/4 | 运行指定步骤｜Run specific step only: 1: 数据预处理｜Data preprocessing 2: 共线性分析｜Synteny analysis 3: 区块提取｜Block extraction 4: 绘图｜Plotting |
| `--log-level` | `INFO` | DEBUG/INFO/WARNING/ERROR/CRITICAL | 日志级别｜Log level (default: INFO) |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- **JCVI**（conda 环境 jcvi_v.1.5.7，负责步骤 1–3 的共线性分析）
- **LAST**（JCVI 的底层比对依赖，需能在 JCVI 环境 bin 或 PATH 中找到 lastdb）
- **pyCirclize**（conda 环境 pycirclize_v.1.10.1，负责步骤 4 绘图）
- Python 3（封装脚本自身，运行在各环境 python 下）

## 常见问题 | FAQ

**Q1：报「缺少基因组文件 Ccu.fa/.fasta 或 Ccu.gff/.gff3」？**
程序按基因清单里的物种 ID 去基因组文件夹找 {ID}.fa/.fasta 与 {ID}.gff/.gff3。物种 ID 必须与文件名前缀完全一致（含大小写）。

**Q2：断点续传怎么用？换了参数为什么没重算？**
每步完成会写 .step_N.done 标记；重跑时已标记的步骤跳过。想用新参数重算某一步，先删掉对应的 .step_N.done（及该步产物），或用 --step 单独指定。

**Q3：报「未找到 LAST 工具」？**
LAST（lastdb 等）是 JCVI 的必需底层比对工具。在 JCVI 环境里 `conda install -c bioconda last`，或保证 PATH 能 which 到 lastdb。

**Q4：共线性几乎没锚点怎么办？**
多半是 --cscore 默认 0.99 太严，或物种亲缘关系太远。先下调 --cscore 到 0.7–0.9，重跑步骤 2；同时确认 GFF 的基因 ID 与清单一致。

**Q5：为什么跳过 dotplot？**
程序在 JCVI ortholog 步骤加了 --no_dotplot，是为了规避某些环境下 zlib 版本导致的绘图失败。dotplot 只是诊断图，不影响 anchors 结果，无需担心。
