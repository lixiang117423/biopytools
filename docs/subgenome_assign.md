# 亚基因组归属 | Subgenome Assignment

一句话理解：**把目标多倍体基因组的每一条染色体，逐个比对到各个亲本参考上，看它跟哪个亲本更像，就把它判给哪个亚基因组**，解决「多倍体里哪条染色体来自哪个祖先」的问题。

## 功能概述 | Overview

- 通过 minimap2 把目标多倍体基因组比对到各亲本参考，按每条染色体的匹配碱基数判定亚基因组来源
- 完整流程：合并亲本 → 比对 → 归属判定 → 拆分 FASTA
- 支持任意多个亲本（至少 2 个），每个亲本可提供多个 hap（单倍型）FASTA
- 输出归属表（TSV）+ 按亚基因组拆分的 FASTA，低置信度染色体显式标记
- 断点续传：已生成的 PAF 非空则跳过比对

## 快速开始 | Quick Start

```bash
biopytools subgenome-assign -i Cf.chr.fa --parent Ca:Ca_hap1.fa,Ca_hap2.fa --parent Ch:Ch_hap1.fa,Ch_hap2.fa -o results -t 16
```

最小输入：一个目标多倍体 FASTA、至少两个亲本（每个亲本一个或多个 hap FASTA）。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| 多倍体 | 基因组里有来自多个祖先的多套染色体（如异源四倍体有 A、B 两套） |
| 亚基因组 | 多倍体里来自某个亲本/祖先的那套染色体 |
| 亲本（parent） | 提供参考序列的祖先物种；归属判定就是「跟谁像就归谁」 |
| hap（单倍型） | 一个亲本的两套等位序列（hap1/hap2），合并后作为该亲本的参考 |
| 比对（alignment） | 把目标序列贴到亲本参考上，看哪里能对上 |
| PAF | 比对结果的文本格式，记录每条比对匹配了多少碱基 |
| 置信度（confidence） | 「这条染色体更像谁」的把握程度，最高分亲本占总分的比例 |
| UNASSIGNED | 没比上任何亲本、无法归属的染色体 |

## 输入 | Input

### 目标多倍体基因组

`-i` / `--target`，FASTA 格式。

### 亲本配置

`--parent NAME:hap1.fa,hap2.fa`，可重复指定多个亲本。格式：亲本名 + 冒号 + 逗号分隔的 hap FASTA 路径（一个亲本可给多个 hap）。

```text
biopytools subgenome-assign -i Cf.chr.fa     --parent Ca:Ca_hap1.fa,Ca_hap2.fa     --parent Ch:Ch_hap1.fa,Ch_hap2.fa     -o results -t 16
```

## 参数说明 | Parameters

### 必需参数 | Required

**通俗理解|In plain words:** 目标基因组（-i）和亲本配置（--parent）缺一不可。亲本至少要给 2 个，否则无法判定归属；每个亲本的名字随意起（会出现在输出文件名和归属表里）。

### 比对参数 | Alignment

**通俗理解|In plain words:** `--preset` 是 minimap2 的 -x 预设，对应目标与亲本的序列差异大小：asm5 适合差异 <1%、asm10 适合 1-5%、asm20 适合 5-15%。亲缘关系近用 asm5，远用 asm20，默认 asm10 覆盖大多数情况。`--minimap2-secondary` 默认关闭次要比对（只留最佳比对），一般不用开。

### 归属判定 | Assignment

**通俗理解|In plain words:** `--min-conf` 是置信度门槛：某染色体「最高分亲本占总分的比例」低于它，就标记为 LOW_CONFIDENCE 提醒你人工复核。默认 0.65，值越大越严格（更多染色体被判低置信），值越小越宽松。

### FASTA 拆分 | FASTA splitting

**通俗理解|In plain words:** 默认会把目标基因组按归属拆成「每个亲本一个 FASTA」并保留未归属的 FASTA。`--no-split` 不拆分；`--no-keep-unassigned` 不保留未归属染色体。

### 工具路径 | Tool paths

**通俗理解|In plain words:** minimap2、samtools 默认取环境里的路径，安装位置不同时再指定。

## 分析流程 | Pipeline

```text
输入目标FASTA + 亲本hap列表
  -> Step 1: 检查依赖(minimap2/samtools)
  -> Step 2: 合并每个亲本的hap -> 01_alignment/<parent>.combined.fa
  -> Step 3: minimap2比对 -> 01_alignment/target_vs_<parent>.paf
  -> Step 4: 解析PAF、逐染色体判定归属 -> 02_assignment/subgenome_assignment.tsv
  -> Step 5: 按归属拆分FASTA -> 03_split_fastas/subgenome_<parent>.fa
  -> 写出00_pipeline_info元数据
```

## 输出 | Output

```text
subgenome_assign_output/
├── 00_pipeline_info/
│   ├── software_versions.yml        # 软件版本信息
│   └── pipeline_params.yaml         # 参数与归属结果汇总
├── 01_alignment/
│   ├── <parent>.combined.fa         # 每个亲本合并后的参考
│   └── target_vs_<parent>.paf       # 目标 vs 每个亲本的比对
├── 02_assignment/
│   └── subgenome_assignment.tsv     # 归属表(核心结果)
├── 03_split_fastas/
│   ├── subgenome_<parent>.fa        # 每个亚基因组的染色体FASTA
│   └── subgenome_UNASSIGNED.fa      # 未归属染色体(默认保留)
└── 99_logs/
    └── subgenome_assign.log
```

关键文件说明：

- **subgenome_assignment.tsv**：每条染色体一行，含各亲本匹配碱基数（`<parent>_match`）、最高分、总分、归属（assigned）、置信度（confidence）、状态（OK/LOW_CONFIDENCE）。
- **subgenome_<parent>.fa**：按归属拆好的各亚基因组 FASTA，可直接用于下游。

## 结果解读 | Interpreting Results

**通俗理解|In plain words:** 看归属表，重点是「归属是否明确、有多少低置信/未归属」。

- **assigned 列**：每条染色体归属到哪个亲本（对应亚基因组）。
- **confidence 列**：最高分亲本占总匹配分的比例；接近 1 说明归属明确，接近 0.5 说明两个亲本势均力敌、归属不可靠。
- **status 列**：LOW_CONFIDENCE 表示置信度低于 --min-conf，需人工复核（日志里会列出这些染色体及分数）。
- **UNASSIGNED**：完全没比上任何亲本的染色体，可能是物种特有序列或比对失败。
- 归属结果质量好坏的直观判据：LOW_CONFIDENCE 和 UNASSIGNED 的染色体越少越好；若大量染色体低置信，说明亲本参考与目标亲缘关系过远或选错了亲本。

## 参数选择建议 | Parameter Guidance

- **--preset**：亲缘近（同种内）用 asm5；近缘种用 asm10（默认）；较远用 asm20/asm25。选太严会比对不上（更多 UNASSIGNED），选太松会误比到错误位置。
- **亲本选择**：务必用真实的祖先/亲本物种参考；亲本选错，归属结果会整体失真。
- **--min-conf**：结果里低置信多就调低（放宽）；想更保守就调高。
- **--no-split**：只需要归属表、不需要拆 FASTA 时用，节省磁盘。

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --target` | 必填 |  | 目标多倍体基因组 FASTA｜Target polyploid genome FASTA |
| `--parent` | 必填 |  | 亲本配置（可重复指定多个亲本）｜Parent spec (can be repeated for multiple parents). 格式｜Format: NAME:hap1.fa,hap2.fa |
| `-o, --output-dir` | `./subgenome_assign_output` |  | 输出目录｜Output directory |
| `--preset` | `asm10` | asm5/asm10/asm20/asm25 | minimap2 -x 预设｜minimap2 preset (asm5=<1%%, asm10=1-5%%, asm20=5-15%%) |
| `-t, --threads` | `12` | int | 线程数｜Number of threads |
| `--minimap2-secondary` | — |  | 保留次要比对（默认 --secondary=no）｜Keep secondary alignments |
| `--min-conf` | `0.65` | float | 置信度阈值｜Confidence threshold for LOW_CONFIDENCE flag |
| `--no-split` | — |  | 不输出拆分的 FASTA｜Do not output split FASTAs |
| `--no-keep-unassigned` | — |  | 不输出未归属染色体的 FASTA｜Do not output unassigned FASTA |
| `--minimap2-path` | `~/miniforge3/envs/align/bin/minimap2` |  | minimap2 二进制路径｜minimap2 binary path |
| `--samtools-path` | `~/.local/bin/samtools` |  | samtools 二进制路径｜samtools binary path |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --target` | 必填 |  | 目标多倍体基因组 FASTA｜Target polyploid genome FASTA |
| `--parent` | 必填 | append | 亲本名:hap1,hap2,...（可重复指定多个亲本）｜Parent name:hap1,hap2,... (can be repeated for multiple parents) |
| `-o, --output-dir` | `./subgenome_assign_output` |  | 输出目录｜Output directory |
| `--preset` | `asm10` | asm5/asm10/asm20/asm25 | minimap2 -x 预设（asm5=<1%%差异, asm10=1-5%%, asm20=5-15%%）｜minimap2 preset |
| `-t, --threads` | `12` | int | 线程数｜Number of threads |
| `--minimap2-secondary` | — | store_true | 保留次要比对（默认 --secondary=no）｜Keep secondary alignments |
| `--min-conf` | `0.65` | float | 置信度阈值，低于此值标记 LOW_CONFIDENCE｜Confidence threshold for LOW_CONFIDENCE flag |
| `--no-split` | — | store_true | 不输出拆分的 FASTA｜Do not output split FASTAs |
| `--no-keep-unassigned` | — | store_true | 不输出未归属染色体的 FASTA｜Do not output unassigned FASTA |
| `--minimap2-path` | — |  | minimap2 二进制路径(默认域环境自动解析)｜minimap2 binary path (default: auto domain env) |
| `--samtools-path` | — |  | samtools 二进制路径(默认域环境自动解析)｜samtools binary path (default: auto domain env) |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- minimap2（自动解析 align 域环境并经 conda run 调用；可用 --minimap2-path 或环境变量 MINIMAP2_PATH 覆盖；域环境缺失时回退 PATH 直接调用）
- samtools（自动解析 align 域环境并经 conda run 调用；可用 --samtools-path 或环境变量 SAMTOOLS_PATH 覆盖；域环境缺失时回退 PATH 直接调用）

## 常见问题 | FAQ

**Q1：断点续传怎么生效？**
比对阶段会检查 PAF 是否存在且非空，存在则跳过比对。但归属判定和拆分每次都会重跑（很快）。想完全重跑请删除 01_alignment 里的 PAF。

**Q2：报「亲本数量必须 >= 2」？**
归属判定需要至少两个亲本对比。确认 --parent 至少指定了 2 个。

**Q3：--parent 格式报错？**
必须写成 `名称:hap1.fa,hap2.fa` 形式（冒号分隔名称与文件、逗号分隔多个 hap）。名称不能为空，至少给一个 hap 文件。

**Q4：大量染色体 UNASSIGNED 或 LOW_CONFIDENCE 怎么办？**
多半是亲本参考与目标亲缘关系过远，或 --preset 选得太严。先放宽 preset（如 asm20），再检查亲本物种是否选对。

**Q5：BAM 还是 FASTA？**
本工具输入是 FASTA（目标基因组 + 亲本 hap FASTA），不接收 BAM。比对由内部 minimap2 完成。

