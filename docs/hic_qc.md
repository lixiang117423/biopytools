# Hi-C 数据质量评估 | Hi-C Data Quality Assessment

一句话理解：**给 Hi-C 测序数据「体检打分」——读取出比对率、有效配对比例、顺式/反式比例等关键指标，逐项对照阈值给出通过/不通过结论**，解决「这批 Hi-C 数据到底能不能用、差在哪」的问题。

## 功能概述 | Overview

- 支持两种输入模式：HiC-Pro 输出目录（默认）或 pairtools 的 pairs/BAM 文件
- HiC-Pro 模式：解析比对（.mmapstat）、配对（.mpairstat）、有效配对过滤（.mRSstat）三个统计文件
- pairs/BAM 模式：调用 pairtools 计算统计，BAM 输入自动转 pairs（需提供 chrom.sizes）
- 输出逐项 PASS/FAIL 的评估报告，任一项不达标则整体判定未通过（退出码 1）
- 全部阈值可自定义，便于按项目标准收紧或放宽

## 快速开始 | Quick Start

```bash
biopytools hic-qc -i /path/to/hicpro_output
```

最小输入：一个 HiC-Pro 输出目录（默认模式），程序自动检测样本名并评估。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| 比对率 | 读段里有多少能贴回基因组；太低说明数据或基因组有问题 |
| 唯一比对率 | 只贴到唯一位置的读段比例；多贴（重复区）太多会干扰 |
| valid pairs | 通过「两头都正常、方向正确」过滤后的有效配对 |
| dangling ends | 「一头悬空」的异常配对，反映酶切或连接问题 |
| self-ligation | 自己跟自己连上的配对（自连），越低越好 |
| religation | 两个酶切末端重新连上（重连），越低越好 |
| cis / trans | cis=同一染色体内的配对，trans=跨染色体的配对 |
| cis/trans 比例 | 顺式配对与反式配对的比值，越高说明信号越集中在染色体内（Hi-C 正常应偏高） |
| PCR 重复率 | 因扩增产生的重复读段比例，过高浪费数据 |

## 输入 | Input

### HiC-Pro 模式（默认）

一个 HiC-Pro 输出目录，需包含 `bowtie_results/` 和 `hic_results/` 两个子目录。程序自动从 `.mmapstat` / `.mpairstat` / `.mRSstat` 检测样本名（也可用 `-s` 手动指定）。

### pairs / BAM 模式

`--input-type pairs` 用 `.pairs.gz` / `.pairs` 文件；`--input-type bam` 用 BAM/SAM 文件（**必须**同时用 `-c` 提供 chrom.sizes 文件）。

```text
# HiC-Pro 目录
biopytools hic-qc -i /path/to/hicpro_output

# pairs 文件
biopytools hic-qc -i sample.pairs.gz --input-type pairs

# BAM 文件（需 chrom.sizes）
biopytools hic-qc -i sample.bam --input-type bam -c genome.chrom.sizes
```

## 参数说明 | Parameters

### 输入与输出 | Input and output

**通俗理解|In plain words:** `--input-type` 决定走哪套评估逻辑：hicpro（读 HiC-Pro 目录）、pairs（读 .pairs 文件）、bam（读 BAM，需 chrom.sizes）。`-s` 在 hicpro 模式指定样本名，不指定会自动检测。`-o` 是日志和中间统计文件的存放目录。

### HiC-Pro 模式阈值 | HiC-Pro thresholds

**通俗理解|In plain words:** 这组阈值只在 hicpro 模式下用，是各项指标的「及格线」。比对率、唯一比对率、valid pairs 比例是「越高越好」（低于阈值判失败）；dangling ends、self-ligation、religation、PCR 重复率是「越低越好」（高于阈值判失败）。默认值来自 Hi-C 社区常用经验，一般不用改。

### pairtools 模式阈值 | Pairtools thresholds

**通俗理解|In plain words:** 这组阈值只在 pairs/bam 模式下用。未比对比例、单端比对比例、PCR 重复率「越低越好」；双端比对率「越高越好」。`--pairtools-path` 指定 pairtools 可执行文件，`-c/--chroms-path` 是 BAM 输入时的染色体大小文件。

### 通用阈值 | Common thresholds

**通俗理解|In plain words:** `--min-cis-trans-ratio` 是顺式/反式比例的最低线，两种模式共用；正常 Hi-C 该比值明显高于阈值（如 5），若过低说明染色体内信号弱、数据可能有问题。

## 分析流程 | Pipeline

```text
输入(HiC-Pro目录 或 pairs/BAM文件)
  -> [hicpro] 解析 .mmapstat / .mpairstat / .mRSstat
  -> [pairs/bam] BAM先转pairs -> 运行pairtools stats -> 解析统计
  -> 计算各质量指标(比对率/valid pairs/各异常比例/cis-trans)
  -> 逐项对照阈值 -> 输出PASS/FAIL报告 + 整体结论
```

## 输出 | Output

评估结果主要打印到终端/日志，同时落盘少量文件：

```text
hic_qc_output/
├── hicpro_qc.log        # hicpro模式运行日志(含评估报告)
└── pairtools_qc.log     # pairs/bam模式运行日志
    # pairs/bam模式额外产出:
    #   <样本>_stats.tsv      pairtools stats 原始统计
    #   <样本>.pairs.gz       BAM转出的pairs文件(仅BAM输入)
```

关键输出说明：

- **评估报告**：逐项打印 `[通过|PASS] / [未通过|FAIL]` 与实测值、阈值，最后给出整体结论。
- **<样本>_stats.tsv**：pairtools 的原始统计，便于自行核对或二次分析。

## 结果解读 | Interpreting Results

**通俗理解|In plain words:** 每一项都对应一个「数据质量问题」，未通过的那项就是最该关注的地方。

- **比对率低（<70%）**：读段质量差、基因组与样本不匹配、或比对参数问题。
- **唯一比对率低（<60%）**：重复序列多（如多倍体、转座子富集），读段大量多贴，影响后续分区。
- **valid pairs 低（<50%）**：大量读段落在异常类别（dangling/自连/重连），建库质量差。
- **dangling ends 高（>15%）**：酶切效率低或片段太短。
- **self-ligation 高（>5%）/ religation 高（>10%）**：建库连接步骤异常。
- **cis/trans 比例低（<5）**：染色体内信号弱，可能是跨染色体噪声多或样本有污染。
- **PCR 重复率高（>30%）**：文库复杂度低，扩增过度。

整体任一指标 FAIL 即判未通过，此时应优先排查对应问题再决定是否继续下游挂载。

## 参数选择建议 | Parameter Guidance

- **输入类型**：跑过 HiC-Pro 就用 hicpro（默认）；只有比对后 pairs/BAM 就用 pairs/bam。
- **阈值调整**：不同物种（尤其是高重复基因组）唯一比对率会偏低，可按项目经验放宽个别阈值；但比对率、valid pairs、cis/trans 建议保持默认。
- **BAM 输入**：务必提供正确的 chrom.sizes（两列：染色体名 + 长度），否则转 pairs 会失败。

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input, -i` | 必填 |  | 输入文件或目录｜Input file or directory (hicpro: HiC-Pro输出目录; pairs/bam: pairs或BAM文件) |
| `--input-type` | `hicpro` | hicpro/pairs/bam | 输入类型｜Input type (hicpro: HiC-Pro输出目录; pairs: .pairs.gz文件; bam: .bam文件) |
| `--output-dir, -o` | `./hic_qc_output` |  | 输出目录｜Output directory |
| `--sample-name, -s` | — |  | 样本名称（仅hicpro模式，可选，默认自动检测）｜Sample name (hicpro mode only, optional, auto-detect by default) |
| `--min-mapping-rate` | `70.0` | float | 最低比对率阈值%%(仅hicpro模式)｜Minimum mapping rate threshold%% (hicpro mode only) |
| `--min-unique-rate` | `60.0` | float | 最低唯一比对率阈值%%(仅hicpro模式)｜Minimum unique mapping rate threshold%% (hicpro mode only) |
| `--min-valid-pairs-rate` | `50.0` | float | 最低valid pairs比例阈值%%(仅hicpro模式)｜Minimum valid pairs rate threshold%% (hicpro mode only) |
| `--max-dangling-ends-rate` | `15.0` | float | 最高dangling ends比例阈值%%(仅hicpro模式)｜Maximum dangling ends rate threshold%% (hicpro mode only) |
| `--max-self-ligation-rate` | `5.0` | float | 最高self-ligation比例阈值%%(仅hicpro模式)｜Maximum self-ligation rate threshold%% (hicpro mode only) |
| `--max-religation-rate` | `10.0` | float | 最高religation比例阈值%%(仅hicpro模式)｜Maximum religation rate threshold%% (hicpro mode only) |
| `--pairtools-path, -p` | `~/miniforge3/envs/hic/bin/pairtools` |  | Pairtools可执行文件路径（仅pairs/bam模式）｜Pairtools executable path (pairs/bam mode only) |
| `--chroms-path, -c` | — |  | Chromosome sizes文件路径（BAM输入时必需）｜Chromosome sizes file path (required for BAM input) |
| `--max-unmapped-rate` | `20.0` | float | 未比对reads比例阈值%%(仅pairs/bam模式)｜Threshold for unmapped reads rate%% (pairs/bam mode only) |
| `--max-single-sided-rate` | `10.0` | float | 单端比对比例阈值%%(仅pairs/bam模式)｜Threshold for single-sided mapping rate%% (pairs/bam mode only) |
| `--min-mapped-rate` | `80.0` | float | 双端比对率阈值%%(仅pairs/bam模式)｜Threshold for paired mapping rate%% (pairs/bam mode only) |
| `--max-dup-rate` | `30.0` | float | PCR重复率阈值%%(仅pairs/bam模式)｜Threshold for PCR duplication rate%% (pairs/bam mode only) |
| `--min-cis-trans-ratio` | `5.0` | float | 最低cis/trans比例阈值｜Minimum cis/trans ratio threshold |
| `--max-duplication-rate` | `30.0` | float | 最高PCR重复率阈值%%｜Maximum PCR duplication rate threshold%% |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- Python 3（HiC-Pro 模式纯解析，无额外软件依赖）
- pairtools（pairs/bam 模式；默认 `~/miniforge3/envs/hic/bin/pairtools`，通过 conda run 自动检测环境）
- HiC-Pro 输出目录（hicpro 模式，需 bowtie_results/ 与 hic_results/ 子目录）

## 常见问题 | FAQ

**Q1：支持断点续传吗？**
不支持。评估是轻量的一次性计算，重跑即重新评估；换阈值重跑只需重新执行命令。

**Q2：hicpro 模式样本名检测错了怎么办？**
用 `-s <样本名>` 手动指定。样本名需与 `.mmapstat` / `.mpairstat` / `.mRSstat` 文件名的前缀一致。

**Q3：报「缺少必要子目录 bowtie_results / hic_results」？**
说明输入目录不是完整的 HiC-Pro 输出。确认 `-i` 指向 HiC-Pro 的顶层输出目录（含这两个子目录）。

**Q4：BAM 输入报「需要提供 chrom.sizes 文件」？**
BAM 转 pairs 必须知道染色体长度。用 `-c genome.chrom.sizes` 提供（格式：`chr1  12345678` 两列）。

**Q5：为什么有些指标显示为 0（如 cis/trans、PCR 重复率）？**
HiC-Pro 模式下，cis/trans 和 PCR 重复率在当前版本未从统计文件中提取，会显示为 0 并自动按「通过」处理；这两个指标在 pairs/bam 模式下才会真实计算。

