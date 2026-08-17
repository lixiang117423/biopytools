# SNP Index 计算与分析 | SNP Index Calculation and Analysis

一句话理解：**从 BSA 混池的 VCF 里，算出每个位点的 SNP index 和两个混池之间的 ΔSNP index，自动挑出极端位点与候选区域并出图**，用来定位控制性状的基因组区间。

## 功能概述 | Overview

- 计算 SNP index = Alt深度 / (Ref深度 + Alt深度)，ΔSNP index = 混池1的SNP index − 混池2的SNP index
- 三重过滤：最小测序深度、最小质量值(QUAL)、最小 mapping 质量(MQ)
- 四层分析：基本统计 → 极端位点 → 潜在目标区域 → 滑动窗口 + 置信区间
- 多张可视化：曼哈顿图、ΔSNP index 分布、样本对比、相关性、滑动窗口折线图
- 三种模式：完整流程（默认）/ 只计算 / 只分析已有结果

## 快速开始 | Quick Start

```bash
biopytools snp-index -i variation.filtered.snp.vcf.gz -o results/
```

最小输入：一个含 AD 字段、至少 2 个样本的 VCF。默认取 VCF 里的前两个样本作为两个「混池」。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| BSA（混池分离分析） | 把「性状相反的个体」各自混成一个池测序，比较两池在每个位点的差异，定位控制性状的基因 |
| SNP index | 某个位点上「替代等位基因」占多大比例；越接近 1，这个位点越偏替代型 |
| ΔSNP index | 两个池的 SNP index 之差；绝对值越接近 1，说明两池在该位点差异越大，越像与性状连锁 |
| 测序深度 | 这个位点被测到多少次；太浅 = 计数不可靠 |
| 滑动窗口 | 把染色体切成一段段「格子」再统计，抹平单点噪音，看趋势 |
| 置信区间 | 用统计方法给出「正常波动范围」，超过这个范围才算异常 |

## 输入 | Input

### VCF 文件

标准 VCF（`.vcf` / `.vcf.gz`），要求：

- 至少 2 个样本（默认取前两个；可用 `--sample-names pool1 pool2` 指定）
- FORMAT 列必须含 **AD**（参考/替代等位基因各自的深度），如 `GT:AD` → `0/1:12,30`

```text
#CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO  FORMAT    bulk1      bulk2
chr1    100  .   A    G    60    PASS    .     GT:AD     0/1:12,30  1/1:2,40
```

## 分析流程 | Pipeline

```text
输入 VCF
    │
    ▼
步骤1: 计算 SNP index（按深度/质量/MQ 过滤）
    │   产出 {prefix}_results.tsv
    ▼
步骤2: 结果分析（基本统计 + 极端位点 + 潜在区域）
    │   产出 {prefix}_extreme_sites.tsv / {prefix}_potential_regions.tsv
    ▼
步骤3: 滑动窗口分析 + 置信区间
    │   产出 {prefix}_sliding_windows.tsv / {prefix}_confidence_intervals.txt
    │   / {prefix}_candidate_regions.tsv
    ▼
步骤4: 可视化（曼哈顿/分布/相关性/滑动窗口折线图）
        产出多张 .png
```

## 输出 | Output

输出目录为平铺结构（所有文件直接放在 `-o` 指定的目录下，无子目录）：

```text
results/
├── snp_index_results.tsv              # 逐位点结果（核心）
├── snp_index_extreme_sites.tsv        # 极端 ΔSNP index 位点
├── snp_index_potential_regions.tsv    # 潜在目标区域
├── snp_index_sliding_windows.tsv      # 滑动窗口结果
├── snp_index_confidence_intervals.txt # 置信区间
├── snp_index_candidate_regions.tsv    # 候选区域
├── snp_index_comprehensive.png        # 综合分析图（4 联）
├── snp_index_manhattan.png            # 曼哈顿图
├── snp_index_delta_distribution.png   # ΔSNP index 分布
├── snp_index_snp_comparison.png       # 两样本 SNP index 对比
├── snp_index_correlation.png          # 两样本相关性
└── snp_index_sliding_window.png       # 滑动窗口折线图
```

`{prefix}` 默认是 `snp_index`，可用 `-p` 改名。`--create-multi-chrom-plot` 会额外生成 `{prefix}_multi_chrom_sliding.png`。

## 结果解读 | Interpreting Results

- **`{prefix}_results.tsv`（核心表）**：列依次为 染色体、位置、参考/替代碱基、样本1的 Ref/Alt 深度与 SNP index、样本2的 Ref/Alt 深度与 SNP index、ΔSNP_index
- **ΔSNP index 怎么读**：接近 0 = 两池在该位点无差异；接近 +1 = 池1 更偏替代型；接近 −1 = 池2 更偏替代型
- **`{prefix}_extreme_sites.tsv`**：|ΔSNP index| > `--extreme-threshold`（默认 0.8）的位点，是最可疑的单点
- **`{prefix}_potential_regions.tsv`**：连续 ≥ `--min-region-snps`（默认 5）个 |ΔSNP index| > `--region-threshold`（默认 0.5）的位点组成的目标区间
- **滑动窗口折线图**：看整条染色体的 ΔSNP index 走势，明显隆起并越过置信区间的山头就是候选区域；曼哈顿图里越高的点越显著
- **好坏判据**：候选区域里 ΔSNP index 峰值接近 ±1、且窗口内 SNP 数足够（不是孤点），才值得往下做功能验证

## 参数选择建议 | Parameter Guidance

**通俗理解|In plain words:** 默认参数即可跑出结果；只有深度阈值和区域判定阈值在数据质量特殊时才需要动。

- **`--min-depth / --min-quality / --min-mapping-quality`（过滤）**：默认 10 / 20 / 20，一般不用动；测序深度低的数据可把 `--min-depth` 降到 5，反之提高到 20 更严格
- **`--sample-names`**：VCF 里样本超过 2 个、且默认取前两个不对时，用它显式指定两个池的样本名
- **`--extreme-threshold / --region-threshold`**：默认 0.8 / 0.5 一般不用动；想更严格收严，想多看点放宽
- **`--window-size / --step-size`（滑动窗口）**：默认 1Mb 窗口 / 100kb 步长，一般不用动
- **`--calculate-only / --analyze-only`**：只想算不想画图用 `--calculate-only`；已有 `results.tsv` 只想重画图用 `--analyze-only -r 旧结果.tsv`
- **`--skip-visualization`**：没装 matplotlib 或只要表格时用

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input, -i` | — |  | 输入VCF文件路径｜Input VCF file path |
| `--output, -o` | `./snp_index_output` |  | 输出目录路径｜Output directory path |
| `--result-file, -r` | — |  | 已有结果文件路径(用于分析模式)｜Existing result file path (for analysis mode) |
| `--prefix, -p` | `snp_index` |  | 输出文件前缀｜Output file prefix |
| `--min-depth` | `10` | int | 最小测序深度｜Minimum sequencing depth |
| `--min-quality` | `20` | int | 最小质量值｜Minimum quality value |
| `--min-mapping-quality` | `20` | int | 最小mapping质量｜Minimum mapping quality |
| `--sample-names` | — | str | 指定样本名称(sample1 sample2)｜Specify sample names (sample1 sample2) |
| `--extreme-threshold` | `0.8` | float | 极端DeltaSNP index阈值｜Extreme DeltaSNP index threshold |
| `--region-threshold` | `0.5` | float | 区域检测阈值｜Region detection threshold |
| `--min-region-snps` | `5` | int | 区域最少SNP数量｜Minimum SNPs for region |
| `--max-region-gap` | `10000` | int | 区域最大gap(bp)｜Maximum gap in region (bp) |
| `--window-size` | `1000000` | int | 滑动窗口大小(bp)｜Sliding window size in bp |
| `--step-size` | `100000` | int | 滑动步长(bp)｜Sliding step size in bp |
| `--min-window-snps` | `5` | int | 窗口最少SNP数｜Minimum SNPs per window |
| `--confidence-level` | `0.95` | float | 置信水平｜Confidence level |
| `--disable-sliding-window-plot` | — |  | 禁用滑动窗口折线图｜Disable sliding window line plot |
| `--create-multi-chrom-plot` | — |  | 创建多染色体分离图｜Create multi-chromosome separated plot |
| `--calculate-only` | — |  | 只计算SNP index，不分析｜Calculate SNP index only, no analysis |
| `--analyze-only` | — |  | 只分析已有结果，不计算｜Analyze existing results only, no calculation |
| `--skip-visualization` | — |  | 跳过可视化｜Skip visualization |
| `--threads, -t` | `12` | int | 线程数｜Number of threads |
| `--force, -f` | — |  | 强制覆盖已存在文件｜Force overwrite existing files |
| `--log-file` | — | str | 日志文件路径｜Log file path |
| `--verbose, -v` | — |  | 详细输出模式｜Verbose mode (-v: INFO, -vv: DEBUG) |
| `--quiet` | — |  | 静默模式｜Quiet mode (ERROR only) |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | — |  | 输入VCF文件路径｜Input VCF file path (required for calculation) |
| `-o, --output` | `./snp_index_output` |  | 输出目录路径｜Output directory path |
| `-p, --prefix` | `snp_index` |  | 输出文件前缀｜Output file prefix |
| `-r, --result-file` | — |  | 已有结果文件路径｜Existing result file path (for analysis only) |
| `--min-depth` | `10` | int | 最小测序深度｜Minimum sequencing depth |
| `--min-quality` | `20` | int | 最小质量值｜Minimum quality value |
| `--min-mapping-quality` | `20` | int | 最小mapping质量｜Minimum mapping quality |
| `--sample-names` | — |  | 指定要分析的样本名称｜Specify sample names to analyze |
| `--extreme-threshold` | `0.8` | float | 极端ΔSNP index阈值｜Extreme ΔSNP index threshold |
| `--region-threshold` | `0.5` | float | 区域检测阈值｜Region detection threshold |
| `--min-region-snps` | `5` | int | 区域最少SNP数量｜Minimum SNPs for region |
| `--max-region-gap` | `10000` | int | 区域最大gap(bp)｜Maximum gap in region |
| `--window-size` | `1000000` | int | 滑动窗口大小(bp)｜Sliding window size in bp |
| `--step-size` | `100000` | int | 滑动步长(bp)｜Sliding step size in bp |
| `--min-window-snps` | `5` | int | 窗口最少SNP数｜Minimum SNPs per window |
| `--confidence-level` | `0.95` | float | 置信水平｜Confidence level |
| `--disable-sliding-window-plot` | — | store_true | 禁用滑动窗口折线图｜Disable sliding window line plot |
| `--create-multi-chrom-plot` | — | store_true | 创建多染色体分离图｜Create multi-chromosome separated plot |
| `--calculate-only` | — | store_true | 只计算SNP index，不分析｜Calculate SNP index only, no analysis |
| `--analyze-only` | — | store_true | 只分析已有结果，不计算｜Analyze existing results only, no calculation |
| `--skip-visualization` | — | store_true | 跳过可视化｜Skip visualization |
| `-v, --verbose` | `0` | count | 详细输出模式｜Verbose mode (-v: INFO, -vv: DEBUG) |
| `--quiet` | — | store_true | 静默模式｜Quiet mode (ERROR only) |
| `--log-file` | — |  | 日志文件路径｜Log file path |
| `-f, --force` | — | store_true | 强制覆盖已存在文件｜Force overwrite existing files |
| `-t, --threads` | `1` | int | 线程数｜Number of threads |
| `-V, --version` | — | version |  |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- Python 库：numpy（分析统计）、matplotlib（可视化；缺失时自动跳过绘图，不影响表格输出）

## 常见问题 | FAQ

**Q1：报「样本数不足 2 个」？**
VCF 的 FORMAT 后至少要 2 个样本列。默认取前两个样本，若你真正要分析的池不是前两个，请用 `--sample-names` 指定。

**Q2：为什么很多位点没进结果？**
位点会被三重过滤丢弃：任一池深度 < `--min-depth`、QUAL < `--min-quality`、MQ < `--min-mapping-quality`。日志会统计低质量/低深度位点数，据此可判断阈值是否过严。

**Q3：没生成任何图片？**
matplotlib 未安装时会跳过可视化（日志有告警），但 TSV 结果不受影响。安装 matplotlib 后重跑，或用 `--analyze-only -r 结果.tsv` 只重画图。

**Q4：`--analyze-only` 为什么报错？**
该模式必须有 `-r` 指定已有的 `results.tsv`（TSV 格式、表头含 `_SNP_index` 列），否则无法加载数据。
