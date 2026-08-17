# 选择性扫荡检测 | Selective Sweep Detection

一句话理解：**输入一份 VCF 和群体分组信息，自动计算多种「扫荡信号」统计量并合成一个综合分，把最像受过自然/人工选择的基因组区域排出来**。

## 功能概述 | Overview

- 一条命令跑完：过滤 → 按群体拆分 → 逐群体/两两群体算统计 → 合并打分 → 输出候选区域
- 集成 6 类统计量：π、Tajima's D（群体内）、Fst、XP-CLR（群体间）、RAiSD μ、SweeD CLR
- 各统计量做百分位排名后合成 composite_score，top 分位数窗口自动合并为候选区域
- 样本量偏少的群体，其噪声大的分量（μ/CLR/XP-CLR）默认自动从打分中排除
- 断点续传：已完成的步骤自动跳过

## 快速开始 | Quick Start

```bash
biopytools selective-sweep -i variants.vcf.gz -p pops.txt -o sweep_output
```

最小输入：一个 VCF（`.vcf` / `.vcf.gz`）+ 一个群体分组文件（每行「样品ID<TAB>分组」，无表头）。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| 选择性扫荡 | 某个有利变异被选择后，它附近整段 DNA 一起「搭顺风车」固定下来，多样性被拉低，留下可检测的信号 |
| π（核苷酸多样性） | 群体内个体之间平均有多大差异；扫荡区 π 会异常低 |
| Tajima's D | 用「位点差异数」与「位点种类数」是否失衡来判断是否偏离中性；扫荡区明显为负 |
| Fst | 两个群体分化程度；扫荡区两群差异被拉开，Fst 升高 |
| RAiSD μ | 一种快速扫描统计量，综合多样性、SFS、LD 三种信号 |
| SweeD CLR | 复合似然比(CLR)，衡量某区域「经历过扫荡」比「中性」的可能性高出多少 |
| XP-CLR | 跨群体的 CLR，像 Fst 一样做两两群体比较 |
| 窗口/步长 | 把染色体切成固定大小的「格子」来统计；窗口=格子宽度，步长=每格往前挪多少 |
| composite_score | 把上面所有信号各自排名后取平均的综合分，越高越像扫荡区 |

## 输入 | Input

### VCF 文件

标准 VCF（`.vcf` / `.vcf.gz`），样本名须与群体分组文件一致。

### 群体分组文件

每行两列，TAB 分隔：第一列样品 ID，第二列分组名，**无表头**：

```text
sample1    wild
sample2    wild
sample3    cultivated
sample4    cultivated
```

- 样品名拼写错误、带表头行会被自动跳过并告警；某群体全部样本都不在 VCF 中时该群体被跳过
- 至少 1 个有效群体才能运行；群体数 ≥2 时才会计算 Fst 和 XP-CLR

## 分析流程 | Pipeline

```text
输入 VCF + 群体分组文件
    │
    ▼
步骤1: bcftools 过滤(双等位 SNP + MAF + 缺失率)
    │   产出 01_filter/filtered.vcf.gz
    ▼
步骤2: 按群体拆分样本列表与子 VCF
    │   产出 01_filter/{pop}.samples.txt / {pop}.vcf.gz
    ▼
步骤3: 每群体统计 π / Tajima's D / RAiSD μ / SweeD CLR
    │   产出 02_stats/{pop}.windowed.pi 等
    ▼
步骤4: 两两群体 Fst / XP-CLR(逐染色体)
    │   产出 02_stats/{a}_{b}.windowed.weir.fst / XPCLR_{a}_{b}.{chr}.tsv
    ▼
步骤5: 各统计量百分位排名 → composite_score
    │   产出 03_sweep/genome_wide_stats.tsv
    ▼
步骤6: top 分位数窗口合并 → 候选区域
        产出 03_sweep/candidate_regions.tsv
```

## 输出 | Output

```text
sweep_output/
├── 00_pipeline_info/
│   ├── pop_summary.tsv              # 各群体样本数
│   └── software_versions.yml        # 软件版本与参数
├── 01_filter/
│   ├── filtered.vcf.gz              # 过滤后的 VCF（含索引）
│   ├── {pop}.samples.txt            # 每群体的样本列表
│   └── {pop}.vcf.gz                 # 每群体的子 VCF
├── 02_stats/
│   ├── {pop}.windowed.pi            # vcftools 窗口 π
│   ├── {pop}.Tajima.D               # vcftools 窗口 Tajima's D
│   ├── RAiSD_Report.{pop}.{chr}     # RAiSD μ 报告
│   ├── SweeD_Report.{pop}.{chr}     # SweeD CLR 报告
│   ├── {a}_{b}.windowed.weir.fst    # 两两 Fst
│   └── XPCLR_{a}_{b}.{chr}.tsv      # 两两 XP-CLR（逐染色体）
├── 03_sweep/
│   ├── genome_wide_stats.tsv        # 全基因组逐窗口综合表（含 composite_score）
│   └── candidate_regions.tsv        # 候选扫荡区域（最终结果）
├── 99_logs/
│   └── *.log                        # 运行日志
└── tmp/                             # 中间文件（运行结束自动清理）
```

## 结果解读 | Interpreting Results

- **`candidate_regions.tsv`（最终结果）**：列为 CHR / START / END / MAX_SCORE / N_WIN。composite_score 越高越可能是扫荡区；一行就是一个被合并后的候选区间
- **`genome_wide_stats.tsv`（明细）**：每个窗口各统计量的原始值与排名列（PI_/TajD_/MU_/CLR_/XPCLR_/Fst_ 前缀），以及 composite_score 和 n_stats_supporting（有几个分量支持该窗口）。支持分量越多、分数越高越可信
- **π 越低、Tajima's D 越负、Fst/μ/CLR/XP-CLR 越高**，都是「扫荡」的一致方向；composite_score 综合了这些
- **好坏判据**：候选区域若同时被多个统计量支持（n_stats_supporting 高），且落在已知功能基因附近，结论更可信；单统计量孤峰要谨慎

## 参数选择建议 | Parameter Guidance

**通俗理解|In plain words:** 默认参数即可直接跑；大多数参数都不用动，只有窗口大小和候选阈值偶尔需要按数据规模调整。

- **`--win / --step`（窗口/步长）**：默认 50kb 适合中等密度数据；SNP 稀疏可加大到 100kb，SNP 极密可减小到 20kb。步长默认等于窗口（无重叠）
- **`--top-quantile`（候选阈值）**：默认 0.01 = 取分数最高的 1% 窗口做候选；想多看放宽到 0.05，想更严格收紧到 0.005
- **`--merge-gap`**：相邻候选窗口间隔 ≤ 此值就合并成一个区域，默认 100kb 一般不用动
- **`--min-maf / --max-missing`（过滤）**：默认 0.05 / 0.10，一般不用动；过滤后 VCF 为空会自动报错提示放宽
- **低样本排除（--raisd/sweed/xpclr-min-samples + --include-*-low-n）**：默认样本量 <15 的群体其对应分量自动排除；这是合理默认，一般不要强制加回
- **`--sweed-unfolded`**：默认折叠 SFS（无需祖先状态，正确）；只有确认有祖先状态时才用未折叠

## 依赖 | Dependencies

- bcftools（conda 环境 align）、vcftools（pop）、RAiSD（pop）、XP-CLR（pop）、SweeD（`~/software/sweed/SweeD-P`）
- Python 库：numpy、pandas（用于统计合并打分）

## 常见问题 | FAQ

**Q1：换参数重跑，为什么结果没变？**
断点续传按各步骤输出文件是否存在判断。换过滤参数（如 `--min-maf`）后，需删除 `01_filter/filtered.vcf.gz` 及下游的 `02_stats/`、`03_sweep/` 旧产物，否则会复用旧结果。

**Q2：pop_info 里的样本被报「不在 VCF 中」？**
bcftools 对不存在的样本会直接报错退出。请核对分组文件里的样本名是否与 VCF 表头（`#CHROM` 行）完全一致，且没有混入表头行。

**Q3：为什么某些群体的 μ/CLR/XP-CLR 没进 composite_score？**
样本量 < 阈值的群体，这些噪声大的分量默认被排除（日志会告警）。这是合理默认；确实需要时才用 `--include-mu-low-n` 等开关强制加入。

**Q4：只有一个群体能跑吗？**
能。单群体只算 π / Tajima's D / RAiSD μ / SweeD CLR；Fst 和 XP-CLR 需要 ≥2 个群体，会自动跳过。

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input, -i` | 必填 |  | 输入VCF文件｜Input VCF file |
| `--pop-info, -p` | 必填 |  | 群体分组文件(样品ID<TAB>分组,无表头)｜Population info file |
| `--output-dir, -o` | `./selective_sweep_output` | Path | 输出目录｜Output directory |
| `--threads, -t` | `12` | int | 线程数｜Number of threads |
| `--win` | `50000` | int | 窗口大小｜Window size |
| `--step` | `50000` | int | 窗口步长｜Window step |
| `--top-quantile` | `0.01` | float | 候选阈值分位数｜Candidate threshold quantile |
| `--merge-gap` | `100000` | int | 候选窗口合并最大间隔｜Max gap for merging candidate windows |
| `--min-maf` | `0.05` | float | 过滤MAF阈值｜Filter MAF threshold |
| `--max-missing` | `0.1` | float | 过滤缺失率阈值｜Filter missing rate threshold |
| `--raisd-window` | `50` | int | RAiSD SNP窗口｜RAiSD SNP window |
| `--raisd-min-samples` | `15` | int | 低样本量阈值(低于此值MU分量默认排除)｜Low sample threshold (MU excluded below) |
| `--include-mu-low-n` | `False` |  | 低样本群体也加入MU分量｜Include MU component for low-n pops |
| `--sweed-grid` | `10000` | int | SweeD CLR网格点数｜SweeD CLR grid points |
| `--sweed-min-samples` | `15` | int | 低样本量阈值(低于此值CLR分量默认排除)｜Low sample threshold (CLR excluded below) |
| `--include-sweed-low-n` | `False` |  | 低样本群体也加入CLR分量｜Include CLR component for low-n pops |
| `--sweed-unfolded` | `False` |  | 使用未折叠SFS(需祖先状态,默认折叠)｜Unfolded SFS (requires ancestral state; folded by default) |
| `--xpclr-maxsnps` | `200` | int | XP-CLR窗口最大SNP数｜XP-CLR max SNPs per window |
| `--xpclr-minsnps` | `10` | int | XP-CLR窗口最小SNP数｜XP-CLR min SNPs per window |
| `--xpclr-ld` | `0.95` | float | XP-CLR LD加权截断｜XP-CLR LD cutoff |
| `--xpclr-min-samples` | `15` | int | 低样本量阈值(任一群体低于此值XP-CLR分量默认排除)｜Low sample threshold (XP-CLR excluded if either pop below) |
| `--include-xpclr-low-n` | `False` |  | 低样本群体对也加入XP-CLR分量｜Include XP-CLR component for low-n pairs |
| `--log-level` | `INFO` |  | 日志级别｜Log level |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 输入VCF文件｜Input VCF file (支持.gz｜.gz supported) |
| `-p, --pop-info` | 必填 |  | 群体分组文件(样品ID<TAB>分组,无表头)｜Population info file (sample<TAB>group, no header) |
| `-o, --output-dir` | `./selective_sweep_output` |  | 输出目录｜Output directory |
| `-t, --threads` | `12` | int | 线程数(默认12)｜Thread count (default 12) |
| `--win` | `50000` | int | 窗口大小(默认50000)｜Window size (default 50000) |
| `--step` | `50000` | int | 窗口步长(默认50000)｜Window step (default 50000) |
| `--top-quantile` | `0.01` | float | 候选阈值分位数(默认0.01)｜Candidate threshold quantile (default 0.01) |
| `--merge-gap` | `100000` | int | 候选窗口合并最大间隔(默认100000)｜Max gap for merging candidate windows (default 100000) |
| `--min-maf` | `0.05` | float | 过滤MAF阈值(默认0.05)｜Filter MAF threshold (default 0.05) |
| `--max-missing` | `0.1` | float | 过滤缺失率阈值(默认0.10)｜Filter missing rate threshold (default 0.10) |
| `--raisd-window` | `50` | int | RAiSD SNP窗口(默认50)｜RAiSD SNP window (default 50) |
| `--raisd-min-samples` | `15` | int | 低样本量阈值,低于此值MU分量默认排除(默认15)｜Low sample threshold; MU excluded below (default 15) |
| `--include-mu-low-n` | — | store_true | 低样本群体也加入MU分量｜Include MU component for low-n pops |
| `--sweed-grid` | `10000` | int | SweeD CLR网格点数(默认10000)｜SweeD CLR grid points (default 10000) |
| `--sweed-min-samples` | `15` | int | 低样本量阈值,低于此值CLR分量默认排除(默认15)｜Low sample threshold; CLR excluded below (default 15) |
| `--include-sweed-low-n` | — | store_true | 低样本群体也加入CLR分量｜Include CLR component for low-n pops |
| `--sweed-unfolded` | — | store_true | 使用未折叠SFS(需祖先状态,默认折叠)｜Unfolded SFS (requires ancestral state; folded by default) |
| `--xpclr-maxsnps` | `200` | int | XP-CLR窗口最大SNP数(默认200)｜XP-CLR max SNPs per window (default 200) |
| `--xpclr-minsnps` | `10` | int | XP-CLR窗口最小SNP数(默认10)｜XP-CLR min SNPs per window (default 10) |
| `--xpclr-ld` | `0.95` | float | XP-CLR LD加权截断(默认0.95)｜XP-CLR LD cutoff (default 0.95) |
| `--xpclr-min-samples` | `15` | int | 低样本量阈值,任一群体低于此值XP-CLR分量默认排除(默认15)｜Low sample threshold; XP-CLR excluded if either pop below (default 15) |
| `--include-xpclr-low-n` | — | store_true | 低样本群体对也加入XP-CLR分量｜Include XP-CLR component for low-n pairs |
| `--log-level` | `INFO` |  | 日志级别(默认INFO)｜Log level (default INFO) |

<!-- END PARAMS:auto -->
