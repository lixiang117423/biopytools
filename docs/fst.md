# Fst 遗传分化计算 | Fst Calculation (PLINK)

一句话理解：**计算群体之间的遗传分化系数(Fst)**，回答「这两个(或多个)群体在遗传上有多不一样」，是衡量种群结构、群体历史的基础指标。

## 功能概述 | Overview

- 基于 PLINK 的 Weir & Cockerham 方法计算群体间 Fst
- 输入 VCF + 群体文件，一键完成 VCF→PLINK→Fst 全流程
- 支持可选质控(MAF/缺失率/HWE)、LD 剪枝、SNP 抽稀
- 支持两群体与多群体两两比较
- 样本量不均衡时自动启用 Bootstrap 抽样，估计 Fst 均值与置信区间
- 输出长格式表格、矩阵格式与 PLINK 原始结果

## 快速开始 | Quick Start

```bash
biopytools fst -i variants.vcf -p population.txt -o fst_output
```

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| Fst(固定指数) | 群体间遗传差异的占比打分，0=完全一样，1=完全分化 |
| 遗传分化 | 两个群体各自「抱团」、基因库不再互通的程度 |
| MAF | 「少数派」等位基因在群体里的频率；太低=这个位点没区分度 |
| 缺失率 | 多少样本在这个位点「没测出来」；太高=判断依据不足 |
| HWE 平衡 | 群体随机交配时基因型应有的比例；偏离太远提示分型错误 |
| LD(连锁不平衡) | 相邻位点「绑定遗传」，信息重复，降维时每个团留一个代表 |
| Bootstrap | 从现有样本里反复「抓阄重抽样」多次，看结果稳不稳 |

## 输入 | Input

### VCF 文件

标准 VCF 格式，含基因型(GT)：

```text
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1	sample2	sample3	sample4
chr1	1000	.	A	G	.	PASS	.	GT	0/0	0/1	1/1	0/0
```

### 群体文件

两列：`样本ID` `群体标签`，支持制表符、逗号、空格分隔：

```text
sample1	POP1
sample2	POP1
sample3	POP2
sample4	POP2
```

样本名必须与 VCF 头部的样本名完全一致。

## 分析流程 | Pipeline

```text
VCF + 群体文件
    │
    ▼
步骤1: VCF → PLINK 格式(.bed/.bim/.fam)
    │
    ▼
步骤2: 质控(默认跳过；--enable-qc 才做 MAF/缺失/HWE 过滤)
    │
    ▼
步骤3: LD 剪枝(默认开启；SNP < 50万时自动跳过)
    │
    ▼
步骤4: PLINK 计算 Fst(两群体单次；多群体两两比较)
    │
    ├─ 直接模式 → 长格式表 + 矩阵表
    └─ 样本量不均衡 → Bootstrap 抽样 + 汇总均值/置信区间
```

## 输出 | Output

直接计算模式：

```text
fst_output/
├── fst_result.fst                # PLINK 原始 Fst 输出(逐位点)
├── population_filtered.plink     # 过滤后的群体文件
├── fst_long_format.txt           # 长格式(Population1 Population2 Fst)
├── fst_matrix.txt                # 矩阵格式(群体×群体)
├── 00_intermediate/              # 中间 PLINK 文件(默认保留)
│   ├── mydata.bed/bim/fam        # VCF 转换结果
│   ├── mydata_qc.bed/bim/fam     # 质控后(启用质控时)
│   └── mydata_pruned.bed/bim/fam # LD 剪枝后(SNP≥50万且启用时)
└── 99_logs/fst_calculation.log   # 运行日志
```

Bootstrap 模式(样本量不均衡自动触发，或 `--enable-bootstrap`)：

```text
fst_output/
├── bootstrap_iterations/         # 每次迭代的中间结果
│   ├── fst_iter1.fst
│   └── ...
├── all_fst_results.txt           # 汇总(Mean_Fst/Std/Min/Max/Median，核心结果)
└── 99_logs/fst_calculation.log
```

- `fst_long_format.txt`：三列 `Population1 Population2 Fst`，一行一对群体，便于统计
- `fst_matrix.txt`：对角线为 0 的对称矩阵，便于画热图
- `all_fst_results.txt`：Bootstrap 汇总，`Mean_Fst` 是主要结果，`Std` 是置信区间

## 结果解读 | Interpreting Results

**通俗理解|In plain words:** Fst 越接近 0 说明两个群体「基因库几乎混在一起」，越接近 1 说明「彻底分家」。粗略参考(物种间有差异)：

| Fst 范围 | 分化程度 |
|----------|----------|
| < 0.05 | 很小，几乎没分化 |
| 0.05–0.15 | 中等分化 |
| 0.15–0.25 | 分化较大 |
| > 0.25 | 分化很大 |

- `fst_result.fst`：每个 SNP 位点的 Fst，画成曼哈顿图可找分化热点区域
- `all_fst_results.txt` 的 `Mean_Fst` 用于汇报「某两群体 Fst = 0.XX」；`Std` 越小结果越稳
- Fst 受群体历史、迁移率、选择影响，是相对指标，比较时保持样本量口径一致

## 参数选择建议 | Parameter Guidance

- **质控默认关闭**：只有加 `--enable-qc` 才会做 MAF/缺失/HWE 过滤；阈值 `--maf 0.05`/`--geno 0.1`/`--mind 0.1`/`--hwe 1e-6` 一般不用动
- **LD 剪枝默认开启**，但 SNP 少于 50 万时自动跳过；窗口 `--ld-window 50`/`--ld-step 10`/`--ld-r2 0.2` 一般不用动
- **Bootstrap**：样本量不均衡会自动启用；想让结果带置信区间可手动 `--enable-bootstrap`，`--bootstrap-iterations` 默认 100 已够用
- **`--min-samples 10`**：样本少于 10 的群体自动排除；想保留小群体可调小，想排除某群体用 `--exclude-pops POP1,POP2`
- **`--thin`**：SNP 特别多时抽稀(0-1 比例)，一般用不到
<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --vcf-file` | 必填 |  | VCF文件路径｜VCF file path |
| `-p, --pop-file` | 必填 |  | 群体文件路径（样本ID + 群体标签）｜Population file path (sample ID + population label) |
| `-o, --output-dir` | `./fst_output` | Path | 输出目录｜Output directory |
| `--plink-path` | — |  | PLINK软件路径｜PLINK software path (default: auto-detect) |
| `--maf` | `0.05` | float | 最小等位基因频率阈值｜Minor allele frequency threshold |
| `--geno` | `0.1` | float | 位点缺失率阈值｜Genotype missing rate threshold |
| `--mind` | `0.1` | float | 样本缺失率阈值｜Sample missing rate threshold |
| `--hwe` | `1e-06` | float | Hardy-Weinberg平衡p值阈值｜Hardy-Weinberg equilibrium p-value threshold |
| `--enable-qc` | — |  | 启用质控过滤（默认禁用）｜Enable quality control filtering |
| `--no-keep-intermediate` | — |  | 不保留中间文件｜Do not keep intermediate files |
| `--enable-bootstrap` | — |  | 启用bootstrap抽样｜Enable bootstrap sampling |
| `--bootstrap-iterations` | `100` | int | Bootstrap迭代次数｜Bootstrap iterations |
| `--min-samples` | `10` | int | 最小样本数阈值｜Minimum sample count threshold |
| `--exclude-pops` | — |  | 手动排除群体（逗号分隔）｜Manually exclude populations (comma-separated) |
| `--no-ld-prune` | — |  | 禁用LD pruning｜Disable LD pruning |
| `--ld-window` | `50` | int | LD pruning窗口大小｜LD pruning window size |
| `--ld-step` | `10` | int | LD pruning步长｜LD pruning step size |
| `--ld-r2` | `0.2` | float | LD pruning R2阈值｜LD pruning R2 threshold |
| `--thin` | — | float | SNP抽稀比例｜SNP thinning ratio |
| `-t, --threads` | `12` | int | 线程数｜Number of threads |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --vcf-file` | 必填 |  | VCF文件路径｜VCF file path |
| `-p, --pop-file` | 必填 |  | 群体文件路径（样本ID + 群体标签）｜Population file path (sample ID + population label) |
| `-o, --output-dir` | `./fst_output` |  | 输出目录｜Output directory (default: ./fst_output) |
| `--plink-path` | — |  | PLINK软件路径｜PLINK software path (default: auto-detect) |
| `--enable-qc` | — | store_true | 启用质控过滤（默认禁用）｜Enable quality control filtering (disabled by default) |
| `--maf` | `0.05` | float | 最小等位基因频率阈值｜Minor allele frequency threshold (default: 0.05) |
| `--geno` | `0.1` | float | 位点缺失率阈值｜Genotype missing rate threshold (default: 0.1) |
| `--mind` | `0.1` | float | 样本缺失率阈值｜Sample missing rate threshold (default: 0.1) |
| `--hwe` | `1e-06` | float | Hardy-Weinberg平衡p值阈值｜Hardy-Weinberg equilibrium p-value threshold (default: 1e-6) |
| `--no-keep-intermediate` | — | store_true | 不保留中间文件｜Do not keep intermediate files |
| `--enable-bootstrap` | — | store_true | 启用bootstrap抽样｜Enable bootstrap sampling |
| `--bootstrap-iterations` | `100` | int | Bootstrap迭代次数｜Bootstrap iterations (default: 100) |
| `--min-samples` | `10` | int | 最小样本数阈值，排除样本数少于此值的群体｜Minimum sample count threshold, exclude populations with fewer samples (default: 10) |
| `--exclude-pops` | — |  | 手动指定要排除的群体（逗号分隔）｜Manually specify populations to exclude (comma-separated) |
| `--no-ld-prune` | — | store_true | 禁用LD pruning｜Disable LD pruning |
| `--ld-window` | `50` | int | LD pruning窗口大小｜LD pruning window size (default: 50) |
| `--ld-step` | `10` | int | LD pruning步长｜LD pruning step size (default: 10) |
| `--ld-r2` | `0.2` | float | LD pruning R2阈值｜LD pruning R2 threshold (default: 0.2) |
| `--thin` | — | float | SNP抽稀比例（0-1之间）｜SNP thinning ratio (between 0-1) |
| `--threads` | `12` | int | 线程数｜Number of threads (default: 12) |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- PLINK(默认 `~/miniforge3/envs/pop/bin/plink`，可用 `--plink-path` 或环境变量 `PLINK_PATH` 覆盖)

## 常见问题 | FAQ

**Q1：质控为什么没生效？**
本模块默认**不做质控**。只有加 `--enable-qc` 才会启用 MAF/缺失/HWE 过滤；不加时 `--maf` 等阈值被忽略。

**Q2：LD 剪枝没执行？**
SNP 数少于 50 万时自动跳过 LD 剪枝(此时剪枝意义不大)。想强制关闭用 `--no-ld-prune`。

**Q3：样本名对不上会怎样？**
群体文件里写了但 VCF 里没有的样本无法匹配，可能被忽略或报错。先核对两边样本名一致(可用 `bcftools query -l variants.vcf` 对照)。

**Q4：换参数重跑要删旧文件吗？**
本模块无断点续传，每次运行重算并覆盖同名输出。换 `--maf`、`--enable-qc` 等重跑同一 `-o` 即可；想保留多组结果请换输出目录。

**Q5：为什么自动进了 Bootstrap 模式？**
当各群体样本量严重不均衡(直接算 Fst 偏差大)时自动启用 Bootstrap 抽样。不想用可核对群体样本数，或改用 `--exclude-pops` 剔除极小的群体。