# rMVP 批量 GWAS 分析 | rMVP Batch GWAS Analysis

一句话理解：**把「哪个基因位点决定了某个性状」这件事，用 R 包 rMVP 的三个统计模型(GLM/MLM/FarmCPU)自动批量扫出来**。输入一个 VCF(群体基因型)和一个表型文件(每行一个样本的性状值)，输出每个位点的显著性排行表 + 曼哈顿图/QQ图。

## 功能概述 | Overview

- 封装 R 包 rMVP，支持三种 GWAS 模型：`GLM`(快但粗糙)、`MLM`(校正亲缘关系，稳)、`FarmCPU`(快而准，大规模首选)
- 支持多表型批量分析：表型文件可以有多个性状列，一次全跑
- 内置 LD 去连锁(PLINK)：亲缘矩阵/PCA 在去连锁后的 SNP 上算，GWAS 用全部 SNP
- 断点续传：已完成的步骤自动跳过(换参数重跑需删旧产物，见 FAQ)
- 自动整合结果：多模型多表型的 P 值合并成一张总表，并导出 ldblockshow 专用 TSV

## 快速开始 | Quick Start

```bash
biopytools rmvp -i input.vcf.gz -p phenotype.txt -o output
```

最小输入：一个 VCF(.vcf/.vcf.gz，压缩格式会自动解压)+ 一个 TSV 表型文件(第一列样本名，后面是性状值)。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| GWAS | 全基因组扫描，找出「哪些位点与性状显著相关」 |
| GLM | 一般线性模型，最朴素：只比「有这个碱基 vs 没这个碱基」的差异，快但不看「亲戚关系」 |
| MLM | 混合线性模型：额外考虑亲缘关系，像「把表兄弟之间的相似性先扣除」，假阳性更少 |
| FarmCPU | 固定/随机效应循环，兼顾速度和准确度，大数据集首选 |
| 亲缘关系矩阵(K) | 个体间基因相似度表格，用来校正「长得像是因为有血缘」 |
| 主成分(PCA/PC) | 群体结构的「主轴」，作为协变量放进模型排除分层干扰 |
| LD 去连锁(pruning) | 相邻位点往往「绑定遗传」，信息重复；去连锁=每个「团」只留一个代表来算 K/PCA，省计算 |
| MAF | 「少数派碱基」占比，太低=位点没信息量 |
| 曼哈顿图 / QQ图 | 前者把每位点 -log10(P) 画成天际线(塔高=显著)，后者诊断整体信号是否正常 |
| 显著性阈值 | 越过这条线的位点才算「显著」 |

## 输入 | Input

### VCF 文件

标准 VCF 格式，支持 `.vcf` 和 `.vcf.gz`。注意：rMVP 底层的 VCF 解析**不支持 gzip**，所以 .gz 会先用 `zcat` 解压到输出目录再分析(无需手动处理)。样本名须与表型文件一致。

### 表型文件

TSV 制表符分隔，**第一列样本名，后面可以有多列性状**(支持多表型批量分析)，必须有表头：

```text
SampleID    Height    Yield
Sample1     0.85      1.2
Sample2     0.72      0.9
Sample3     0.91      NA
```

- 表型值为数值，缺失用 `NA`。
- 有几列性状，就会跑「性状数 × 模型数」次分析(默认 3 模型)。

## 参数说明 | Parameters

### 必需参数 | Required

| 参数 | 说明 |
|------|------|
| `-i, --vcf` | 输入 VCF 文件 |
| `-p, --pheno` | 表型文件(TSV: 样本ID + 性状列) |
| `-o, --output` | 输出目录 |

### 模型与计算 | Models & computing

**通俗理解|In plain words:** `--models` 选跑哪几个模型(默认三个全跑，结果互相对照)。`--ncpus` 是并行核数；`--maxLine` 是每次读多少个 SNP——**内存不够时调小它**(如 5000)，一般不用动。

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--models` | `GLM MLM FarmCPU` | 分析模型(可多选) |
| `--ncpus` | `12` | CPU 核心数 |
| `--maxLine` | `10000` | 每次读取的 SNP 数(影响内存) |
| `--output-prefix` | `RMVP_Result` | 输出前缀 |

### PCA 协变量 | PCA covariates

**通俗理解|In plain words:** 每个模型用几个主成分当协变量来校正群体结构。**默认 3 个，一般不用动**；群体结构特别复杂(如多地区混合)时可适当加。

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--n-pc-glm` | `3` | GLM 使用的 PC 数 |
| `--n-pc-mlm` | `3` | MLM 使用的 PC 数 |
| `--n-pc-farmcpu` | `3` | FarmCPU 使用的 PC 数 |

### 模型内部参数 | Model internals

**通俗理解|In plain words:** 这些是 rMVP 模型内部的「引擎参数」，**默认值经过实践验证，一般不用动**。`--vc-method` 决定 MLM 怎么估方差组分；`--max-loop` 和 `--method-bin` 是 FarmCPU 的迭代与分箱方式。

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--vc-method` | `BRENT` | MLM 方差组分方法(BRENT/EMMA/HE) |
| `--max-loop` | `10` | FarmCPU 最大迭代次数 |
| `--method-bin` | `static` | FarmCPU bin 方法(static/fast-lmm) |

### 位点过滤 | Marker filtering

**通俗理解|In plain words:** 默认**不过滤**(`--maf`/`--miss` 为空)。数据质量差时再给值，如 `--maf 0.05 --miss 0.1` 表示删掉「少数派占比 <5%」和「缺失 >10%」的位点。

### LD 去连锁 | LD pruning

**通俗理解|In plain words:** 默认开启，用 PLINK 去连锁后算亲缘矩阵和 PCA，GWAS 仍用全部 SNP。**默认参数(3000kb 窗口、r²=0.2)一般不用动**；位点特别密时可调小窗口。

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--ld-pruning/--no-ld-pruning` | 开 | LD 去连锁开关 |
| `--ld-window` | `3000kb` | LD 修剪窗口 |
| `--ld-step` | `1` | LD 修剪步长 |
| `--ld-r2` | `0.2` | LD r² 阈值 |
| `--plink-path` | 自动 | PLINK 路径(默认 conda env pop) |

### 输出控制 | Output control

**通俗理解|In plain words:** 控制图的格式和显著性线。`--file-type` 选图片格式，`--dpi` 设清晰度，`--threshold` 是「显著」的判定线(默认 0.05，可设 1e-6 等更严)。

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--file-type` | `jpg` | 图片格式 jpg/pdf/tiff |
| `--dpi` | `300` | 图片分辨率 |
| `--threshold` | `0.05` | 显著性阈值 |

### R 环境 | R environment

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--r-env` | `rMVP` | R conda 环境名或路径 |
| `--r-path` | 无 | R 可执行文件路径(直调模式) |

## 分析流程 | Pipeline

**通俗理解|In plain words:** 先检查依赖 → 用 PLINK 去连锁 → 转成 rMVP 格式并算亲缘矩阵/PCA → 逐个模型跑 GWAS → 把结果整合成总表。

```text
输入 VCF + 表型
    │
    ▼
检查依赖(R/rMVP 包/PLINK) → 验证输入
    │
    ▼
LD 去连锁(PLINK，默认开): 去连锁 SNP 用于算 K/PCA
    │
    ▼
数据转换(rMVP MVP.Data): VCF → 基因型矩阵 + K + PCA
    │
    ▼
批量 GWAS: 逐表型逐模型跑 GLM/MLM/FarmCPU
    │
    ▼
结果整合: 合并 P 值总表 + 汇总报告 + ldblockshow TSV
```

## 输出 | Output

输出全部平铺在 `--output` 目录下(无子目录)。记 `{prefix}`=输出前缀(默认 `RMVP_Result`)、`{trait}`=性状名、`{model}`=GLM/MLM/FarmCPU：

```text
output/
├── {prefix}.log                        # 主日志
├── {prefix}.geno.desc / .phe / .geno.map   # rMVP 基因型矩阵+表型+位点图
├── {prefix}_pruned.kin.desc / .pc.desc     # 亲缘矩阵/PCA(去连锁 SNP 上算)
├── {prefix}_convert.R / {prefix}_batch.R   # 生成的 R 脚本
├── {trait}.{model}.{trait}.csv             # 每位点 P 值(全量结果)
├── {trait}.{model}_signals.{trait}.csv     # 显著位点
├── {trait}.{model}.{PlotType}.{trait}.jpg  # 曼哈顿图/QQ图/环形图等
├── {prefix}_glm_integrated.txt             # GLM 整合表(每性状×模型一列 P 值)
├── {prefix}_mlm_integrated.txt / {prefix}_farmcpu_integrated.txt
├── {prefix}_all_models_integrated.txt      # 三模型合并总表
├── {prefix}_merged_{model}.csv             # 多表型时按显著合并(仅 >1 性状)
├── {trait}.{model}.tsv                     # ldblockshow 专用 3 列 TSV
└── {prefix}_summary_report.txt             # 汇总报告
```

## 结果解读 | Interpreting Results

### 1. 全量 P 值表(`{trait}.{model}.{trait}.csv`)

**通俗理解|In plain words:** 最核心的结果，一行一个位点。含 `SNP/CHROM/POS` 列，**最后一列是 P 值**(列名=性状.模型)。P 越小越显著；画成曼哈顿图，塔越高越像真信号。

### 2. 整合表(`{prefix}_all_models_integrated.txt`)

**通俗理解|In plain words:** 把三个模型、所有性状的 P 值拼成一张大表，方便横向对比同一个位点在不同模型/性状下是否都显著——**多模型多性状都显著的点最可信**。

### 3. 汇总报告(`{prefix}_summary_report.txt`)

**通俗理解|In plain words:** 一句话总结：跑了几个性状、几个模型、各出了多少显著位点和图。

### 4. 图(曼哈顿/QQ/环形)

**通俗理解|In plain words:** 曼哈顿图找「冒尖的塔」，QQ 图看「尾巴翘不翘」。QQ 图右上角翘起=有真信号；全程贴对角线=基本没信号。

## 参数选择建议 | Parameter Guidance

| 场景 | 建议 |
|------|------|
| 常规分析 | 全部默认，三模型一起跑互相对照 |
| 内存不足 | `--maxLine 5000`(减小单次读入的 SNP 数) |
| 样本大(>1万) | 用 `--models FarmCPU`，兼顾速度和准确 |
| 只做快速初筛 | `--models GLM`，最快 |
| 数据质量一般 | 加 `--maf 0.05 --miss 0.1` 先过滤 |
| 位点极密 | 调小 `--ld-window`(如 1000kb) |
| 想要 PDF 高清图 | `--file-type pdf` |

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--vcf, -i` | 必填 |  | 输入VCF文件｜Input VCF file |
| `--pheno, -p` | 必填 |  | 输入表型文件｜Input phenotype file |
| `--output, -o` | 必填 |  | 输出目录｜Output directory |
| `--output-prefix` | `RMVP_Result` |  | 输出前缀｜Output prefix |
| `--models` | `['GLM', 'MLM', 'FarmCPU']` | GLM/MLM/FarmCPU | 分析模型｜Analysis models (可多选｜can select multiple) |
| `--r-env` | `~/miniforge3/envs/rMVP` |  | R conda环境路径｜R conda environment path (e.g., ~/miniforge3/envs/rMVP) |
| `--r-path` | — |  | R可执行文件路径｜R executable path |
| `--ncpus` | `12` | int | CPU核心数｜Number of CPU cores |
| `--maxLine` | `10000` | int | 每次读取的SNP数量｜Number of SNPs to read at once (较小值减少内存｜smaller uses less memory) |
| `--n-pc-glm` | `3` | int | GLM模型使用的PC数量｜Number of PCs for GLM |
| `--n-pc-mlm` | `3` | int | MLM模型使用的PC数量｜Number of PCs for MLM |
| `--n-pc-farmcpu` | `3` | int | FarmCPU模型使用的PC数量｜Number of PCs for FarmCPU |
| `--vc-method` | `BRENT` | BRENT/EMMA/HE | MLM方差组分分析方法｜MLM variance component method |
| `--max-loop` | `10` | int | FarmCPU最大迭代次数｜FarmCPU max iterations |
| `--method-bin` | `static` | static/fast-lmm | FarmCPU bin方法｜FarmCPU bin method |
| `--maf` | — | float | 最小等位基因频率阈值｜Minor allele frequency threshold |
| `--miss` | — | float | 缺失率阈值｜Missing rate threshold |
| `--file-type` | `jpg` | jpg/pdf/tiff | 图片格式｜Figure format |
| `--dpi` | `300` | int | 图片分辨率｜Figure DPI |
| `--threshold` | `0.05` | float | 显著性阈值｜Significance threshold |
| `--ld-pruning/--no-ld-pruning` | `True` |  | LD去连锁（默认开启）：kinship/PCA在去连锁SNP上计算，GWAS用全部SNP｜LD pruning (default on): K/PCA on pruned SNPs, GWAS uses all SNPs |
| `--ld-window` | `3000kb` |  | LD修剪窗口｜LD pruning window (e.g. 3000kb or 500) |
| `--ld-step` | `1` | int | LD修剪步长｜LD pruning step size |
| `--ld-r2` | `0.2` | float | LD r2阈值｜LD r2 threshold |
| `--plink-path` | — |  | PLINK可执行文件路径｜PLINK executable path (default: conda env Population_genetics) |
| `--log-level` | `INFO` | DEBUG/INFO/WARN/ERROR | 日志级别｜Log level |
| `--quiet` | — |  | 静默模式｜Quiet mode |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- R + rMVP 包(conda 环境默认 `rMVP`，用 `--r-env` 指定；另需 `bigmemory` 包)
- PLINK(默认 `~/miniforge3/envs/pop/bin/plink`，可用 `PLINK_PATH` 或 `--plink-path` 覆盖)
- conda(调用 R 时用 `conda run` 包装)、zcat(解压 gzip VCF)

## 常见问题 | FAQ

**Q1：换参数重跑，结果为什么没变？**
断点续传按输出文件存在性判断。换过滤参数(如 `--maf`、`--ld-window`)重跑旧目录前，先删对应旧产物(如 `{prefix}.geno.desc`、`{prefix}_pruned.vcf`)，否则会复用旧结果。

**Q2：VCF 是 gz 压缩的能直接用吗？**
能。工具会自动 `zcat` 解压到输出目录再分析(因为 rMVP 底层不支持 gzip)，无需手动处理。

**Q3：表型文件能放多个性状吗？**
能。除第一列样本名外，每列一个性状，工具会「性状数 × 模型数」全部跑完，并自动合并多表型显著结果。

**Q4：kin/PCA 是在哪些 SNP 上算的？**
默认 LD 去连锁后，在去连锁的 SNP 上算亲缘矩阵和 PCA；GWAS 本身仍用全部 SNP。关掉去连锁(`--no-ld-pruning`)则全部 SNP 上算。

**Q5：结果文件名里 `{trait}` 重复出现正常吗？**
正常。rMVP 的文件名规则是 `{trait}.{model}.{trait}.csv`(memo 用性状名)，例如 `Height.GLM.Height.csv`。
