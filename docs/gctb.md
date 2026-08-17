# GCTB 复杂性状贝叶斯分析 | GCTB Bayesian Analysis

一句话理解：**从 VCF 出发，自动完成格式转换、质控、LD 矩阵构建，再用贝叶斯模型估算「每个位点对性状贡献多大、整个性状的遗传力有多高」**。支持个体基因型（个体水平）和 GWAS 汇总统计（汇总水平）两种输入。

## 功能概述 | Overview

- 一键式流程：VCF 转换 → 质控 → 频率 →（汇总模式）LD 矩阵 → 贝叶斯分析
- 三种贝叶斯模型：BayesS、BayesR、BayesC（--bayes-type）
- 两种分析模式：individual（个体基因型+表型）、summary（GWAS 汇总统计+LD 矩阵）
- 三种 LD 矩阵：sparse（稀疏，默认）、block（分块）、eigen（特征值分解）
- 批量模式（默认开启）自动拆分多表型，逐表型分析并汇总遗传力等参数
- 断点续传：已完成的步骤自动跳过

## 快速开始 | Quick Start

```bash
biopytools gctb -i variants.vcf -p phenotype.txt -o results
```

最小输入：一个 VCF + 一个表型文件（个体水平分析，FID/IID/PHENO 三列）。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| 遗传力(h2) | 性状差异里「能遗传」的比例，像「身高有多少是天生的」 |
| 贝叶斯模型 | 先假设「大多数位点效应很小、少数位点效应大」，再据此估计每个位点的真实效应 |
| BayesS / R / C | 三种贝叶斯模型，区别在「允许效应非零的位点比例」和「效应大小的分布」假设不同 |
| LD 矩阵 | 位点之间「绑定遗传」程度的表，汇总水平分析用来还原单个位点的信号 |
| 稀疏/分块 LD | 只存邻近位点或按块存，省内存；特征值分解(eigen)适合后续预测 |
| PIP | 后验包含概率，某位点「确实有非零效应」的概率，越大越可信 |
| PVE | 解释的表型方差比例，衡量位点/模型解释了性状多少变异 |

## 输入 | Input

### VCF 文件

标准 VCF（.vcf / .vcf.gz），程序转成 PLINK bed 格式。

### 表型文件（个体水平）

制表符分隔，FID、IID、PHENO 三列；也可用「sample + 多列性状」格式，批量模式会自动拆分：

```text
FID    IID    height
S1     S1     12.5
S2     S2     15.2
```

### GWAS 汇总统计文件（汇总水平，--analysis-mode summary）

包含 SNP、A1、A2、频率、BETA、SE、P、N 等列的汇总统计表，与 LD 矩阵配套使用。

## 分析流程 | Pipeline

```text
输入 VCF + 表型/汇总统计
    |
    v
1. VCF 转 PLINK bed(01_vcf_to_plink)
    |
    v
2. 质控(--maf/--geno/--hwe → 02_quality_control)
    |
    v
3. 计算等位基因频率(input_qc.frq)
    |
    v
4. [仅汇总模式] 构建 LD 矩阵(03_ld_matrix)
    |
    v
5. 贝叶斯分析(individual: --bayes; summary: --sbayes)
    |
    v
6. [批量模式] 汇总各表型的遗传力/候选位点(summary/)
```

## 输出 | Output

```text
results/
+-- 01_vcf_to_plink/input.bed/.bim/.fam   # PLINK 格式基因型
+-- 02_quality_control/input_qc.bed/.bim/.fam + input_qc.frq  # 质控后数据+频率
+-- 03_ld_matrix/ld_matrix.ldm.sparse      # LD 矩阵(仅汇总模式)
+-- 04_gctb_analysis/
|   +-- bayess.snpRes / bayess.parRes / bayess.log   # 个体水平结果
|   +-- (批量模式) 每个表型一个子目录,内含同名结果
+-- summary/                               # 批量汇总(仅批量模式)
|   +-- all_traits_parameters.tsv          # 各表型遗传力/多基因性等参数
|   +-- all_traits_top_snps.tsv            # 各表型按 PIP 排序的候选位点
|   +-- all_traits_all_snps.tsv            # 合并的全部位点结果
|   +-- analysis_report.md                 # 可读汇总报告
+-- converted_phenotypes/                  # 拆出的单表型文件(批量模式)
+-- 99_logs/                               # 日志
```

## 结果解读 | Interpreting Results

### 1. bayess.snpRes（每个位点的结果）

核心列：Name（SNP 名）、Chrom（染色体）、Position（位置）、A1Effect（A1 等位基因效应）、PIP（后验包含概率）。**PIP 越接近 1 越可能是真实效应位点**；按 PIP 排序取前几十个就是候选位点清单。

### 2. bayess.parRes（模型参数结果）

包含遗传力(hsq)、有效 SNP 数(NnzSnp)、遗传方差(GenVar)、残差方差(ResVar)等。**hsq 是核心结论**：性状遗传力高低；NnzSnp 反映「多少位点在起作用」，越大说明越「多基因」。

### 3. summary/（批量汇总）

批量模式下自动生成：all_traits_parameters.tsv 一张表对比所有性状的遗传力；all_traits_top_snps.tsv 直接给出每个性状的 top 候选位点；analysis_report.md 是给人读的总结。

## 参数选择建议 | Parameter Guidance

- --bayes-type：默认 S（BayesS，最常用）；R 适合「少数大效应位点」场景；C 假设所有位点共享同一效应分布
- --analysis-mode：有原始个体基因型用 individual；只有 GWAS 汇总统计用 summary（此时必须构建 LD 矩阵）
- --ld-matrix-type：默认 sparse（省内存）；样本/位点很多时仍超内存可试 block；eigen 用于后续预测
- --maf-threshold / --miss-threshold：默认 0.01 / 0.1，一般不用动
- --no-batch：表型文件只有一个性状或想手动控制时禁用批量
- --step：只想补跑某一步（如只重算 LD 矩阵）时指定 convert/qc/freq/ld/analysis

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --vcf-file` | 必填 |  | VCF变异文件｜VCF variant file |
| `-p, --pheno-file` | 必填 |  | 表型文件或GWAS汇总统计文件｜Phenotype file or GWAS summary statistics file |
| `-o, --output-dir` | `./gctb_output` |  | 输出目录｜Output directory |
| `--gctb-path` | `~/miniforge3/envs/gctb/bin/gctb` |  | GCTB软件路径｜GCTB software path |
| `--plink-path` | `~/miniforge3/envs/pop/bin/plink` |  | PLINK软件路径｜PLINK software path |
| `--maf-threshold` | `0.01` | float | MAF阈值｜MAF threshold |
| `--miss-threshold` | `0.1` | float | 缺失率阈值｜Missing rate threshold |
| `--bayes-type` | `S` | S/R/C | 贝叶斯模型类型｜Bayesian model type |
| `--analysis-mode` | `individual` | individual/summary | 分析模式｜Analysis mode |
| `--ld-matrix-type` | `sparse` | sparse/block/eigen | LD矩阵类型｜LD matrix type |
| `-t, --threads` | `12` | int | 线程数｜Number of threads |
| `--batch/--no-batch` | `True` |  | 批量处理多个表型（默认开启）｜Batch process multiple phenotypes (enabled by default) |
| `--step` | — | convert/qc/freq/ld/analysis | 只运行指定步骤｜Run only specified step |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --vcf-file` | 必填 |  | VCF变异文件路径｜VCF variant file path |
| `-p, --pheno-file` | 必填 |  | 表型文件路径（个体水平分析）或GWAS汇总统计文件（汇总水平分析）｜Phenotype file path (individual-level) or GWAS summary statistics file (summary-level) |
| `-o, --output-dir` | `./gctb_output` |  | 输出目录｜Output directory |
| `--gctb-path` | `~/miniforge3/envs/gctb/bin/gctb` |  | GCTB软件路径｜GCTB software path |
| `--plink-path` | `~/miniforge3/envs/pop/bin/plink` |  | PLINK软件路径｜PLINK software path |
| `--maf-threshold` | `0.01` | float | MAF阈值｜MAF threshold |
| `--miss-threshold` | `0.1` | float | 缺失率阈值｜Missing rate threshold |
| `--hwe-p` | `1e-06` | float | HWE p值阈值｜HWE p-value threshold |
| `--bayes-type` | `S` | S/R/C | 贝叶斯模型类型｜Bayesian model type (S/R/C) |
| `--analysis-mode` | `individual` | individual/summary | 分析模式｜Analysis mode |
| `--ld-matrix-type` | `sparse` | sparse/block/eigen | LD矩阵类型｜LD matrix type (sparse/block/eigen) |
| `--threads` | `12` | int | 线程数｜Number of threads |
| `--seed` | — | int | 随机种子｜Random seed |
| `--pi` | — | float | polygenicity参数｜Polygenicity parameter |
| `--sigma-g` | — | float | 遗传方差｜Genetic variance |
| `--rho` | — |  | SNP效应与MAF关系参数｜Parameter for effect-MAF relationship |
| `--step` | — | convert/qc/freq/ld/analysis | 只运行指定步骤｜Run only specified step |
| `--batch` | `True` | store_true | 批量处理多个表型（默认开启，使用--no-batch禁用）｜Batch process multiple phenotypes (enabled by default, use --no-batch to disable) |
| `--no-batch` | — | store_false | 禁用批量处理模式｜Disable batch processing mode |
| `--keep-intermediate` | — | store_true | 保留中间文件｜Keep intermediate files |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- GCTB 2.5+（默认 ~/miniforge3/envs/gctb/bin/gctb，可用 --gctb-path 指定）
- PLINK 1.9/2.0（默认 ~/miniforge3/envs/pop/bin/plink，可用 --plink-path 指定）
- Python 3 + pandas
- conda（通过 conda run 调用 gctb/plink）

## 常见问题 | FAQ

**Q1：换参数重跑，为什么某些步骤被跳过？**
本模块有断点续传，按各步骤输出文件是否存在判断。换 --maf-threshold / --miss-threshold / --bayes-type 等参数重跑旧目录前，先删对应旧产物（如 02_quality_control/input_qc.bed、04_gctb_analysis/ 下的 .snpRes），否则会复用旧参数的结果。

**Q2：汇总模式报「LD matrix not found」？**
summary 模式需要 LD 矩阵，但完整流程里只有 summary 模式才构建 LD 矩阵。确保用 --analysis-mode summary 时从 convert 步骤完整跑（或先用 --step ld 单独构建），并确认 --ld-matrix-type 与矩阵文件对应。

**Q3：批量模式没生成 summary/ 目录？**
summary/ 只在批量模式下所有性状都成功后才生成；若某性状失败，先看 04_gctb_analysis/ 下对应子目录的 .log 排查，失败过多时可用 --no-batch 逐个排查。

**Q4：内存不足？**
优先用默认 sparse LD 矩阵（最省内存），再考虑减少 --threads；个体水平大样本可试着降低位点数（提高 --maf-threshold）。
