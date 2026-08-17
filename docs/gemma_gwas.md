# GEMMA 全基因组关联分析 | GEMMA GWAS (Batch)

一句话理解：**把「群体里每个 DNA 位点跟某个性状（株高、抗病、产量等）是否相关」这件事，用 GEMMA 混合线性模型逐个位点打分**。输入一个基因型文件(VCF)和一份性状记录表(表型文件)，输出每个位点的关联打分和一张「显著位点统计表」。

## 功能概述 | Overview

- 端到端批量分析：一个 VCF + 一张多列表型文件，自动逐个表型跑完整 GWAS
- 内置样本交集筛选（bcftools），只保留 VCF 与表型里都有的样本
- 可选 PLINK 质控（MAF / 缺失率 / HWE），给了参数才启用，默认关闭
- 自动计算 PCA 主成分作为协变量，再用 GEMMA 计算亲缘关系矩阵(K)
- 用 GEMMA 的 LMM（线性混合模型）做关联检验，扣除群体结构与亲缘关系的干扰
- 汇总报告自动统计每个表型在不同显著性阈值下的显著位点数

## 快速开始 | Quick Start

```bash
biopytools gemma-gwas -i genotype.vcf.gz -p phenotype.txt
```

最小输入：一个 VCF（.vcf/.vcf.gz）+ 一个表型文件（TSV，首列为样本 ID 且带表头，其余列为各表型值）。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| SNP | 基因组上的一个「点位」，像高速公路上的里程桩，工具逐个检查这些桩和性状的关系 |
| 表型 | 你要研究的外在特征（株高、抗病、产量），像「成绩单」上的各科分数 |
| MAF | 「少数派等位基因」在群体里占多少；太低=这个位点几乎人人一样，区分不出差异 |
| 亲缘关系矩阵(K) | 样本之间「亲戚有多近」的表格，用来排除「长得像是因为亲戚近」的干扰 |
| PCA 协变量 | 群体结构的「主成分」，像把人群按祖先来源分组后作为背景变量扣除 |
| LMM | 混合线性模型，既看位点效应、又扣除亲缘+群体结构，抗假阳性更强 |
| P 值 | 该位点「与性状无关」的概率；越小越像真的有真实关联 |

## 输入 | Input

### VCF 文件

标准 VCF 格式（支持 .vcf / .vcf.gz），基因型编码 0/0、0/1、1/1。

### 表型文件

TSV 制表符分隔，**必须带表头**，第一列是样本 ID（须与 VCF 样本名一致），其余每一列是一个表型：

```text
sample    height    yield
21-18     12.5      3.2
21-19     15.2      4.1
```

## 分析流程 | Pipeline

```text
输入 VCF + 表型文件
    |
    v
步骤1: 准备表型文件(拆出带表头/不带表头两个版本)
    |
    v
步骤2: 样本交集筛选(bcftools 提取共同样本,转成 PLINK bed 格式)
    |
    v
步骤3: PLINK 质控(可选,仅当给了 --maf/--geno/--mind/--hwe 时才启用)
    |
    v
步骤4: PCA 分析(plink --pca,生成协变量文件)
    |
    v
步骤5: 计算亲缘关系矩阵(GEMMA -gk)
    |
    v
步骤6: 逐表型 GEMMA LMM 关联分析 + 汇总报告
```

## 输出 | Output

```text
gemma_results/
+-- genotype.bed/.bim/.fam         # VCF 转出的 PLINK 格式基因型
+-- pca.eigenvec                   # PCA 结果(前 n_pca 个主成分)
+-- covariate.txt                  # 协变量文件(供 GEMMA 使用)
+-- output/
|   +-- kinship.cXX.txt            # 亲缘关系矩阵
|   +-- 表型名_lmm.assoc.txt       # 每个表型的 GEMMA 关联结果
+-- gwas_summary.txt               # 汇总表:各表型显著位点数
+-- *.log                          # 各步骤日志
```

## 结果解读 | Interpreting Results

### 1. output/表型名_lmm.assoc.txt（主结果）

GEMMA 关联结果，一行一个 SNP。核心列：位点 ID、染色体(chr)、位置(ps)、以及按所选 LMM 方法输出的 P 值列（如 p_wald 对应 Wald 检验）。**P 值越小越显著**；常用 p<1e-5 记为「提示性显著」，p<5e-8 为「全基因组显著」（严格阈值需按位点数做校正）。

### 2. gwas_summary.txt（汇总表）

每行一个表型，列出总位点数与 p<1e-5 / p<1e-6 / p<1e-7 的显著位点数。**一眼看出哪个表型有强信号、信号强弱**。

## 参数选择建议 | Parameter Guidance

- --n-pca：默认 10，一般不用动；群体结构复杂（多个亚群混样）时可加到 20
- --maf / --geno / --mind / --hwe：给了任一参数才启用 PLINK 质控，**默认关闭**；数据较脏时建议给 --maf 0.05 --geno 0.1 先过滤一遍
- --lmm：默认 4（Wald/LRT/Score 全跑），一般不用动；只想看 Wald 检验用 1
- --gk：默认 1（中心化亲缘矩阵），一般不用动
- --notsnp：只要总效应、不关心每个 SNP 的估计值时加上，可提速
- --gemma：默认 ~/.local/bin/gemma，装了 conda 版 GEMMA 时用 --gemma 指定实际路径

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input, -i` | 必填 |  | VCF格式基因型文件｜VCF format genotype file |
| `--pheno, -p` | 必填 |  | 表型文件｜Phenotype file (first column: sample ID, has header) |
| `--output-dir, -o` | `gemma_results` | Path | 输出目录｜Output directory |
| `--n-pca` | `10` | int | PCA主成分数｜Number of PCA components |
| `--threads` | `12` | int | 线程数｜Number of threads |
| `--gemma` | `~/.local/bin/gemma` |  | GEMMA程序路径｜GEMMA program path |
| `--maf` | — | float | 最小等位基因频率阈值｜Minor allele frequency threshold |
| `--geno` | — | float | SNP缺失率阈值｜SNP missing rate threshold |
| `--mind` | — | float | 样本缺失率阈值｜Sample missing rate threshold |
| `--hwe` | — | float | Hardy-Weinberg p值阈值｜Hardy-Weinberg p-value threshold |
| `--no-qc` | — |  | 跳过PLINK质控｜Skip PLINK quality control |
| `--lmm` | — | int | LMM检验方法: 1=Wald, 2=LRT, 3=Score, 4=全部｜LMM test method: 1=Wald, 2=LRT, 3=Score, 4=all |
| `--gk` | — | int | 亲缘关系矩阵方法: 1=中心化, 2=标准化｜Kinship matrix method: 1=centered, 2=standardized |
| `--miss-gemma` | — | float | GEMMA缺失率阈值｜GEMMA missing rate threshold |
| `--maf-gemma` | — | float | GEMMA MAF阈值｜GEMMA MAF threshold |
| `--notsnp` | — |  | 不输出每个SNP的估计值(更快)｜Do not output estimated values for each SNP (faster) |
| `--verbose, -v` | — |  | 详细模式(-v: INFO, -vv: DEBUG)｜Verbose mode (-v: INFO, -vv: DEBUG) |
| `--quiet` | — |  | 静默模式(仅ERROR)｜Quiet mode (only ERROR) |
| `--log-file` | — | Path | 日志文件路径｜Log file path |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 输入VCF文件(可压缩)｜Input VCF file (can be compressed) |
| `-p, --pheno` | 必填 |  | 表型文件(第一列:样本ID,含表头)｜Phenotype file (first column: sample ID, has header) |
| `-o, --output-dir` | `gemma_results` |  | 输出目录｜Output directory |
| `--n-pca` | `10` | int | PCA主成分数｜Number of PCA components |
| `--threads` | `12` | int | 线程数｜Number of threads |
| `--gemma` | `~/.local/bin/gemma` |  | GEMMA程序路径｜GEMMA program path |
| `--maf` | — | float | 最小等位基因频率阈值｜Minor allele frequency threshold |
| `--geno` | — | float | SNP缺失率阈值｜SNP missing rate threshold |
| `--mind` | — | float | 样本缺失率阈值｜Sample missing rate threshold |
| `--hwe` | — | float | Hardy-Weinberg p值阈值｜Hardy-Weinberg p-value threshold |
| `--no-qc` | — | store_true | 跳过PLINK质控｜Skip PLINK quality control |
| `--lmm` | — | 1/2/3/4 | LMM检验方法｜LMM test method: 1=Wald, 2=LRT, 3=Score, 4=all |
| `--gk` | — | 1/2 | 亲缘关系矩阵方法｜Kinship matrix method: 1=centered, 2=standardized |
| `--miss-gemma` | — | float | GEMMA缺失率阈值｜GEMMA missing rate threshold |
| `--maf-gemma` | — | float | GEMMA MAF阈值｜GEMMA MAF threshold |
| `--notsnp` | — | store_true | 不输出每个SNP的估计值(更快)｜Do not output estimated values for each SNP (faster) |
| `-v, --verbose` | `0` | count | 详细模式(-v: INFO, -vv: DEBUG)｜Verbose mode (-v: INFO, -vv: DEBUG) |
| `--quiet` | — | store_true | 静默模式(仅ERROR)｜Quiet mode (only ERROR) |
| `--log-file` | — |  | 日志文件路径｜Log file path |
| `-V, --version` | — | version |  |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- GEMMA（默认 ~/.local/bin/gemma，可用 --gemma 指定）
- PLINK（VCF 转 bed、PCA、质控）
- bcftools（样本交集提取）
- awk（表型文件处理）
- Python 3

## 常见问题 | FAQ

**Q1：为什么默认不做质控？**
默认 plink_qc.enable=False，只有给了 --maf/--geno/--mind/--hwe 任一参数才启用 PLINK 质控。建议先用默认跑一遍，再按数据质量决定要不要加质控参数。

**Q2：GEMMA 报「找不到 output/kinship.cXX.txt」？**
GEMMA 会在输出目录下自动创建 output/ 子目录并写入亲缘矩阵。检查 --gemma 路径是否正确、GEMMA 是否有执行权限、前一步 gemma_kinship.log 里是否有报错。

**Q3：多个表型能一次跑吗？**
能。表型文件里放多列（除首列样本 ID 外的每一列都是一个表型），程序会逐个表型跑 GEMMA 并在 gwas_summary.txt 汇总。

**Q4：样本 ID 不一致会怎样？**
步骤 2 用 bcftools 求 VCF 与表型样本的交集，交集外的样本被剔除；若交集为空会直接报错退出。请核对样本名拼写与大小写。
