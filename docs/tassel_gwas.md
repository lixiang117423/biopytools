# TASSEL 全基因组关联分析 | TASSEL GWAS

一句话理解：**用 TASSEL 软件把「群体里每个位点和多个性状是否相关」这件事，一次批量全跑完**。输入一个 VCF 和一张多列表型文件，自动为每个性状分别做 GWAS，输出每个性状的关联结果和一张批量汇总报告。

## 功能概述 | Overview

- 自动识别表型文件里的所有性状列，逐个性状独立分析，也可 --parallel 并行提速
- 支持 GLM（一般线性模型）、MLM（混合线性模型，默认）、BOTH 三种模式
- MLM 模式自动预计算亲缘关系矩阵(K)与 PCA 协变量，所有性状共享，避免重复计算
- 可复用外部 Q 矩阵（群体结构）和 K 矩阵（亲缘关系）
- 每个性状产出独立的曼哈顿图输入文件(.manht_input)，批量报告汇总成功/失败

## 快速开始 | Quick Start

```bash
biopytools tassel-gwas -i input.vcf.gz -p traits.txt -o results
```

最小输入：一个 VCF + 一个 TASSEL 格式的表型文件（见下）。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| GLM | 一般线性模型，只做位点与性状的回归，快但容易受群体结构干扰、假阳性多 |
| MLM | 混合线性模型，额外扣除「亲戚关系」和「群体分群」的影响，结果更可靠（默认） |
| 亲缘关系矩阵(K) | 样本两两「亲戚有多近」的表，用来排除「长得像是因为亲戚近」 |
| Q 矩阵 | 群体结构的量化表（常由 PCA 或 STRUCTURE 得到），作为背景变量扣除 |
| PCA 协变量 | 主成分，把人群按遗传背景排序，前几列作为协变量 |
| 曼哈顿图输入(.manht_input) | 一行一个位点 + P 值，供 R 等工具画「山峰图」 |

## 输入 | Input

### VCF 文件

标准 VCF（.vcf / .vcf.gz）。

### 表型文件（TASSEL 格式）

制表符分隔，第一行表头、第二行数据类型标记，第一列是样本 ID，其余每列一个性状：

```text
<Trait>    Trait1    Trait2
<Data>     <Data>    <Data>
Sample1    12.5      0.85
Sample2    15.2      0.92
```

## 分析流程 | Pipeline

```text
输入 VCF + 表型文件
    |
    v
自动识别所有性状列
    |
    v
[MLM 模式] 预计算共享矩阵: 排序/过滤 VCF → 亲缘关系矩阵(K) + PCA 协变量
    |
    v
逐个性状(可并行):
    提取单性状表型文件 → 运行 GLM 和/或 MLM → 提取结果 → 统计报告
    |
    v
批量汇总报告(batch_processing_report.txt)
```

## 输出 | Output

```text
results/
+-- batch_processing_report.txt     # 批量汇总:各性状成功/失败、成功率
+-- failed_traits.log               # 失败性状及原因
+-- gwas.log                        # 全局日志
+-- 性状1/                          # 每个性状一个子目录
|   +-- 性状1_GWAS.mlm.manht_input  # 曼哈顿图输入(4列:Chr Pos SNP P-value)
|   +-- 性状1_GWAS.glm.manht_input  # (GLM 模式时)
|   +-- 性状1_GWAS.pipeline.log     # 该性状运行日志
|   +-- 性状1.pheno.txt             # 提取出的单性状表型
|   +-- 性状1_GWAS.stats.txt        # 该性状统计报告
+-- 性状2/
+-- ...
```

## 结果解读 | Interpreting Results

### 1. 性状_GWAS.mlm.manht_input（主结果）

四列：Chr（染色体）、Pos（位置）、SNP（位点名）、P-value（P 值）。**P 值越小越显著**。把这个文件喂给 R 的 qqman/CMplot 就能画曼哈顿图和 QQ 图；关注「山峰」最高的区域，那是最可能的候选区间。

### 2. batch_processing_report.txt（批量汇总）

列出总性状数、成功数、失败数、成功率，以及每个成功性状的输出目录。**先看这里定位哪些性状跑成功、哪些失败了**。

### 3. 性状_GWAS.stats.txt（统计报告）

记录该性状的变异位点数、样本数、模型、显著位点数(P<1e-5)、运行时长等，用于核对分析规模是否正常。

## 参数选择建议 | Parameter Guidance

- --model：默认 MLM，有群体结构时首选；追求速度或群体同质时用 GLM；想对比用 BOTH
- --memory：默认 100g（Java 堆内存），大数据集可加大到 200g~400g，内存不足报错时先调这里
- --pca-components：默认 5（MLM 协变量数量），一般不用动，群体结构复杂可加到 10
- --maf / --miss：给了才在 VCF 排序过滤时启用，默认不额外过滤
- --q-matrix / --kinship：已有现成的群体结构/亲缘矩阵时传入，可跳过内部计算
- --parallel --workers：多性状时开启并行，workers 默认 4，按机器核数调整
- --skip-sort：默认已跳过 VCF 排序（假设输入已过滤好）；输入是原始未过滤 VCF 时去掉它

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--vcf, -i` | 必填 |  | VCF基因型文件｜VCF genotype file |
| `--pheno, -p` | 必填 |  | 表型数据文件｜Phenotype data file |
| `--output, -o` | 必填 |  | 输出目录｜Output directory |
| `--model, -m` | `MLM` | GLM/MLM/BOTH | GWAS模型｜GWAS model (GLM/MLM/BOTH) |
| `--memory` | `100g` |  | Java最大内存｜Java maximum memory (e.g., 100g, 200g) |
| `--threads, -t` | `12` | int | 并行线程数｜Number of parallel threads |
| `--maf` | — | float | 最小等位基因频率过滤｜Minimum allele frequency filter |
| `--miss` | — | float | 最大缺失率过滤｜Maximum missing rate filter |
| `--pca-components` | `5` | int | PCA主成分数(MLM协变量)｜Number of PCA components for MLM covariates |
| `--q-matrix, -q` | — | Path | 群体结构Q矩阵文件｜Population structure Q matrix file |
| `--kinship, -k` | — | Path | 亲缘关系矩阵文件｜Kinship matrix file |
| `--tassel-path` | — | Path | TASSEL安装路径｜TASSEL installation path |
| `--skip-sort` | — |  | 跳过VCF排序｜Skip VCF sorting |
| `--keep-temp` | — |  | 保留临时文件｜Keep temporary files |
| `--parallel` | — |  | 并行处理多个性状｜Parallel process multiple traits |
| `--workers` | `4` | int | 并行工作进程数｜Number of parallel workers |
| `--log-level` | `INFO` | DEBUG/INFO/WARN/ERROR | 日志级别｜Log level |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --vcf` | 必填 |  | VCF基因型文件路径｜VCF genotype file path |
| `-p, --pheno` | 必填 |  | 表型文件路径｜Phenotype file path (可包含多个表型) |
| `-o, --output` | 必填 |  | 输出目录路径｜Output directory path |
| `--model` | `MLM` | GLM/MLM/BOTH | GWAS模型选择｜GWAS model selection |
| `--memory` | `100g` |  | Java最大内存｜Java maximum memory |
| `--threads` | `4` | int | 并行线程数｜Number of parallel threads |
| `--maf` | — | float | 最小等位基因频率过滤｜Minimum allele frequency filter |
| `--miss` | — | float | 最大缺失率过滤｜Maximum missing rate filter |
| `--pca-components` | `5` | int | PCA主成分数量（用作MLM协变量）｜ Number of PCA components for MLM covariates |
| `--q-matrix` | — |  | 群体结构Q矩阵文件｜Population structure Q matrix file |
| `--kinship` | — |  | 亲缘关系K矩阵文件｜Kinship K matrix file |
| `--tassel-path` | — |  | TASSEL安装路径｜TASSEL installation path |
| `--no-skip-sort` | — | store_false | 不跳过VCF排序（默认跳过，适用于已过滤的VCF）｜Do not skip VCF sorting (default: skip, suitable for pre-filtered VCF) |
| `--keep-temp` | — | store_true | 保留临时文件｜Keep temporary files |
| `--parallel` | — | store_true | 并行处理多个表型｜Parallel process multiple traits |
| `--workers` | `4` | int | 并行工作线程数｜Number of parallel workers |
| `--log-level` | `INFO` | DEBUG/INFO/WARN/ERROR | 日志级别｜Log level |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- TASSEL 5.x（perl run_pipeline.pl 调用，可用 --tassel-path 指定）
- Java 8+（TASSEL 运行环境）
- bcftools（VCF 排序/过滤时使用，跳过排序则不需要）
- VCF2PCACluster（MLM 模式自动算 PCA 时使用；提供 --q-matrix 则可不用）
- Python 3

## 常见问题 | FAQ

**Q1：表型文件为什么要 <Trait> 和 <Data> 两行表头？**
这是 TASSEL 的标准输入格式：第一行是性状名、第二行声明数据类型（<Data> 表示数值）。第一列必须是样本 ID。少了这两行 TASSEL 会解析失败。

**Q2：换参数重跑，为什么结果没变？**
批量处理没有断点续传，但会复用已生成的中间矩阵/文件。换 --maf / --miss / --model 等关键参数重跑旧目录前，先删掉旧输出目录（或换新 -o），否则可能读到旧矩阵。

**Q3：MLM 报「未提供 K 矩阵但跳过了排序」？**
跳过排序(--skip-sort)时程序假设 K 矩阵已预计算。要么去掉 --skip-sort 让程序自动算，要么用 --kinship 传入现成矩阵。

**Q4：PCA 计算失败？**
MLM 自动算 PCA 需要 VCF2PCACluster（默认找 ~/software/VCF2PCACluster-1.42/bin/VCF2PCACluster）。找不到时用 --q-matrix 传入外部 PCA 文件即可绕过。

**Q5：某些性状失败怎么排查？**
先看 failed_traits.log 定位失败性状，再看对应子目录里的 性状_GWAS.pipeline.log 查具体报错（常见是内存不足或样本不匹配）。
