# vcf2gwas GWAS 分析 | vcf2gwas GWAS Analysis

一句话理解：**一个轻量封装，把第三方工具 vcf2gwas 用 conda 环境调起来做全基因组关联分析(GWAS)，所有参数原样透传**——输入 VCF + 表型表，输出关联结果和曼哈顿图/QQ图。

## 功能概述 | Overview

- 通过 `conda run` 调用 vcf2gwas 的 conda 环境，所有命令行参数**直接透传**给原始工具
- vcf2gwas 提供线性混合模型(LMM)和广义线性模型(GLM)关联分析，内置群体结构校正和亲缘矩阵计算
- 适合中等到大规模群体的 GWAS 场景
- **本模块没有自己的分析逻辑**：输入/输出/参数全部由 vcf2gwas 决定

## 快速开始 | Quick Start

```bash
biopytools vcf2gwas -v input.vcf.gz -pf phenotype.csv -p 1 -lmm
```

最小输入：一个 VCF + 一个表型文件，用 LMM 模型跑关联分析。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| GWAS | 全基因组扫描，找出「哪些位点与性状显著相关」 |
| LMM / GLM | 两种统计模型：LMM 校正亲缘关系更稳，GLM 更快 |
| 表型列(-p) | 表型文件里第几列是要分析的性状 |
| 协变量(-c) | 需要「扣除」的干扰因素(如性别、群体结构) |
| PCA 主成分(-P) | 校正群体结构的成分数量 |
| 曼哈顿图 / QQ图 | 前者找「冒尖的塔」(显著位点)，后者诊断整体信号是否正常 |

## 输入 | Input

- `-v`：输入 VCF(通常 .vcf.gz)。
- `-pf`：表型文件(CSV/TSV，含样本名列和性状列)。
- `-cf`：协变量文件(可选)。
- 样本名须与 VCF 一致。具体列格式以 vcf2gwas 官方文档为准(见 FAQ)。

## 参数说明 | Parameters

### 包装器选项 | Wrapper options

**通俗理解|In plain words:** 只有一个参数——用哪个 conda 环境里的 vcf2gwas。**默认环境名一般不用动**，除非你自己建了别的环境。

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--vcf2gwas-env` | `vcf2gwas_v.0.8.9` | vcf2gwas conda 环境名 |

### vcf2gwas 透传参数 | Passthrough arguments

**通俗理解|In plain words:** 以下所有参数都**原样交给 vcf2gwas**，不是本包装器自己处理。最常用的几个：

| 参数 | 说明 |
|------|------|
| `-v` | 输入 VCF 文件(.vcf.gz) |
| `-pf` | 表型文件 |
| `-p` | 表型列号 |
| `-cf` | 协变量文件(可选) |
| `-c` | 协变量列号 |
| `-lmm` / `-glm` | 线性混合模型 / 广义线性模型 |
| `-P` | 主成分数量 |

> 完整参数列表请运行 `biopytools vcf2gwas --help`(会透传给 vcf2gwas 显示其原生帮助)。

## 输出 | Output

输出文件由 vcf2gwas 决定，典型包括：

- 关联分析结果(gwas_results 相关的 .txt/.csv)
- 曼哈顿图(manhattan_plot_*.png)
- QQ图(qq_plot_*.png)
- 运行日志

具体文件名以 vcf2gwas 实际产出为准。

## 结果解读 | Interpreting Results

**通俗理解|In plain words:** 关联结果表里一行一个位点，看 P 值——越小越显著；曼哈顿图里冒尖的塔就是候选位点；QQ 图右上角翘起说明有真信号。

## 参数选择建议 | Parameter Guidance

| 场景 | 建议 |
|------|------|
| 常规分析 | `-lmm`(校正亲缘关系，稳) |
| 快速初筛 | `-glm`(更快) |
| 有协变量 | 加 `-cf cov.csv -c 1` |
| 有群体结构 | 加 `-P 3` 用前 3 个主成分校正 |

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--vcf2gwas-env` | `vcf2gwas_v.0.8.9` |  | [STR] conda环境名｜conda env name (default: vcf2gwas_v.0.8.9) |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- vcf2gwas：conda 环境 `vcf2gwas_v.0.8.9`(https://github.com/frankraden/vcf2gwas)
- conda(推荐 miniforge)
- vcf2gwas 内部依赖(GCTA/LIMIX 等由其环境自带)

## 常见问题 | FAQ

**Q1：想看 vcf2gwas 的全部参数怎么办？**
运行 `biopytools vcf2gwas --help`，`--help` 会透传给 vcf2gwas 显示其原生帮助；不带任何参数则显示本包装器的简介。

**Q2：本工具会断点续传吗？**
不会。它只是把参数交给 vcf2gwas 跑一次，没有额外的续传逻辑(是否续传取决于 vcf2gwas 本身)。

**Q3：表型文件什么格式？**
以 vcf2gwas 官方为准，通常是 CSV/TSV，含样本名列和性状列，`-p` 指定用第几列做性状。

**Q4：换个 conda 环境怎么弄？**
`--vcf2gwas-env 你的环境名` 即可。
