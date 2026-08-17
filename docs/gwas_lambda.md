# GWAS 基因组膨胀系数评估 | GWAS Lambda GC Calculator

一句话理解：**批量检查你的 GWAS 结果有没有「整体假阳性膨胀」，给每个结果打个质量分**。输入一堆 GWAS 结果文件，输出每个文件的 Lambda GC 值和质量评级。

## 功能概述 | Overview

- 批量计算多个 GWAS 结果文件的 Lambda GC 值，一次评估全部性状
- Lambda GC = 观测卡方中位数 / 期望卡方中位数，衡量群体分层/假阳性膨胀程度
- 同时统计每个文件的显著位点数(P < 阈值)与总位点数
- 按 Lambda 值自动分级：Ideal / Acceptable / Inflated / Deflated / No Signals
- 没有显著位点的文件不计算 Lambda（记为 NA），避免误导

## 快速开始 | Quick Start

```bash
biopytools gwas-lambda
```

默认按 feture_*/GWAS_Result.mlm.manht_input 模式搜索文件；也可用 -p 指定你自己的文件模式。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| Lambda GC | 基因组膨胀系数，衡量 P 值整体是否「虚高」，1.0 左右最理想 |
| 群体分层 | 样本里混了多个遗传背景不同的亚群，会制造大量假关联 |
| 卡方分布 | 把 P 值转成卡方统计量后的分布，是算 Lambda 的中间量 |
| 膨胀(Inflation) | Lambda > 1.1，像「天平没校准」，很多位点被误判显著 |

## 输入 | Input

通过 -p 指定 glob 模式匹配 GWAS 结果文件（制表符/空格分隔，P 值在指定列）：

```text
Chromosome    Position    Marker    Allele1    Allele2    P-value
1             14370       rs6054257 G          A          2.1e-08
1             17330       .         T          A          1.5e-06
```

P 值列默认第 4 列（0-based 索引 3），用 -c 指定其他列。程序自动跳过表头和非数值行。

## 输出 | Output

输出单个 TSV 文件（默认 Batch_Lambda_Assessment.txt）：

```text
Folder    Lambda_GC    Sig_SNPs(<1e-5)    Total_SNPs    Status
trait_1   1.0234       45                  1057153       Acceptable (可接受)
trait_2   1.2345       128                 1057153       Inflated (膨胀/假阳性高)
trait_3   NA           0                   1057153       No Signals (无显著关联)
```

## 结果解读 | Interpreting Results

| Lambda GC 范围 | 状态 | 怎么理解 |
|------|------|----------|
| 0.95 – 1.05 | Ideal（完美） | 无群体分层，结果可靠 |
| 0.90 – 1.10 | Acceptable（可接受） | 轻微偏差，基本可信 |
| > 1.10 | Inflated（膨胀） | 有群体分层或假阳性，需 MLM/PCA 校正 |
| 0.80 – 0.90 | Deflated（过度校正） | 校正过头，检查协变量设置 |
| < 0.80 | Extreme Deflated | 模型可能失效 |
| NA | No Signals | 无 P<阈值 位点，正常，无需关注 |

**重点看 Inflated 的文件**：Lambda 明显 >1.1 说明该性状的 GWAS 假阳性偏高，需要回炉（加 PCA/亲缘关系校正）。

## 参数选择建议 | Parameter Guidance

- -p / --pattern：默认匹配 TASSEL 的 manht_input 输出；改成你自己的结果文件路径模式即可
- -c / --p-column：0-based 列索引，P 值不在第 4 列时才需要动
- -t / --threshold：默认 1e-5，只影响「显著位点数」统计，不影响 Lambda 计算
- -o / -d：输出文件名与目录，一般不用动

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--pattern, -p` | `feture_*/GWAS_Result.mlm.manht_input` |  | 文件搜索模式｜File search pattern |
| `--output, -o` | `Batch_Lambda_Assessment.txt` |  | 输出文件名｜Output filename |
| `--threshold, -t` | `1e-05` | float | 显著性阈值｜Significance threshold |
| `--p-column, -c` | `3` | int | P值列索引｜P-value column index (0-based) |
| `--output-dir, -d` | `./gwas_lambda_output` | Path | 输出目录｜Output directory |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--pattern, -p` | `feture_*/GWAS_Result.mlm.manht_input` |  | 文件搜索模式｜File search pattern |
| `--output, -o` | `Batch_Lambda_Assessment.txt` |  | 输出文件名｜Output filename |
| `--threshold, -t` | `1e-05` | float | 显著性阈值｜Significance threshold |
| `--p-column, -c` | `3` | int | P值所在列索引｜P-value column index (0-based) |
| `--output-dir, -d` | `./gwas_lambda_output` |  | 输出目录｜Output directory |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- Python 3 + numpy、scipy（卡方分布 ppf 计算）

## 常见问题 | FAQ

**Q1：为什么有些文件的 Lambda 显示 NA？**
该文件里没有 P < 阈值（默认 1e-5）的位点，程序不计算 Lambda（避免只有少数点导致失真），记为 No Signals。这是正常情况。

**Q2：Lambda 明显大于 1 说明什么？**
说明该 GWAS 结果存在群体分层或假阳性膨胀（P 值整体偏小）。建议回炉：用 MLM 混合模型、加 PCA 协变量、或做基因组控制校正。

**Q3：P 值不在第 4 列怎么办？**
用 -c 指定 0-based 列索引（第 1 列=0，第 4 列=3）。例如 P 值在第 6 列就 -c 5。

**Q4：一个文件都匹配不到？**
检查 -p 的 glob 模式是否写对（相对当前目录），或用绝对路径；模式要能匹配到结果文件本体。
