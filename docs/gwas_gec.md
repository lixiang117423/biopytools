# GWAS 多重检验校正 | GWAS Genome-wide Error Correction (GEC)

一句话理解：**算出你这次 GWAS「真正独立测了几次」，从而给出一个更准确、不那么保守的显著性阈值**。输入 GWAS 的 P 值文件和参考 VCF，输出校正后的 P 值阈值。

## 功能概述 | Overview

- 用 GEC 算法（KGGSee 实现）基于参考群体的 LD 结构，估算「有效独立检验数」
- 有效独立检验数显著小于实际 SNP 数（相邻位点往往绑定遗传，不算独立检验）
- 计算校正后阈值 = alpha / 总有效检验数，比传统 Bonferroni 更宽松、更准确
- 汇总报告含各染色体的 LD 块数与有效检验数
- 自动把 Chr01 这类染色体名转成数字格式（KGGSee 需要）

## 快速开始 | Quick Start

```bash
biopytools gwas-gec -i gwas.txt -r input.vcf.gz
```

最小输入：一个 GWAS P 值文件（含染色体、位置、P 值三列）+ 一个参考 VCF。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| 多重检验 | 测几百万个位点，纯靠运气也会有几千个「碰巧显著」，必须校正 |
| LD 块 | 相邻位点往往「绑定遗传」，一个块里的位点信息高度重复，只算一次独立检验 |
| 有效独立检验数 | 去掉重复后「实际独立的检验次数」，是校正阈值的分母 |
| FWER | 家族错误率，控制「全基因组至少出现一个假阳性」的概率 |

## 输入 | Input

### GWAS P 值文件

制表符分隔、带表头，需含染色体、位置、P 值三列（列名可用 --chrom-col/--pos-col/--p-col 指定，默认 CHR/BP/P）：

```text
CHR    BP         P
1      9572       0.9799
1      752566     1.2e-05
```

### 参考 VCF

研究群体的基因型（.vcf / .vcf.gz），用于计算 LD 结构。KGGSee 仅支持 VCF 格式。

## 分析流程 | Pipeline

```text
输入 P 值文件 + 参考 VCF
    |
    v
1. 自动转换染色体格式(Chr01 → 1,可选)
    |
    v
2. 运行 GEC 算法(KGGSee --var-gec,基于 LD 计算有效检验数)
    |
    v
3. 解析结果,计算校正后阈值 = alpha / 总有效检验数
    |
    v
4. 生成汇总报告
```

## 输出 | Output

```text
gec_output/
+-- gwas_gec.effective.size.txt.gz   # 每个 LD 块的有效检验数
+-- gwas_gec_summary.txt             # 汇总报告(总块数/总有效检验数/校正阈值/各染色体统计)
+-- gwas_gec.log                     # 运行日志
```

## 结果解读 | Interpreting Results

### 1. gwas_gec_summary.txt（核心结论）

重点看两行：**总有效检验数**和**校正后阈值**。校正阈值 = 0.05 / 总有效检验数。例如有效检验数 50 万，阈值约 1e-7——用这个阈值去判你 GWAS 结果里的显著位点，比一刀切的 5e-8 更贴合你的数据。

### 2. gwas_gec.effective.size.txt.gz（明细）

每个 LD 块的 Chrom（染色体）、StartPos/EndPos（起止）、Num（块内 SNP 数）、EffectiveNum（有效检验数）。有效检验数越接近 Num 说明该区域 LD 越弱、位点越独立。

## 参数选择建议 | Parameter Guidance

- -t / --threads：默认 12，一般不用动
- -m / --memory：默认 100G（Java 内存），VCF 很大时按需加大；内存不足报错先调这里
- --maf-filter / --p-cutoff：默认 0.05，一般不用动
- --chrom-col / --pos-col / --p-col：P 值文件列名不是 CHR/BP/P 时用它们指定
- --alpha：默认 0.05（FWER），想要更严格的阈值可调小
- --no-convert-chrom：确定染色体格式已是数字时可加，跳过自动转换
- --kggsee-jar：默认 ~/software/kmmsee/kggsee.jar，若未装在此路径需指定

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--pfile, -i` | 必填 |  | GWAS P值汇总统计文件｜GWAS P-value summary statistics file |
| `--reference, -r` | 必填 |  | 参考文件（VCF或PLINK binary前缀）｜Reference file (VCF or PLINK binary prefix) |
| `--output-dir, -o` | `./gec_output` | Path | 输出目录｜Output directory |
| `--threads, -t` | `12` | int | 线程数｜Number of threads |
| `--memory, -m` | `100G` |  | Java内存分配｜Java memory allocation |
| `--maf-filter` | `0.05` | float | MAF过滤阈值｜MAF filter threshold |
| `--p-cutoff` | `0.05` | float | P值阈值｜P-value threshold |
| `--no-keep-ref` | — |  | 不保留参考文件缓存｜Do not keep reference file cache |
| `--no-convert-chrom` | — |  | 禁用自动染色体格式转换｜Disable automatic chromosome format conversion |
| `--chrom` | — |  | 染色体范围(如1-22或1-10,22)｜Chromosome range (e.g., 1-22 or 1-10,22) |
| `--chrom-col` | `CHR` |  | 染色体列名｜Chromosome column name |
| `--pos-col` | `BP` |  | 位置列名｜Position column name |
| `--p-col` | `P` |  | P值列名｜P-value column name |
| `--alpha` | `0.05` | float | 显著性水平(FWER)｜Significance level |
| `--kggsee-jar` | `~/software/kmmsee/kggsee.jar` |  | KGGSee JAR文件路径｜KGGSee JAR file path |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-p, --pfile` | 必填 |  | GWAS P值汇总统计文件｜GWAS P-value summary statistics file |
| `-r, --reference` | 必填 |  | 参考文件（VCF文件或PLINK binary前缀）｜Reference file (VCF or PLINK binary prefix) |
| `-o, --output-dir` | `./gec_output` |  | 输出目录｜Output directory |
| `-t, --threads` | `16` | int | 线程数｜Number of threads |
| `-m, --memory` | `8g` |  | Java内存分配｜Java memory allocation |
| `--maf-filter` | `0.05` | float | MAF过滤阈值｜MAF filter threshold |
| `--p-cutoff` | `0.05` | float | P值阈值｜P-value threshold |
| `--no-keep-ref` | — | store_false | 不保留参考文件缓存｜Do not keep reference file cache |
| `--chrom` | — |  | 染色体范围｜Chromosome range |
| `--no-convert-chrom` | — | store_false | 禁用自动染色体格式转换｜Disable automatic chromosome format conversion |
| `--chrom-col` | `CHR` |  | 染色体列名｜Chromosome column name |
| `--pos-col` | `BP` |  | 位置列名｜Position column name |
| `--p-col` | `P` |  | P值列名｜P-value column name |
| `--alpha` | `0.05` | float | 显著性水平(FWER)｜Significance level (FWER) |
| `--kggsee-jar` | `~/software/kmmsee/kggsee.jar` |  | KGGSee JAR文件路径｜KGGSee JAR file path |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- Java（运行 KGGSee）
- KGGSee JAR（默认 ~/software/kmmsee/kggsee.jar，可用 --kggsee-jar 指定）
- Python 3

## 常见问题 | FAQ

**Q1：报「KGGSee JAR 文件不存在」？**
默认找 ~/software/kmmsee/kggsee.jar。用 --kggsee-jar 指定你机器上 kggsee.jar 的实际路径。

**Q2：参考文件能用 PLINK bed 格式吗？**
不能。KGGSee 只支持 VCF 格式的参考文件，请提供 .vcf 或 .vcf.gz。

**Q3：P 值文件染色体写 Chr01、VCF 写 1，会不匹配吗？**
程序默认会自动把 Chr01→1 转换，无需手动改。若不想转换用 --no-convert-chrom。

**Q4：校正阈值比 5e-8 大，能用吗？**
能。GEC 阈值基于你数据的真实 LD 结构计算，比固定的 5e-8 更准确，通常更宽松，能发现更多真实关联。
