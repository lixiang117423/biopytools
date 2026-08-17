# PLINK 全基因组关联分析 | PLINK GWAS

一句话理解：**用 PLINK 把「群体里每个 DNA 位点和性状是否相关」这件事，从头到尾自动跑一遍**。输入 VCF 和表型文件，自动完成质控、群体结构校正、关联分析和多重检验校正，输出显著位点清单和曼哈顿图/QQ 图。

## 功能概述 | Overview

- 端到端流程：质控 → LD 剪枝 → PCA → 关联分析（含 PCA 校正）→ 多重检验校正 → 汇总与绘图
- 支持质量性状（病例/对照，逻辑回归）与数量性状（连续值，线性回归）两种表型
- 支持加性(additive)、显性(dominant)、隐性(recessive)、全部(all)四种遗传模型
- 三种显著性校正任选或全跑：Bonferroni、提示性阈值、FDR（Benjamini-Hochberg）
- 自动生成曼哈顿图与 QQ 图，多模型时附带模型比较报告

## 快速开始 | Quick Start

```bash
biopytools plink-gwas -i data.vcf.gz -p pheno.txt -o results
```

最小输入：一个 VCF（.vcf/.vcf.gz）+ 一个两列表型文件（样本 ID + 表型值）。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| MAF | 「少数派等位基因」占比，太低=这个位点几乎人人一样，没信息量 |
| 缺失率 | 多少样本在该位点「没测出来」，缺失太多=判断依据不足 |
| HWE | Hardy-Weinberg 平衡检验，像「抽检位点是否符合自然规律」，严重偏离的位点多为测序错误 |
| LD 剪枝 | 相邻位点往往「绑定遗传」、信息重复，剪枝=每个团只留一个代表 |
| PCA | 主成分分析，把人群按遗传背景分组，作为协变量扣除「人群本身分群」带来的假信号 |
| 遗传模型 | 加性=每多一个等位基因效果累加；显性=带一个就够；隐性=要两个都带才显效 |
| Bonferroni / FDR | 多重检验校正，用来回答「这个显著是真的还是运气」（测了几百万次总有几个会碰巧显著） |
| OR | 比值比，质量性状里 >1 表示该等位基因增加患病风险，<1 表示保护 |
| BETA | 效应方向，数量性状里 >0 增大表型值，<0 减小 |

## 输入 | Input

### VCF 文件

标准 VCF（.vcf / .vcf.gz）。程序内部转成 PLINK bed 格式，支持非常规染色体名（如 OV、scaffold）。

### 表型文件

两列，第一列样本 ID、第二列表型值（分隔符自动识别）：

质量性状（默认，0=对照/抗病、1=病例/感病；内部转成 PLINK 的 1/2 编码，缺失记为 -9）：

```text
sample    disease
S1        0
S2        1
S3        0
```

数量性状（用 -T quantitative 指定，数值直接使用）：

```text
sample    height
S1        12.5
S2        15.2
```

## 分析流程 | Pipeline

```text
输入 VCF + 表型文件
    |
    v
1. 复制输入文件到输出目录
    |
    v
2. 转换表型(质量性状 0/1 转 1/2,缺失转 -9)
    |
    v
3. VCF 转 PLINK bed 格式(raw_data)
    |
    v
4. 合并表型(data_with_pheno)
    |
    v
5. 质量控制(--mind/--geno/--maf/--hwe → data_qc1)
    |
    v
6. 群体结构分析(LD 剪枝 + PCA)
    |
    v
7. 关联分析(--logistic 或 --linear,基础 + PCA 校正两套)
    |
    v
8. 结果处理(显著性校正)
    |
    v
9. 生成汇总报告   →   10. 生成曼哈顿图/QQ图
```

## 输出 | Output

所有文件平铺在输出目录（默认 gwas_results/）下：

```text
gwas_results/
+-- gwas_results_ADD.txt            # 主结果(加性模型;显性/隐性为 DOM/REC)
+-- significant_hits_bonferroni_additive.txt   # Bonferroni 显著位点
+-- suggestive_hits_additive.txt               # 提示性显著位点
+-- gwas_results_fdr_additive.txt              # FDR 校正后完整结果(含 q 值)
+-- fdr_significant_hits_additive.txt          # FDR 显著位点
+-- gwas_summary_report.txt         # 汇总报告(样本数/SNP数/各校正显著数)
+-- model_comparison_report.txt     # 模型比较报告(仅 -m all)
+-- manhattan_plot_additive.png     # 曼哈顿图
+-- qq_plot_additive.png            # QQ 图
+-- plink_analysis.log              # 运行日志
+-- data_qc1.bed/.bim/.fam          # 质控后 PLINK 数据(中间产物)
+-- pca_results.eigenvec            # PCA 结果
+-- gwas_basic_*.assoc.* / gwas_adjusted_*.assoc.*  # PLINK 原始关联结果
```

## 结果解读 | Interpreting Results

### 1. gwas_results_ADD.txt（主结果）

一行一个 SNP，关键列：CHR（染色体）、SNP（位点名）、BP（位置）、A1（效应等位基因）、OR/BETA（效应）、P（P 值）。**P 值越小越显著**，质量性状看 OR 判断风险方向，数量性状看 BETA 判断效应方向。

### 2. 校正后显著位点文件

- significant_hits_bonferroni_*.txt：按 0.05/总SNP数 的严格阈值筛出的位点，几乎无假阳性
- suggestive_hits_*.txt：P < 1e-5 的提示性位点，更宽松、适合初步挖掘
- fdr_significant_hits_*.txt：控制假发现率的位点，含 q 值列

### 3. 曼哈顿图与 QQ 图

曼哈顿图横轴染色体、纵轴 -log10(P)，越过虚线的「山峰」即候选位点；QQ 图观察点是否偏离对角线，**明显上扬=群体分层或假阳性**。

## 参数选择建议 | Parameter Guidance

- -T / -m：按表型性质选 quantitative/qualitative、按假设选遗传模型；不确定就 -m all 一次看全
- --maf / --geno：默认 0.05，一般不用动；--hwe / --mind 默认关闭，质量性状建议给 --hwe 1e-6
- --no-strat-corr：只有确定群体高度同质（无分层）时才加，可跳过 LD 剪枝和 PCA 提速；否则保持默认（做校正）
- --correction-method：默认 all（三种都跑）；只想要最严格的用 bonferroni，想多挖掘候选用 fdr
- --pca-use：默认用前 5 个主成分，一般不用动
- --force / --dry-run：换参数重跑旧目录前，用 --force 覆盖，或用 --dry-run 先预览

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--vcf, -i` | 必填 |  | VCF文件路径｜Path to VCF file |
| `--phenotype, -p` | 必填 |  | 表型文件路径｜Path to phenotype file |
| `--trait-type, -T` | `qualitative` | qualitative/quantitative | 表型类型｜Trait type |
| `--genetic-model, -m` | `additive` | additive/dominant/recessive/all | 遗传模型｜Genetic model |
| `--output-dir, -o` | `gwas_results` |  | 输出目录｜Output directory |
| `--no-strat-corr` | — |  | 禁用群体结构校正｜Disable population stratification correction |
| `--mind` | — | float | 个体缺失率阈值｜Individual missing rate threshold |
| `--geno` | `0.05` | float | SNP缺失率阈值｜SNP missing rate threshold |
| `--maf` | `0.05` | float | 最小等位基因频率｜Minor allele frequency |
| `--hwe` | — | float | Hardy-Weinberg平衡P值阈值｜HWE p-value threshold |
| `--ld-window-size` | `50` | int | LD窗口大小｜LD window size in kb |
| `--ld-step-size` | `5` | int | LD步长大小｜LD step size in SNPs |
| `--ld-r2-threshold` | `0.2` | float | LD r²阈值｜LD r² threshold |
| `--pca-components` | `10` | int | 主成分数量｜Number of PCA components |
| `--pca-use` | `5` | int | 关联分析中使用的主成分数量｜Number of PCs to use in association |
| `--correction-method` | `all` | bonferroni/suggestive/fdr/all | 显著性校正方法｜Significance correction method |
| `--bonferroni-alpha` | `0.05` | float | Bonferroni校正alpha水平｜Bonferroni alpha level |
| `--suggestive-threshold` | `1e-05` | float | 提示性关联阈值｜Suggestive threshold |
| `--fdr-alpha` | `0.05` | float | FDR校正q值阈值｜FDR q-value threshold |
| `--threads, -t` | `12` | int | 线程数｜Number of threads |
| `--verbose, -v` | — |  | 详细输出模式｜Verbose mode (-v: INFO, -vv: DEBUG) |
| `--quiet` | — |  | 静默模式｜Quiet mode (only ERROR) |
| `--log-file` | — | Path | 日志文件路径｜Log file path |
| `--force, -f` | — |  | 强制覆盖已存在的输出目录｜Force overwrite existing output directory |
| `--dry-run` | — |  | 模拟运行不实际执行分析｜Dry run without actual analysis |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --vcf` | 必填 |  | VCF文件路径｜Path to VCF file |
| `-p, --phenotype` | 必填 |  | 表型文件路径｜Path to phenotype file |
| `-T, --trait-type` | `qualitative` | qualitative/quantitative | 表型类型｜Trait type |
| `-m, --genetic-model` | `additive` | additive/dominant/recessive/all | 遗传模型｜Genetic model |
| `-o, --output-dir` | `gwas_results` |  | 输出目录｜Output directory |
| `--no-strat-corr` | — | store_false | 禁用群体结构校正（跳过LD剪枝和PCA）｜ Disable population stratification correction (skip LD pruning and PCA) |
| `--mind` | — | float | 个体缺失率阈值｜Individual missing rate threshold |
| `--geno` | `0.05` | float | SNP缺失率阈值｜SNP missing rate threshold |
| `--maf` | `0.05` | float | 最小等位基因频率｜Minor allele frequency |
| `--hwe` | — | float | Hardy-Weinberg平衡P值阈值｜HWE p-value threshold |
| `--ld-window-size` | `50` | int | LD剪枝窗口大小(kb)｜LD window size in kb |
| `--ld-step-size` | `5` | int | LD剪枝步长(SNP数)｜LD step size in SNPs |
| `--ld-r2-threshold` | `0.2` | float | LD剪枝r²阈值｜LD r² threshold |
| `--pca-components` | `10` | int | 计算的主成分数量｜Number of PCA components |
| `--pca-use` | `5` | int | 关联分析中使用的主成分数量｜Number of PCs to use in association |
| `--correction-method` | `all` | bonferroni/suggestive/fdr/all | 显著性校正方法｜Significance correction method |
| `--bonferroni-alpha` | `0.05` | float | Bonferroni校正alpha水平｜Bonferroni alpha level |
| `--suggestive-threshold` | `1e-05` | float | 提示性关联阈值｜Suggestive threshold |
| `--fdr-alpha` | `0.05` | float | FDR校正q值阈值｜FDR q-value threshold |
| `-t, --threads` | `1` | int | 使用的线程数｜Number of threads |
| `-v, --verbose` | `0` | count | 详细输出模式｜Verbose mode (-v: INFO, -vv: DEBUG) |
| `--quiet` | — | store_true | 静默模式，只输出ERROR｜Quiet mode (only ERROR) |
| `--log-file` | — |  | 日志文件路径｜Log file path |
| `-f, --force` | — | store_true | 强制覆盖已存在的输出目录｜Force overwrite existing output directory |
| `--dry-run` | — | store_true | 模拟运行，不实际执行分析｜Dry run without actual analysis |
| `-V, --version` | — | version |  |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- PLINK 1.9 / 2.0（GWAS 计算）
- Python 3 + pandas、numpy、matplotlib（绘图）、statsmodels（FDR 校正用 multipletests）

## 常见问题 | FAQ

**Q1：质量性状的 0/1 会被自动转换吗？**
会。内部把 0 转成 1（对照）、1 转成 2（病例），缺失值转成 -9。请确保表型文件里只有 0 和 1 两种取值，其他值会被当成缺失。

**Q2：为什么换参数重跑结果没变？**
本模块没有断点续传，但会直接覆盖输出目录里的同名文件；若担心旧文件干扰，用 --force 强制覆盖，或换一个新的 -o 目录。

**Q3：FDR 校正没产出结果？**
FDR 依赖 statsmodels 库。若环境里没装，会打 warning 并跳过 FDR（Bonferroni 和 suggestive 仍正常）。安装：pip install statsmodels。

**Q4：曼哈顿图没生成？**
绘图依赖 matplotlib 且要求结果文件里有有效的 CHR/BP/P 列。若 SNP 名不含位置信息且 PLINK 结果缺 BP 列，位置会被填 0，图可能异常，请检查 VCF 是否带标准位置信息。
