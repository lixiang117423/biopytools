# ATOMM 双物种关联分析 | ATOMM Two-Organism Mixed Model Association

一句话理解：**当「宿主感染病原后病得重不重」同时由宿主和病原双方的基因决定时，这个工具把两边的基因组放在一起扫，找出「宿主的哪个 SNP、病原的哪个 SNP、以及它们俩怎么配合」共同影响了感染结果**。输入宿主 VCF、病原 VCF 和一张「宿主 × 病原」交叉感染表型表，输出宿主边际效应、病原边际效应、宿主-病原交互效应的显著性排行。

## 功能概述 | Overview

- 实现 Wang et al. 2018 PNAS 的 ATOMM(双物种混合效应模型)，同时建模宿主 SNP、病原 SNP 及其交互对交叉感染表型的贡献
- 支持两种输入方式：直接给 ATOMM 格式的基因型/表型文件，或给 VCF + 交叉感染表型矩阵由工具自动转换
- 先估计方差组分(宿主/病原/交互/残差)，再做宿主边际、病原边际、宿主×病原交互三套 Score 检验
- 输出为扁平 TSV 表(边际/交互)+ 一个 Excel 汇总(需 openpyxl)；**不生成图**

## 快速开始 | Quick Start

```bash
biopytools atomm --host-vcf host.vcf.gz --pathogen-vcf pathogen.vcf.gz --phenotype-matrix pheno.tsv -o atomm_out
```

最小输入：宿主 VCF + 病原 VCF + 一张「行=宿主、列=病原」的 TSV 表型矩阵。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| 交叉感染(cross-infection) | 拿一批宿主和一批病原两两配对做感染实验，得到一张「宿主 × 病原」的成绩表 |
| SNP | 基因组上的一个「点位」，像高速公路上的里程桩，本工具逐桩扫描 |
| 边际效应(marginal effect) | 只算一方(宿主或病原)某个 SNP 单独对表型的贡献，像「单打独斗」 |
| 交互效应(interaction) | 宿主某 SNP × 病原某 SNP 联合起的额外作用，像「搭档配合」产生 1+1≠2 的效果 |
| 混合效应模型(mixed model) | 把「固定因素」和「随机因素」(如个体间的亲戚关系)一起放进统计模型里算 |
| 亲缘关系矩阵(GRM/kinship) | 个体之间基因相似度的表格；用于校正「长得像是因为有亲戚关系」带来的假信号 |
| 方差组分(variance components) | 表型差异被拆成几份：宿主基因占多少、病原基因占多少、交互占多少、随机噪声占多少 |
| P 值 / -log10(P) | 某效应「纯属巧合」的概率；P 越小证据越强，-log10(P) 越大越显著(3 约等于 P=0.001) |
| MAF | 「少数派碱基」在群体里占多少；太低=这个位点没信息量 |

## 输入 | Input

### 方式一：VCF + 表型矩阵(推荐，自动转换)

- `--host-vcf`：宿主 VCF(支持 .gz)，标准 VCF 格式。多等位位点、任一样本基因型缺失的位点会被跳过。
- `--pathogen-vcf`：病原 VCF，规则同上。
- `--phenotype-matrix`：交叉感染表型矩阵，**Tab 分隔**(注意：源码只认 Tab，逗号 CSV 会解析错)。行=宿主、列=病原：

```text
            pathogen1   pathogen2   pathogen3
host1       0.8         1.2         0.5
host2       NA          0.9         1.1
```

- 第一行是病原 ID，第一列是宿主 ID，单元格是表型值(如病情指数)。
- 缺失值默认用 `NA` 标记(可用 `--missing-value` 改)；`nan`/`NaN`/`.`/`-`/`NA` 恒被视为缺失，缺失的「宿主-病原」组合直接跳过不参与分析。

### 方式二：直接给 ATOMM 格式(进阶)

直接调用模块时可用 `-gh/-gp/-p`(click 包装器未暴露，一般走转换模式即可)：

- 基因型文件(宿主/病原同格式)：无表头，空白分隔，每行一个 SNP：`染色体ID  SNP_ID  样本1基因型  样本2基因型 …`，基因型为 0/1 数值。
- 表型文件：纯数值、空白分隔、无表头，列顺序固定为 `host_id  pathogen_id  [协变量…]  表型值`，ID 从 1 开始，须含截距列(通常全 1)。

## 参数说明 | Parameters

### 必需输入 | Required input

**通俗理解|In plain words:** 三个参数就是「宿主是谁、病原是谁、它们两两配对后结果如何」。用 `--host-vcf/--pathogen-vcf/--phenotype-matrix` 走自动转换最省事。

| 参数 | 说明 |
|------|------|
| `--host-vcf` | 宿主 VCF(支持 .gz) |
| `--pathogen-vcf` | 病原 VCF(支持 .gz) |
| `--phenotype-matrix` | 交叉感染表型矩阵(行=宿主、列=病原，Tab 分隔) |

### VCF 转换与编码 | VCF conversion & encoding

**通俗理解|In plain words:** `--encoding` 决定 VCF 里的基因型(0/1/2)怎么转成模型能用的数字。**默认 `auto` 会自动判断单倍体还是二倍体，一般不用动。** `--convert-maf` 只在 VCF 转换这一步过滤低频位点，和后面的 `--maf` 是两个独立开关。

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--encoding` | `auto` | 基因型编码：auto(自动判定)/haploid(单倍体)/dosage(剂量) |
| `--convert-maf` | `0.05` | VCF 转换时 MAF 阈值 |
| `--missing-value` | `NA` | 表型矩阵缺失值标记 |

### 过滤与检验范围 | Filtering & test ranges

**通俗理解|In plain words:** `--maf` 过滤「没信息量」的位点，默认 0.05 基本不用改。四个 `*-range` 用来限定「只测某一段 SNP」——想缩小范围、减少计算量或聚焦候选区时用，给 `起始 结束` 两个整数(0-based，含端点)。**交互检验必须同时给宿主和病原两个 range 才会跑**，否则只出边际效应结果。

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--maf` | `0.05` | MAF 过滤阈值(亲缘矩阵计算用) |
| `--host-snp-range` | 全部 | 宿主边际检验 SNP 范围(起 止) |
| `--pathogen-snp-range` | 全部 | 病原边际检验 SNP 范围(起 止) |
| `--interaction-host-range` | 全部 | 交互检验宿主 SNP 范围(起 止) |
| `--interaction-pathogen-range` | 全部 | 交互检验病原 SNP 范围(起 止) |

### 优化器 | Optimizer

**通俗理解|In plain words:** 估计方差组分用到一个数值优化器，`--tol` 是「算到多精细才停」、`--maxiter` 是「最多试几次」。**两个都一般不用动**，只有模型不收敛时才考虑放宽(调大 tol、调小 maxiter)。

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--tol` | `1e-6` | 优化容忍度 |
| `--maxiter` | `10000` | 优化器最大迭代次数 |

## 分析流程 | Pipeline

**通俗理解|In plain words:** 先把 VCF 转成模型格式、算出个体间「亲戚关系」，再估算「表型差异各由谁贡献多少」，最后逐个 SNP 做三套显著性检验写结果。

```text
宿主VCF + 病原VCF + 表型矩阵
    │
    ▼
(可选) 转换为 ATOMM 格式(编码 auto/haploid/dosage)
    │
    ▼
读取基因型/表型 → 计算亲缘关系矩阵(按 --maf 过滤)
    │
    ▼
估计方差组分(宿主/病原/交互/残差，Nelder-Mead 优化)
    │
    ▼
宿主边际检验 → 病原边际检验 → (可选)交互检验
    │
    ▼
写 TSV 结果 + Excel 汇总
```

## 输出 | Output

输出为**扁平结构**(无子目录、无图、无日志文件)，全部直接写在 `--output-dir` 下：

```text
atomm_out/
├── estimate.tsv            # 方差组分估计(宿主/病原/交互/残差)
├── marginal_host.tsv       # 宿主 SNP 边际效应
├── marginal_pathogen.tsv   # 病原 SNP 边际效应
├── interaction.tsv         # 宿主×病原交互效应(仅给两个交互 range 时生成)
├── atomm_results.xlsx      # 汇总 Excel(需 openpyxl，未装则跳过)
├── host_genotype.txt       # 转换模式生成的宿主基因型(ATOMM 格式)
├── pathogen_genotype.txt   # 转换模式生成的病原基因型(ATOMM 格式)
└── phenotype.txt           # 转换模式生成的表型(长格式 4 列)
```

## 结果解读 | Interpreting Results

### 1. 方差组分(`estimate.tsv`)

**通俗理解|In plain words:** 回答「感染结果里，宿主基因、病原基因、双方配合、随机噪声各占多少」。`Host`/`Pathogen`/`Interaction`/`Residual` 四行，数值越大贡献越大。交互组分接近 0 说明「配合」不重要，可能交互效应弱。

### 2. 边际效应(`marginal_host.tsv` / `marginal_pathogen.tsv`)

**通俗理解|In plain words:** 这是「单打独斗」的成绩单——宿主的每个 SNP(或病原的每个 SNP)单独对表型有没有作用。

列：`Host_Chr  Host_POS  Statistic  Pvalue  NegLog10P`(病原表同理)。`Pvalue` 越小、`NegLog10P` 越大越显著；NegLog10P ≥ 3(即 P ≤ 0.001)通常视为候选信号。`Statistic` 是 Score 检验统计量，越大效应越强。

### 3. 交互效应(`interaction.tsv`)

**通俗理解|In plain words:** 这是「搭档配合」的成绩单——宿主的某个 SNP 和病原的某个 SNP 凑在一起，会不会产生「1+1 远大于 2」的效果。这类信号往往是「宿主的某个抗病基因恰好遇上病原的某个效应因子」这类机制的候选。

列：`Host_Chr  Host_POS  Pathogen_Chr  Pathogen_POS  Statistic  Pvalue  NegLog10P`。交互检验是全基因组两两组合、计算量巨大，**建议先用 `--interaction-host-range`/`--interaction-pathogen-range` 限定候选区**再跑。

### 4. Excel 汇总(`atomm_results.xlsx`)

4 个 Sheet：`Heritability`(方差组分+占比百分比)、`Host_Marginal`、`Pathogen_Marginal`、`Interaction`，内容与上面 TSV 一致，方便 Excel 用户筛选排序。

## 参数选择建议 | Parameter Guidance

| 场景 | 建议 |
|------|------|
| 常规分析 | 全部默认，直接跑 |
| 只关心单边效应(不查交互) | 不给 `--interaction-*-range` 即可，交互检验自动跳过 |
| 想查交互但位点太多 | 先用边际效应结果圈候选区，再用 `--interaction-host-range`/`--interaction-pathogen-range` 限定范围 |
| 单倍体病原/宿主 | `--encoding haploid`；不确定就默认 `auto` |
| 表型矩阵用逗号分隔 | **不支持**，须转成 Tab 分隔 |
| 想省内存/加速 | 用 `--host-snp-range`/`--pathogen-snp-range` 只测一段 |

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--host-vcf` | 必填 |  | 宿主VCF文件(支持.gz)｜Host VCF file (.gz supported) |
| `--pathogen-vcf` | 必填 |  | 病原VCF文件(支持.gz)｜Pathogen VCF file (.gz supported) |
| `--phenotype-matrix` | 必填 |  | 交叉感染表型矩阵(行=宿主,列=病原)｜Cross-infection phenotype matrix (rows=hosts, cols=pathogens) |
| `-o, --output-dir` | `./output` |  | 输出目录｜Output directory |
| `--maf` | `0.05` | float | MAF过滤阈值｜MAF filter threshold |
| `--encoding` | `auto` | auto/haploid/dosage | 基因型编码方式｜Genotype encoding mode |
| `--convert-maf` | `0.05` | float | VCF转换时MAF阈值｜MAF threshold for VCF conversion |
| `--missing-value` | `NA` |  | 表型缺失值标记｜Missing value marker in phenotype matrix |
| `--host-snp-range` | — | int | 宿主边际检验SNP范围｜Host marginal test SNP range |
| `--pathogen-snp-range` | — | int | 病原边际检验SNP范围｜Pathogen marginal test SNP range |
| `--interaction-host-range` | — | int | 交互检验宿主SNP范围｜Interaction test host SNP range |
| `--interaction-pathogen-range` | — | int | 交互检验病原SNP范围｜Interaction test pathogen SNP range |
| `--tol` | `1e-06` | float | 优化容忍度｜Optimizer tolerance |
| `--maxiter` | `10000` | int | 最大迭代次数｜Max iterations |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-gh, --host-genotype` | — |  | 宿主基因型文件(ATOMM格式)｜Host genotype file (ATOMM format) |
| `-gp, --pathogen-genotype` | — |  | 病原基因型文件(ATOMM格式)｜Pathogen genotype file (ATOMM format) |
| `-p, --phenotype` | — |  | 表型文件(ATOMM格式)｜Phenotype file (ATOMM format) |
| `--host-vcf` | — |  | 宿主VCF文件(自动转换为ATOMM格式)｜Host VCF file (auto-convert to ATOMM format) |
| `--pathogen-vcf` | — |  | 病原VCF文件(自动转换为ATOMM格式)｜Pathogen VCF file (auto-convert to ATOMM format) |
| `--phenotype-matrix` | — |  | 交叉感染表型矩阵(行=宿主,列=病原)｜Cross-infection phenotype matrix |
| `--encoding` | `auto` | auto/haploid/dosage | 基因型编码方式｜Genotype encoding mode (default: auto) |
| `--convert-maf` | `0.05` | float | VCF转换时MAF阈值｜MAF threshold for VCF conversion (default: 0.05) |
| `--missing-value` | `NA` |  | 表型缺失值标记｜Missing value marker (default: NA) |
| `-o, --output-dir` | `./output` |  | 输出目录｜Output directory |
| `--maf` | `0.05` | float | MAF过滤阈值｜MAF filter threshold (default: 0.05) |
| `--host-snp-range` | — | int | 宿主边际检验SNP范围(0-based)｜Host marginal test SNP range |
| `--pathogen-snp-range` | — | int | 病原边际检验SNP范围(0-based)｜Pathogen marginal test SNP range |
| `--interaction-host-range` | — | int | 交互检验宿主SNP范围｜Interaction test host SNP range |
| `--interaction-pathogen-range` | — | int | 交互检验病原SNP范围｜Interaction test pathogen SNP range |
| `--tol` | `1e-06` | float | 优化容忍度｜Optimizer tolerance (default: 1e-6) |
| `--maxiter` | `10000` | int | 优化器最大迭代次数｜Max optimizer iterations (default: 10000) |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- Python 3，第三方库：`numpy`(必需)、`scipy`(必需，优化+卡方检验)、`openpyxl`(可选，仅 Excel 输出)
- 无外部软件依赖(不调用 bcftools 等，VCF 纯 Python 解析)

## 常见问题 | FAQ

**Q1：为什么没有生成 interaction.tsv？**
交互检验只有在**同时**提供 `--interaction-host-range` 和 `--interaction-pathogen-range` 时才执行；两个都给才会生成交互结果。

**Q2：表型矩阵是 CSV(逗号分隔)能直接用吗？**
不能。源码只按 Tab 分隔解析，逗号分隔会解析错。请转成 Tab 分隔的 .tsv。

**Q3：`--encoding dosage` 是不是把基因型 0/1/2 除以 2？**
不是。dosage 模式下 VCF 转换保留 0/1/2(REF 计数)；只有 `auto` 判定为二倍体时才会除以 2 得 0/0.5/1。若对编码有疑问，直接看转换产物 `host_genotype.txt`。

**Q4：输出在哪？有没有图？**
输出是输出目录下的扁平 TSV 文件 + `atomm_results.xlsx`，**不生成任何图**(无曼哈顿图)，亲缘矩阵也不落盘。

**Q5：重跑会跳过已完成步骤吗？**
不会。本工具无断点续传，重跑会完整重算并覆盖所有输出文件。
