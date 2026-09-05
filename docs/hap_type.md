# 单倍型提取 | Haplotype Extraction (geneHapR-compatible)

一句话理解：**从 VCF 的某个区间里，把"每个样本在这段 DNA 上实际携带的碱基串"提取出来，把携带相同碱基串的样本归为一组（单倍型）**，输出 geneHapR 兼容的表格。

## 功能概述 | Overview

- 从 VCF 指定区间提取单倍型，支持单区间（`chr:start-end`）或 BED 文件批量处理
- 基于 geneHapR 的 vcf2hap 算法：每个样本的基因型拼成字符串，模式相同的样本归为同一单倍型，按频率排序编号 H001/H002…
- 支持去杂合位点（`--hetero-remove`）、去缺失位点（`--na-drop`，默认开）
- 输出 hapResult / hapSummary / sampleHap 三种表格（txt + xlsx 各一份）
- 用 cyvcf2 读取 VCF；无 tabix 索引时自动回退全文件遍历

## 快速开始 | Quick Start

```bash
biopytools hap-type -i sample.vcf -r chr1:1000-5000 -o result
```

最小输入：一个 VCF + 一个区间（或 BED 文件）。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| 单倍型(haplotype) | 一条染色体上多个位点碱基的"组合"，像一段固定的"字母串" |
| VCF | 记录所有样本在每个位点变异（基因型）的标准文件 |
| 基因型(genotype) | 一个位点上两条染色体的构成，如 0/0、0/1、1/1 |
| 杂合(heterozygous) | 两条染色体在该位点不同（0/1） |
| 缺失(missing) | 该位点没测出来（./.） |
| geneHapR | 下游单倍型分析 R 包，本工具输出它的标准输入格式 |

## 输入 | Input

- `-i/--vcf`：VCF 变异文件（支持 .vcf / .vcf.gz）
- `-r/--region`：单个区间 `chr:start-end`（1-based 闭区间），或一个 BED 文件（BED3/4 列自动识别，start 为 0-based 自动转 1-based）

BED 文件示例（第 4 列是可选名称）：

```text
chr1    1000    5000    geneA
chr1    8000    12000   geneB
```

## 参数说明 | Parameters

### 必需参数 | Required

**通俗理解|In plain words:** `-i` 是变异来源，`-r` 决定"看哪一段"。单个区间直接写 `chr:start-end`；要一次看很多段就写一个 BED 文件路径（每行一段）。**VCF 需要能按区间随机访问，最好有 tabix 索引**（没有也能跑，只是更慢）。

### 输出参数 | Output

**通俗理解|In plain words:** `-o` 是输出文件前缀（或目录）。不填时写到当前目录。单区间模式输出文件名自带 `chr_start_end` 前缀，BED 模式每个区间各生成一组文件、日志共用一份。**文件名由程序自动生成，一般不用操心**。

### 位点过滤 | Site filtering

**通俗理解|In plain words:** 这两个开关决定"哪些位点不参与单倍型判定"。`--na-drop`（默认开）去掉有缺失的位点——某个样本在该位点没测出来，整列都信不过；`--hetero-remove`（默认关）去掉含杂合的位点——杂合位点无法确定这条染色体的"唯一字母"，去掉后单倍型更"纯"。**一般保持默认**；样本高度杂合、想得到更清晰的单倍型分组时再开 `--hetero-remove`。

### 编号参数 | Numbering

**通俗理解|In plain words:** `--hap-prefix` 是单倍型 ID 的前缀（默认 H），`--pad` 是编号位数（默认 3，即 H001）。**一般不用动**；单倍型超过 999 个时编号会自然变长，不受影响。

## 分析流程 | Pipeline

```text
VCF + 区间(或 BED)
    │
    ▼
cyvcf2 读取区间内变异位点（无索引则回退全文件遍历）
    │
    ▼
基因型编码（纯合=单字母，杂合=A|T，缺失=N）
    │
    ▼
可选过滤：去杂合位点(--hetero-remove) / 去缺失位点(--na-drop)
    │
    ▼
每个样本的基因型拼成字符串，按模式分组
    │
    ▼
按频率降序排序，分配 H001/H002… 编号
    │
    ▼
导出 hapResult / hapSummary / sampleHap（txt + xlsx）
```

## 输出 | Output

单区间模式（`-r chr1:1000-5000`）：

```text
chr1_1000_5000_hapResult.txt      # 每行=一个样本：单倍型ID + 各位点碱基 + 样本名
chr1_1000_5000_hapResult.xlsx
chr1_1000_5000_hapSummary.txt     # 每行=一个唯一单倍型：ID + 碱基串 + 样本名单 + 频率
chr1_1000_5000_hapSummary.xlsx
chr1_1000_5000_sampleHap.txt      # 样本→单倍型映射表
chr1_1000_5000_sampleHap.xlsx
chr1_1000_5000.log                # 运行日志
```

BED 批量模式：每个区间生成如上的一组文件（前缀为 `chr_start_end_name`），日志共用 `hap_type.log`。

## 结果解读 | Interpreting Results

### 1. hapResult（每行=样本）

**通俗理解|In plain words:** 一行为一个样本，第一列是它所属的单倍型 ID，中间是它在每个位点的碱基，最后一列是样本名。同 ID 的样本 = 这段 DNA 携带完全相同的碱基串。

### 2. hapSummary（每行=唯一单倍型）

**通俗理解|In plain words:** 把 hapResult 去重、按频率排序后的"单倍型清单"，最后一列 freq 是该单倍型的样本数。**这是最常用的一张表**——频率最高的几个单倍型就是这段区间的主流基因型。

### 3. sampleHap（样本→单倍型映射）

**通俗理解|In plain words:** 两列：样本名、单倍型 ID。用来快速查"某个样本属于哪个单倍型"，可直接喂给 geneHapR 做下游分析。

### 4. 好坏判据

- 单倍型数（n_haplotypes）远小于样本数：群体结构清晰，主流单倍型明显
- 过滤后变异位点数为 0：区间太窄或变异太稀，换更大区间或关掉过滤
- 某区间被跳过（日志警告"无变异位点"）：该区间没有可用的变异，属正常

## 参数选择建议 | Parameter Guidance

- 一般场景：全默认，先看 hapSummary 的单倍型分布
- 想要"更纯"的单倍型：开 `--hetero-remove`（去除含杂合位点，编号更稳定）
- 缺失位点较多：保持 `--na-drop` 默认开，避免缺失污染分组
- 想自定义编号：`--hap-prefix H --pad 3` 控制前缀与位数
- 一次分析多个基因区间：写 BED 文件，比逐个跑更省事

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --vcf` | 必填 |  | VCF变异文件｜VCF variant file |
| `-r, --region` | 必填 |  | 基因组区间(chr:start-end)或BED文件｜Genomic interval or BED file |
| `-o, --output` | — |  | 输出文件前缀(可选，默认自动生成)｜Output file prefix (auto-generated if omitted) |
| `--hetero-remove/--no-hetero-remove` | `False` |  | 去除杂合位点｜Remove heterozygous sites |
| `--na-drop/--no-na-drop` | `True` |  | 去除缺失位点｜Remove missing sites |
| `--hap-prefix` | `H` |  | 单倍型ID前缀｜Haplotype ID prefix |
| `--pad` | `3` | int | 单倍型ID位数｜Haplotype ID padding digits |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- Python 3
- cyvcf2（读取 VCF，按区间随机访问）
- openpyxl（写 xlsx 文件）
- 可选：bcftools（tabix 索引缺失时加速；无索引也可运行；自动解析 align 域环境并经 conda run 调用，可用环境变量 BCFTOOLS_PATH 覆盖；域环境缺失时回退 PATH 直接调用）

## 常见问题 | FAQ

**Q1：区间没有变异怎么办？**
日志会警告"区间无变异位点，跳过"，该区间不产生文件。换更大的区间或换一个变异更密集的区域。

**Q2：跑得很慢？**
大概率是 VCF 没有 tabix 索引，程序回退到全文件遍历。对 VCF 建索引（`bcftools index sample.vcf.gz`）后按区间随机访问会快很多。

**Q3：输出写到哪了？**
未指定 `-o` 时写到当前工作目录，文件名以 `chr_start_end` 开头。BED 模式日志统一在输出目录下的 `hap_type.log`。

**Q4：为什么 hapSummary 里的频率加起来不等于样本数？**
频率=携带该单倍型的样本数，每个样本恰好归入一个单倍型，所以频率之和=样本总数。若对不上，检查是否有样本被过滤（如缺失过多）。

**Q5：能不能去重/去杂合后保留原位点信息？**
可以。hapResult 和 hapSummary 的前几行（CHR/POS/INFO/ALLELE）记录了过滤后保留位点的坐标和 REF/ALT，对照即可。
