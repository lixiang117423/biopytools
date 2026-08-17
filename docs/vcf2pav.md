# VCF 转 PAV 矩阵 | VCF to PAV Matrix

一句话理解：**把合并后的结构变异 VCF 转成一张 0/1 表——行是每个 SV，列是每个样本，1=有这个变异，0=没有**，方便后续做群体比较、聚类或进化分析。

## 功能概述 | Overview { #overview }

- 把 VCF 转成 PAV（Presence/Absence Variation，存在/缺失变异）矩阵
- 输出矩阵：行=SV，列=样本，值=0/1；附带每个样本的统计摘要
- 按 GT 解码存在/缺失：0/0、./. 记为 0，含替代等位基因记为 1
- 支持多等位位点，每个 ALT 等位基因单独成一行
- 纯 Python，无外部依赖

## 快速开始 | Quick Start { #quick-start }

```bash
biopytools vcf2pav -i pan_sv.survivor.vcf -o output_dir
```

最小输入：一个 VCF 文件（如 SURVIVOR 合并后的 SV VCF）。

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语<br>Term | 通俗理解 |
|---|---|
| PAV | 存在/缺失变异，某段序列这个样本有、那个样本没有 |
| SVTYPE | 结构变异类型：DEL 缺失 / INS 插入 / INV 倒位 / DUP 重复 / BND 易位 |
| GT | 基因型字段，决定样本是否「拥有」该变异 |
| SURVIVOR VCF | 常见 SV 合并工具的输出，带 SVTYPE/END/SVLEN 等 INFO 字段 |

## 输入 | Input { #input }

标准 VCF（含 `#CHROM` 头、FORMAT/GT、样本列），优先读 INFO 中的 `SVTYPE`、`END`、`SVLEN`、`CHR2`。典型为 SURVIVOR 合并产物：

```text
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1	sample2
Chr01	100	.	A	<DEL>	.	PASS	SVTYPE=DEL;END=1000;SVLEN=-900	GT:...	0/0	1/1
```

## 参数说明 | Parameters { #parameters }

### 必需参数 | Required

**通俗理解|In plain words:** `-i` 输入 VCF；`-o` 输出目录（默认 `./vcf2pav_output`）。

相关参数：`-i/--input`、`-o/--output-dir`（默认 `./vcf2pav_output`）。

### 日志 | Logging

**通俗理解|In plain words:** `--log-level` 控制日志级别，一般不用动。

相关参数：`--log-level`（默认 INFO）。

## 分析流程 | Pipeline { #pipeline }

```text
VCF 文件
    |
    v
第一遍: 解析 #CHROM 头收集样本名
    |
    v
第二遍: 逐行解析 SV + GT,构建 0/1 矩阵
    |
    v
写 pav_matrix.tsv + pav_summary.tsv
```

## 输出 | Output { #output }

```text
vcf2pav_output/
├── pav_matrix.tsv      # PAV 矩阵(行=SV,列=样本,值=0/1)
├── pav_summary.tsv     # 每样本的 SV 类型计数与存在/缺失统计
└── vcf2pav.log         # 运行日志
```

`pav_matrix.tsv` 列：`SV_ID、SVTYPE、SVLEN、CHROM、START、END、样本列…`。

- SV_ID 格式：跨染色体 BND/TRA 为 `{CHROM}_{START}-{CHR2}_{END}_{SVTYPE}`，其余为 `{CHROM}_{START}-{END}_{SVTYPE}`
- 每个 ALT 等位基因单独成一行

## 结果解读 | Interpreting Results { #interpreting-results }

- `pav_matrix.tsv`：1=样本含该 SV，0=不含。可直接转成聚类/系统发育输入
- `pav_summary.tsv`：每样本 total_sv、present、absent 及各 SVTYPE（DEL/INS/INV/DUP/BND/TRA/UNK）计数，可核对不同样本的变异负载差异
- 无法解析的 GT 记为 `NaN`，日志会给出数量

## 参数选择建议 | Parameter Guidance { #parameter-guidance }

- 标准 SURVIVOR 合并 VCF：直接默认参数
- 只关心特定类型：下游按 `SVTYPE` 列过滤即可（本工具不做类型过滤）

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 输入 VCF 文件｜Input VCF file |
| `-o, --output-dir` | `./vcf2pav_output` |  | 输出目录(默认./vcf2pav_output)｜Output directory |
| `--log-level` | `INFO` |  | 日志级别(默认INFO)｜Log level (default INFO) |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 输入 VCF 文件｜Input VCF file |
| `-o, --output-dir` | `./vcf2pav_output` |  | 输出目录｜Output directory (default ./vcf2pav_output) |
| `--log-level` | `INFO` |  | 日志级别｜Log level (default INFO) |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies { #dependencies }

- Python 3 标准库（csv）
- 无外部生信软件、无 conda 环境依赖

## 常见问题 | FAQ { #faq }

**Q1：会断点续传吗？**
不会。每次运行重新解析并覆盖输出。

**Q2：什么样的 VCF 能用？**
需含 `#CHROM` 头和 GT 字段；SV 类型优先读 INFO 的 SVTYPE，缺失时对跨染色体 BND/TRA 用 CHR2 构造 ID。

**Q3：缺失基因型怎么算？**
`./.`、`.`、`0/0`、`0|0` 都记为 0（缺失），`1/1`、`0/1`、`1/0` 记为 1（存在）。
