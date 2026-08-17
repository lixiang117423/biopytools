# VCF 基因型提取 | VCF Genotype Extraction

一句话理解：**把 VCF 里每个变异位点、每个样本的基因型（0/0、0/1、1/1）抽出来整理成一张表格**，可按样本、变异长度、双等位等条件筛选，输出 txt 或 csv。

## 功能概述 | Overview { #overview }

- 从 VCF 提取每个样本的基因型，输出行=变异、列=样本的表格
- 附带每位点的纯合/杂合比例、各基因型计数列
- 支持样本筛选、双等位过滤、变异长度过滤、按染色体拆分输出
- 输出 txt 或 csv；两遍扫描流式写出，内存占用低
- 优先用 cyvcf2 加速，缺失时自动回退原生 Python 解析

## 快速开始 | Quick Start { #quick-start }

```bash
biopytools vcf2genotype -i variants.vcf -o genotypes
```

最小输入：一个 VCF 文件。

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语<br>Term | 通俗理解 |
|---|---|
| 基因型(GT) | 一个位点上两条染色体的构成，0/0=纯合参考，0/1=杂合，1/1=纯合替代 |
| 双等位位点 | 该位点只有一种替代等位基因（ALT 不含逗号） |
| 变异长度 | REF 与 ALT 的长度差；SNP 为 0，插入/缺失为正数 |
| 纯合比例 | 该位点样本里「两个等位相同」的占比 |
| 缺失基因型 | 该样本这个位点没测出来，记为 `./.` |

## 输入 | Input { #input }

标准 VCF 文件（`.vcf` 或 `.vcf.gz`），需含 FORMAT/GT 和样本列。

## 参数说明 | Parameters { #parameters }

### 输入输出 | Input & output

**通俗理解|In plain words:** `-i` 输入 VCF；`-o` 输出前缀，`-t` 选 txt/csv，`--output-dir` 指定输出目录。注意默认输出前缀实际是 `vcf_genotype`（见 FAQ）。

相关参数：`-i/--input`、`-o/--output`、`-t/--output-type`（默认 txt）、`--output-dir`（默认 `./`）。

### 样本与过滤 | Samples & filtering

**通俗理解|In plain words:** `-s` 选样本（默认 all 全选，或逗号分隔若干样本）；`--biallelic-only` 只要双等位位点；`--min-length`/`--max-length` 按变异长度过滤（默认不限）。**一般只动 `-s`，其余按需。**

相关参数：`-s/--samples`（默认 all）、`--biallelic-only`、`--min-length`、`--max-length`。

### 输出组织 | Output organization

**通俗理解|In plain words:** `-e/--each` 设为 yes/y 时按染色体拆成多个文件（一个染色体一个文件），否则全部写进一个 `_all` 文件。

相关参数：`-e/--each`（默认 n）。

### 日志与调试 | Logging & debug

**通俗理解|In plain words:** `--dry-run` 只扫描统计不写数据文件（快速预览有多少变异、分布在哪些染色体）；`-v` 提高日志详细度、`--quiet` 静默、`--log-file`/`--log-level` 控制日志落盘与级别。**平时不用动。**

相关参数：`--dry-run`、`-v/--verbose`、`--quiet`、`--log-file`、`--log-level`（默认 INFO）。

## 分析流程 | Pipeline { #pipeline }

```text
VCF 文件
    |
    v
第一遍: 扫描收集基因型类型 + 过滤统计
    |
    v
第二遍: 流式逐行写出基因型表
    |
    v
写汇总文件
```

## 输出 | Output { #output }

```text
vcf_genotype_all.txt          # 全部变异(默认前缀 vcf_genotype,后缀 _all)
vcf_genotype_summary.txt      # 汇总统计
99_logs/vcf_extraction_YYYYMMDD_HHMMSS.log   # 日志
```

- 按染色体拆分（`-e yes`）时为 `vcf_genotype_{染色体}.txt`
- `-t csv` 时扩展名为 `.csv`

主表列顺序：`CHROM、POS、ID、REF、ALT、QUAL、Homozygous_Ratio、Heterozygous_Ratio、GT_*（各基因型计数）、样本列…`

## 结果解读 | Interpreting Results { #interpreting-results }

- 每行一个变异位点，样本列给出该样本基因型（`0/0`、`0/1`、`1/1`、`./.`）
- `GT_0/0`、`GT_0/1`、`GT_1/1`、`GT_./.` 等列是该位点各基因型的样本数
- `Homozygous_Ratio`/`Heterozygous_Ratio` 是纯合/杂合比例（0~1）
- `summary.txt` 给出总变异数、染色体数与各染色体变异数，可快速核对数据规模

## 参数选择建议 | Parameter Guidance { #parameter-guidance }

- 拿全量基因型做关联/GWAS：默认即可
- 只要部分样本：`-s sample1,sample2`
- 只做 SNP 分析：`--biallelic-only`（或再加长度过滤）
- 快速预览数据规模：`--dry-run` 不落数据文件

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input, -i` | 必填 | Path | 输入VCF文件路径｜Input VCF file path |
| `--output, -o` | `vcf2genotype` | str | 输出文件前缀｜Output file prefix |
| `--samples, -s` | `all` | str | 样本选择：all(所有样本)或逗号分隔的样本名称｜Sample selection: all or comma-separated sample names |
| `--biallelic-only` | — |  | 只保留双等位位点｜Keep only biallelic sites |
| `--min-length` | — | int | 最小变异长度（包含），默认不限制｜Minimum variant length (inclusive), default no limit |
| `--max-length` | — | int | 最大变异长度（包含），默认不限制｜Maximum variant length (inclusive), default no limit |
| `--each, -e` | `n` | yes/y/no/n | 按染色体拆分输出文件：yes/y或no/n｜Split output files by chromosome: yes/y or no/n |
| `--output-type, -t` | `txt` | txt/csv | 输出文件格式｜Output file format |
| `--output-dir` | `./` | Path | 输出目录｜Output directory |
| `--verbose, -v` | — |  | 增加输出详细程度｜Increase output verbosity |
| `--quiet` | — |  | 静默模式｜Quiet mode |
| `--log-file` | — | str | 日志文件路径｜Log file path |
| `--log-level` | `INFO` | DEBUG/INFO/WARNING/ERROR/CRITICAL | 日志级别｜Log level |
| `--dry-run` | — |  | 试运行模式｜Dry run mode |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | VCF文件路径(支持.gz压缩格式)｜VCF file path (supports .gz compressed format) |
| `-o, --output` | `vcf_genotype` |  | 输出文件前缀｜Output file prefix |
| `-t, --output-type` | `txt` | txt/csv | 输出文件格式｜Output file format |
| `--output-dir` | `./` |  | 输出目录｜Output directory |
| `-s, --samples` | `all` |  | 样本选择：all（所有样本）或逗号分隔的样本名称｜Sample selection: all (all samples) or comma-separated sample names |
| `--biallelic-only` | — | store_true | 只保留双等位位点｜Keep only biallelic sites |
| `--min-length` | — | int | 最小变异长度（包含），默认不限制｜Minimum variant length (inclusive), default no limit |
| `--max-length` | — | int | 最大变异长度（包含），默认不限制｜Maximum variant length (inclusive), default no limit |
| `-e, --each` | `n` | yes/y/no/n | 按染色体拆分输出文件：yes/y（是）或no/n（否）｜Split output files by chromosome: yes/y or no/n |
| `--dry-run` | — | store_true | 试运行模式，不实际执行操作｜Dry run mode, no actual operations performed |
| `-v, --verbose` | `0` | count | 增加输出详细程度 (-v, -vv, -vvv)｜Increase output verbosity (-v, -vv, -vvv) |
| `--quiet` | — | store_true | 静默模式，仅输出错误信息｜Quiet mode, only output error messages |
| `--log-file` | — | str | 日志文件路径｜Log file path |
| `--log-level` | `INFO` | DEBUG/INFO/WARNING/ERROR/CRITICAL | 日志级别｜Log level |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies { #dependencies }

- Python 3 标准库
- cyvcf2（可选，加速解析；缺失时自动回退原生解析器）
- 无外部生信软件、无 conda 环境依赖

## 常见问题 | FAQ { #faq }

**Q1：会断点续传吗？**
不会。每次运行都重新扫描并覆盖输出。

**Q2：为什么默认输出文件叫 `vcf_genotype_all.txt` 而不是 `vcf2genotype_all.txt`？**
CLI 包装器的 `-o` 默认值是 `vcf2genotype`，但模块 main.py 的默认前缀是 `vcf_genotype`，且包装器只在值非默认时才透传 `-o`，所以不显式传 `-o` 时实际用的是 `vcf_genotype`。想自定义就用 `-o 你的前缀`。

**Q3：没装 cyvcf2 会怎样？**
自动回退到原生 Python 解析器，结果一致，只是速度慢一些。

**Q4：缺失基因型怎么表示？**
记为 `./.`，且会计入 `GT_./.` 列。
