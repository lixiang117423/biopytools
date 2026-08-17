# VCF 筛选 | VCF Filtering

一句话理解：**从一个 VCF 里按染色体/区间/质量/样本等条件，快速把不要的变异和样本筛掉**，核心用 bcftools 高速执行。

## 功能概述 | Overview { #overview }

- 按染色体、坐标区间筛选变异（支持多染色体逗号分隔）
- 质量过滤：MAF、最大缺失率、QUAL 阈值、双等位、去 indel
- 样本筛选：keep-samples / remove-samples
- 核心用 bcftools 引擎（view + filter），速度快
- 支持 .vcf/.vcf.gz（.gz 缺索引时自动生成 tabix 索引）

## 快速开始 | Quick Start { #quick-start }

```bash
biopytools vcf-filter -i input.vcf -c chr1 -s 1000 -e 2000 -o filtered.vcf
```

最小输入：VCF + 染色体名（`-c` 必填）。

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语<br>Term | 通俗理解 |
|---|---|
| MAF | 次要等位基因频率，太低说明该位点在群体里几乎没差异 |
| 缺失率 | 多少样本在这个位点没测出来（F_MISSING） |
| QUAL | 位点质量的打分，越低越不可信 |
| 双等位位点 | 只有一种替代等位基因的位点 |
| indel | 插入/缺失变异（相对 SNP 而言） |

## 输入 | Input { #input }

标准 VCF（`.vcf` 或 `.vcf.gz`）。`-c/--chr` 必填，指定染色体名（与 VCF 中 CHROM 列完全一致，含前缀），多个用逗号分隔。

## 参数说明 | Parameters { #parameters }

### 输入输出 | Input & output

**通俗理解|In plain words:** `-i` 输入 VCF；`-o` 输出文件（不传则自动命名为 `{输入文件名去掉扩展名}_filtered.vcf`）。

相关参数：`-i/--input`、`-o/--output`。

### 位置筛选 | Region

**通俗理解|In plain words:** `-c` 指定染色体（必填，可多个逗号分隔）；`-s`/`-e` 是区间起止（可选，只给一个也可）。三者合起来就是「取哪条染色体的哪一段」。

相关参数：`-c/--chr`、`-s/--start`、`-e/--end`。

### 质量过滤 | Quality

**通俗理解|In plain words:** `--maf` 设最小等位基因频率（低于就删）；`--max-missing` 设最大缺失率（高于就删）；`--quality-threshold` 设 QUAL 最低分；`--biallelic-only` 只要双等位；`--remove-indels` 删掉插入缺失只留 SNP。**都不传就是只做区域筛选，不做质量过滤。**

相关参数：`--maf`、`--max-missing`、`--quality-threshold`、`--biallelic-only`、`--remove-indels`。

### 样本筛选 | Sample filtering

**通俗理解|In plain words:** `--keep-samples`/`--remove-samples` 按样本名保留/删除（逗号分隔）。

相关参数：`--keep-samples`、`--remove-samples`。

### 性能与验证 | Performance & validation

**通俗理解|In plain words:** `--skip-validation`（默认开）跳过输入预验证提速；`--force-validation` 强制验证；`-v` 显示详细信息。**一般保持默认即可。**

相关参数：`--skip-validation`（默认开）、`--force-validation`、`-v/--verbose`。

## 分析流程 | Pipeline { #pipeline }

```text
VCF 文件
    |
    v
bcftools view(区域/样本筛选)
    |
    v
bcftools filter(--include 质量表达式)
    |
    v
输出筛选后的 VCF
```

## 输出 | Output { #output }

- 筛选后的 VCF 文件（`-o` 指定，或默认 `{输入文件名去掉扩展名}_filtered.vcf`）
- 日志：输出目录下的 `vcf_filter.log`
- 样本筛选时可能在输出目录旁生成临时 `.keep_samples.txt` / `.remove_samples.txt`

## 结果解读 | Interpreting Results { #interpreting-results }

- 输出 VCF 即筛选后保留的变异，可直接做下游分析
- 若筛选后为空，多半是染色体名不匹配（与 VCF 的 CHROM 列逐字对照，注意前缀）或过滤条件过严
- 日志中会打印执行的具体 bcftools 命令，可用于核对过滤逻辑

## 参数选择建议 | Parameter Guidance { #parameter-guidance }

- 只取某区域：`-c chr1 -s 1000 -e 2000`（不加质量过滤参数）
- 群体分析前常规质控：`--maf 0.05 --max-missing 0.1 --biallelic-only`
- 只做 SNP：`--remove-indels`
- 只留若干样本：`--keep-samples s1,s2`

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input, -i` | 必填 |  | 输入VCF文件路径｜Input VCF file path |
| `--output, -o` | — | Path | 输出VCF文件路径｜Output VCF file path |
| `--chr, -c` | 必填 | str | 染色体名称(支持逗号分隔多个)｜Chromosome name(s) (comma-separated for multiple) |
| `--start, -s` | — | int | 起始位置｜Start position |
| `--end, -e` | — | int | 结束位置｜End position |
| `--convert-format` | — |  | 使用PLINK进行格式转换｜Use PLINK for format conversion |
| `--plink-path` | `plink` | str | PLINK可执行文件路径｜PLINK executable path |
| `--allow-extra-chr` | `True` |  | 允许额外染色体｜Allow extra chromosomes |
| `--maf` | — | float | 最小等位基因频率｜Minimum allele frequency |
| `--max-missing` | — | float | 最大缺失率｜Maximum missing rate |
| `--quality-threshold` | — | float | 质量阈值｜Quality threshold |
| `--min-depth` | — | int | 最小深度｜Minimum depth |
| `--max-depth` | — | int | 最大深度｜Maximum depth |
| `--keep-samples` | — | str | 保留样本名称(逗号分隔)｜Sample names to keep (comma-separated) |
| `--remove-samples` | — | str | 移除样本名称(逗号分隔)｜Sample names to remove (comma-separated) |
| `--keep-ids` | — | str | 保留变异位点ID(逗号分隔)｜Variant IDs to keep (comma-separated) |
| `--remove-ids` | — | str | 移除变异位点ID(逗号分隔)｜Variant IDs to remove (comma-separated) |
| `--biallelic-only` | — |  | 只保留双等位基因位点｜Keep only biallelic sites |
| `--remove-indels` | — |  | 移除插入缺失变异｜Remove indel variants |
| `--skip-validation` | `True` |  | 跳过输入验证以提高速度｜Skip input validation for speed |
| `--force-validation` | — |  | 强制执行输入验证｜Force input validation |
| `--verbose, -v` | — |  | 显示详细信息｜Show verbose information |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 输入VCF文件路径｜Input VCF file path |
| `-o, --output` | — |  | 输出VCF文件路径｜Output VCF file path |
| `-c, --chr, --chromosome` | 必填 |  | 染色体名称 (支持逗号分隔的多个染色体)｜Chromosome name(s) (comma-separated for multiple) |
| `-s, --start` | — | int | 起始位置｜Start position |
| `-e, --end` | — | int | 结束位置｜End position |
| `--convert-format` | — | store_true | 使用PLINK进行格式转换｜Use PLINK for format conversion |
| `--plink-path` | `plink` |  | PLINK可执行文件路径｜PLINK executable path |
| `--allow-extra-chr` | `True` | store_true | 允许额外染色体｜Allow extra chromosomes |
| `--maf` | — | float | 最小等位基因频率｜Minimum allele frequency |
| `--max-missing` | — | float | 最大缺失率｜Maximum missing rate |
| `--quality-threshold` | — | float | 质量阈值｜Quality threshold |
| `--min-depth` | — | int | 最小深度｜Minimum depth |
| `--max-depth` | — | int | 最大深度｜Maximum depth |
| `--keep-samples` | — |  | 保留的样本名称 (逗号分隔)｜Sample names to keep (comma-separated) |
| `--remove-samples` | — |  | 移除的样本名称 (逗号分隔)｜Sample names to remove (comma-separated) |
| `--keep-ids` | — |  | 保留的变异位点ID (逗号分隔)｜Variant IDs to keep (comma-separated) |
| `--remove-ids` | — |  | 移除的变异位点ID (逗号分隔)｜Variant IDs to remove (comma-separated) |
| `--biallelic-only` | — | store_true | 只保留双等位基因位点｜Keep only biallelic sites |
| `--remove-indels` | — | store_true | 移除插入缺失变异｜Remove indel variants |
| `--skip-validation` | `True` | store_true | 跳过输入验证以提高速度（默认开启）｜ Skip input validation for speed (default enabled) |
| `--force-validation` | — | store_true | 强制执行输入验证｜Force input validation |
| `--verbose, -v` | — | store_true | 显示详细信息｜Show verbose information |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies { #dependencies }

- bcftools（必需，核心引擎）
- tabix（可选，.gz 输入建索引时用；推荐 conda 环境 `align`）
- bgzip（可选）
- Python 3 标准库

## 常见问题 | FAQ { #faq }

**Q1：会断点续传吗？**
不会。每次运行重新筛选并覆盖输出。

**Q2：`--min-depth`/`--max-depth` 为什么没效果？**
请求深度过滤会切换到较慢的 Python 模式，但当前代码里深度过滤逻辑尚未实现（源码注释注明省略），等于没有过滤。需要按深度过滤请用 bcftools 的 DP 表达式自行处理。

**Q3：`--convert-format` 能转 PLINK 格式吗？**
当前主流程不再调用 PLINK，该参数与 `--plink-path`/`--allow-extra-chr` 只是保留接口，实际不生效。

**Q4：`--keep-ids`/`--remove-ids` 有用吗？**
当前主流程未把它们传入筛选核心，实际不生效，请改用样本或质量条件筛选。

**Q5：.gz 输入报错或很慢？**
工具会自动为 .gz 生成 tabix 索引；若系统缺 tabix，性能会下降，建议安装（conda 环境 `align`）。

**Q6：MAF 过滤是怎么算的？**
用的是双侧判据 `AF > maf 且 AF < 1-maf`，即同时排除过稀有和过主流的等位基因。
