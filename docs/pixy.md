# Pixy群体遗传统计 | Pixy Population Genetics Statistics

**基于pixy计算pi/dxy/fst/Watterson theta/Tajima's D, 无偏且支持不变位点 | Unbiased population genetics statistics supporting invariant sites**

## 功能概述 | Overview

pixy 模块封装了 pixy 工具, 用于计算群体遗传学统计量。相比传统工具(如VCFtools), pixy 通过包含不变位点(invariant sites)避免"除以零"导致的fst高估问题, 是当前计算pi、dxy、fst的金标准工具。支持滑窗和全基因组两种模式, 适用于群体遗传结构分析、群体分化评估、选择信号检测等场景。

## 快速开始 | Quick Start

```bash
# 滑窗分析(100kb窗口)
biopytools pixy -i variants.vcf.gz -p populations.txt -o pixy_output -w 100000

# 指定统计量
biopytools pixy -i variants.vcf.gz -p pops.txt -o out --stats pi,fst,dxy
```

## 参数说明 | Parameters

### 必需参数 | Required

| 参数 | 描述 |
|------|------|
| `-i, --vcf-file` | VCF文件(需bgzip压缩并建立tabix索引) |
| `-p, --pop-file` | 群体文件(两列:样本ID 群体名) |

### 常用可选参数 | Common Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-o, --output-dir` | `./pixy_output` | 输出目录 |
| `--stats` | `pi,dxy,fst,watterson,tajima` | 要计算的统计量,逗号分隔 |
| `-w, --window-size` | `None` | 窗口大小bp(不设则全基因组) |
| `-b, --bed-file` | `None` | BED文件定义窗口 |
| `-s, --sites-file` | `None` | 位点文件(只计算特定位点) |
| `--min-samples` | `0` | 每个群体最小样本数 |
| `--max-missing` | `1.0` | 最大缺失率 |
| `--min-maf` | `0.0` | 最小等位基因频率 |
| `-c, --chromosomes` | `None` | 指定染色体列表,逗号分隔 |
| `-t, --threads` | `12` | 线程数 |
| `--pixy-path` | `pixy` | pixy可执行文件路径 |
| `--conda-env` | `~/miniforge3/envs/pixy_v.2.0.0` | conda环境路径 |
| `--bypass-invariant-check` | `False` | 强制绕过不变位点检查(默认自动检测) |
| `--keep-intermediate` | `False` | 保留中间文件 |

(运行 `biopytools pixy -h` 查看完整参数列表)

> **注意**: pixy必须指定`-w`(窗口)、`-b`(BED文件)或`-s`(位点文件)之一。

## 输出 | Output

- `pixy_{stat}_{pop}_{window}.txt`: 每个统计量的滑窗/全基因组结果
- `pixy.log`: 运行日志

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --vcf-file` | 必填 |  | VCF文件路径（需用bgzip压缩并建立tabix索引）｜VCF file path (must be bgzip-compressed and tabix-indexed) |
| `-p, --pop-file` | 必填 |  | 群体文件路径（两列：样本ID 群体名）｜Population file path (two columns: sample_id population_name) |
| `-o, --output-dir` | `./pixy_output` | Path | 输出目录｜Output directory |
| `--stats` | `pi,dxy,fst,watterson,tajima` |  | 要计算的统计量，逗号分隔｜Statistics to calculate, comma-separated |
| `--calc-pi` | — |  | 计算pi（核苷酸多样性）｜Calculate pi (nucleotide diversity) |
| `--calc-dxy` | — |  | 计算dxy（群体间核苷酸差异）｜Calculate dxy (nucleotide divergence) |
| `--calc-fst` | — |  | 计算fst（遗传分化系数）｜Calculate fst (genetic differentiation) |
| `--calc-watterson-theta` | — |  | 计算Watterson's theta｜Calculate Watterson's theta |
| `--calc-tajima-d` | — |  | 计算Tajima's D｜Calculate Tajima's D |
| `-w, --window-size` | — | int | 窗口大小bp（不设置则全基因组计算）｜Window size in bp (null for genome-wide) |
| `-b, --bed-file` | — |  | BED文件定义窗口｜BED file defining windows |
| `-s, --sites-file` | — |  | 位点文件（只计算特定位点）｜Sites file (calculate only specific sites) |
| `--min-samples` | `0` | int | 每个群体最小样本数｜Minimum samples per population |
| `--max-missing` | `1.0` | float | 最大缺失率｜Maximum missing rate |
| `--min-maf` | `0.0` | float | 最小等位基因频率｜Minor allele frequency |
| `--zscore-window` | — | int | Z-score过滤窗口大小｜Z-score filtering window size |
| `-c, --chromosomes` | — |  | 指定染色体列表，逗号分隔｜List of chromosomes, comma-separated |
| `--pixy-path` | `pixy` |  | pixy可执行文件路径｜pixy executable path |
| `--conda-env` | `~/miniforge3/envs/pixy_v.2.0.0` |  | conda环境路径｜conda environment path |
| `-t, --threads` | `12` | int | 线程数｜Number of threads |
| `--bypass-invariant-check` | — |  | 强制绕过不变位点检查（默认自动检测VCF并自动绕过）｜Force bypass invariant sites check (default: auto-detect VCF and bypass if needed) |
| `--keep-intermediate` | — |  | 保留中间文件｜Keep intermediate files |
| `-v, --verbose` | — |  | 详细输出｜Verbose output |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --vcf-file` | 必填 |  | 输入VCF文件（需用bgzip压缩并建立tabix索引）｜Input VCF file (must be bgzip-compressed and tabix-indexed) |
| `-p, --pop-file` | 必填 |  | 群体文件（两列：样本ID 群体名）｜Population file (two columns: sample_id population_name) |
| `-o, --output-dir` | 必填 |  | 输出目录｜Output directory |
| `--stats` | `pi,dxy,fst,watterson,tajima` |  | 要计算的统计量，逗号分隔（默认: pi,dxy,fst,watterson,tajima）｜Statistics to calculate, comma-separated (default: pi,dxy,fst,watterson,tajima) |
| `--calc-pi` | — | store_true | 计算pi（核苷酸多样性）｜Calculate pi (nucleotide diversity) |
| `--calc-dxy` | — | store_true | 计算dxy（群体间核苷酸差异）｜Calculate dxy (nucleotide divergence) |
| `--calc-fst` | — | store_true | 计算fst（遗传分化系数）｜Calculate fst (genetic differentiation) |
| `--calc-watterson-theta` | — | store_true | 计算Watterson's theta｜Calculate Watterson's theta |
| `--calc-tajima-d` | — | store_true | 计算Tajima's D｜Calculate Tajima's D |
| `-w, --window-size` | — | int | 窗口大小bp（不设置则全基因组计算）｜Window size in bp (null for genome-wide) |
| `-b, --bed-file` | — |  | BED文件定义窗口（自定义大小窗口）｜BED file defining windows (custom-sized windows) |
| `-s, --sites-file` | — |  | 位点文件（只计算特定位点）｜Sites file (calculate only specific sites) |
| `--min-samples` | `0` | int | 每个群体最小样本数（默认: 0=不限制）｜Minimum samples per population (default: 0=no limit) |
| `--max-missing` | `1.0` | float | 最大缺失率（默认: 1.0=不限制）｜Maximum missing rate (default: 1.0=no limit) |
| `--min-maf` | `0.0` | float | 最小等位基因频率（默认: 0.0=不限制）｜Minor allele frequency (default: 0.0=no limit) |
| `--zscore-window` | — | int | Z-score过滤窗口大小（不设置则不过滤）｜Z-score filtering window size (null=no filter) |
| `-c, --chromosomes` | — |  | 指定染色体列表，逗号分隔（不设置则全部）｜List of chromosomes, comma-separated (null for all) |
| `--pixy-path` | `pixy` |  | pixy可执行文件路径（默认: pixy）｜pixy executable path (default: pixy) |
| `--conda-env` | `~/miniforge3/envs/pixy_v.2.0.0` |  | conda环境路径｜conda environment path |
| `-t, --threads` | `12` | int | 线程数（默认: 12）｜Number of threads (default: 12) |
| `--bypass-invariant-check` | — | store_true | 强制绕过不变位点检查（默认自动检测VCF并自动绕过）｜Force bypass invariant sites check (default: auto-detect VCF and bypass if needed) |
| `--keep-intermediate` | — | store_true | 保留中间文件｜Keep intermediate files |
| `-v, --verbose` | — | store_true | 详细输出｜Verbose output |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- **pixy**: 无偏群体遗传统计工具 (https://github.com/ksamuk/pixy)
- **tabix/bgzip**: VCF索引 (htslib)
- **conda环境**: 推荐在pixy专用环境运行

## 引用 | Citation

- Korunes K.L., Samuk K. (2021) pixy: Unbiased estimation of nucleotide diversity within and between populations. Molecular Ecology Resources. 21(3):878-889.

## 相关链接 | References

- [项目主页](https://github.com/lixiang117423/biopytools)
