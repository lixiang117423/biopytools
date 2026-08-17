# HiFiasm 基因组组装 | HiFiasm Genome Assembly

**基于 PacBio HiFi reads 的基因组组装自动化流程 | Automated genome assembly pipeline based on PacBio HiFi reads**

## 功能概述 | Overview

HiFiasm 是一款主流的高质量基因组组装软件，擅长利用 PacBio HiFi 高准确度长读长数据进行单倍型分相组装。本模块对 hifiasm 进行了完整封装，支持纯 HiFi 组装、HiFi + ONT 整合组装、HiFi + Hi-C 染色体级别挂载等多种组合方案。

流程覆盖从输入数据校验、组装、purge_dup 冗余清理、GFA 转 FASTA，到 BUSCO/QUAST 质量评估、单倍型差异分析的完整链路，并自动整理日志、配置文件和结果报告。适用于二倍体、多倍体物种的参考基因组或泛基因组构建。

## 快速开始 | Quick Start

```bash
# 基本用法：纯 HiFi 组装
biopytools hifiasm -i sample.hifi.fq.gz -o ./hifiasm_results -p sample

# HiFi + Hi-C 染色体级别组装
biopytools hifiasm -i sample.hifi.fq.gz \
    --hi-c-1 hic_R1.fq.gz --hi-c-2 hic_R2.fq.gz \
    -o ./hifiasm_hic -p sample --hg-size 3g

# HiFi + ONT 整合组装并跳过评估
biopytools hifiasm -i sample.hifi.fq.gz --ont-reads ont.fq.gz \
    --skip-busco --skip-quast -t 32
```

## 参数说明 | Parameters

### 必需参数 | Required

| 参数 | 描述 |
|------|------|
| `-i, --input-reads` | HiFi 测序数据文件（FASTQ，可压缩）|

### 常用可选参数 | Common Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-o, --output-dir` | `./hifiasm_output` | 输出目录 |
| `-p, --prefix` | `sample` | 输出文件前缀 |
| `-t, --threads` | `12` | 线程数 |
| `--hg-size` | `auto` | 基因组大小估计（如 `1.4g`、`500m`）|
| `-l, --purge-level` | `3` | Purge 级别（0-3）|
| `--purge-max` | `65` | 最大 purge 覆盖度 |
| `-s, --similarity-threshold` | `0.75` | 单倍型相似性阈值 |
| `--ont-reads` | - | ONT 长读长数据文件（整合组装）|
| `--hi-c-1 / --hi-c-2` | - | Hi-C 双端数据文件（需同时提供）|
| `--extra-hifiasm-args` | `''` | 额外透传给 hifiasm 的参数 |
| `--skip-busco / --skip-quast` | 关 | 跳过对应质量评估 |
| `--busco-lineage` | `auto` | BUSCO 谱系数据集 |
| `--reference-genome` | - | 参考基因组（用于 QUAST 比对）|
| `--min-contig-length` | `1000` | 最小 contig 长度过滤 |
| `--assembly-type` | `auto` | 组装类型（auto/diploid/triploid/polyploid）|
| `--output-formats` | `both` | 输出格式（fasta/gfa/both，可多选）|
| `--memory` | `100` | 内存限制（GB）|
| `--max-runtime` | `48` | 最大运行时间（小时）|
| `--resume` | 关 | 恢复中断的分析 |
| `--dry-run` | 关 | 只打印命令不执行 |
| `--config-file` | - | YAML 配置文件路径 |

（运行 `biopytools hifiasm -h` 查看完整参数列表，包括各外部工具路径）

## 输出 | Output

```
hifiasm_output/
├── 01_assembly/        # hifiasm 原始 GFA 组装结果
├── 02_fasta/           # 转换后的单倍型 FASTA 文件
├── 03_purge/           # purge_dups 清理结果
├── 04_busco/           # BUSCO 完整性评估
├── 05_quast/           # QUAST 组装统计
├── 06_haplotype/       # 单倍型差异分析（可选）
├── logs/               # 运行日志
└── config.yaml         # 本次运行配置快照
```

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input-reads, -i` | 必填 |  | HiFi测序数据文件｜Input HiFi sequencing data file |
| `--output-dir, -o` | `./hifiasm_output` | Path | 输出目录｜Output directory |
| `--prefix, -p` | `sample` | str | 输出文件前缀｜Output file prefix |
| `--threads, -t` | `12` | int | 线程数｜Number of threads |
| `--hg-size` | `auto` | str | 基因组大小估计(如1.4g, 2.1g)｜Genome size estimation (e.g., 1.4g, 2.1g) |
| `--purge-level, -l` | `3` | int | Purge级别(0-3)｜Purge level (0-3) |
| `--purge-max` | `65` | int | 最大purge覆盖度｜Maximum purge coverage |
| `--similarity-threshold, -s` | `0.75` | float | 相似性阈值｜Similarity threshold |
| `--ont-reads` | — |  | ONT长读长数据文件｜ONT long-read data file |
| `--hi-c-1` | — |  | Hi-C第一端数据文件｜Hi-C first-end data file |
| `--hi-c-2` | — |  | Hi-C第二端数据文件｜Hi-C second-end data file |
| `--extra-hifiasm-args` | `` | str | 额外的HiFiasm参数｜Additional HiFiasm arguments |
| `--skip-busco` | — |  | 跳过BUSCO质量评估｜Skip BUSCO quality assessment |
| `--busco-lineage` | `auto` | str | BUSCO谱系数据集｜BUSCO lineage dataset |
| `--busco-mode` | `genome` | genome/proteins/transcriptome | BUSCO评估模式｜BUSCO assessment mode |
| `--skip-quast` | — |  | 跳过QUAST质量评估｜Skip QUAST quality assessment |
| `--reference-genome` | — |  | 参考基因组文件(用于QUAST)｜Reference genome file (for QUAST) |
| `--analyze-haplotypes` | — |  | 分析单倍型差异｜Analyze haplotype differences |
| `--min-contig-length` | `1000` | int | 最小contig长度过滤｜Minimum contig length filter |
| `--generate-plots` | — |  | 生成可视化图表｜Generate visualization plots |
| `--assembly-type` | `auto` | auto/diploid/triploid/polyploid | 组装类型｜Assembly type |
| `--keep-intermediate` | — |  | 保留中间文件｜Keep intermediate files |
| `--compress-output` | — |  | 压缩输出文件｜Compress output files |
| `--output-formats` | `['both']` | fasta/gfa/both | 输出格式选择｜Output format selection |
| `--memory` | `100` | int | 内存大小(GB)｜Memory size (GB) |
| `--tmp-dir` | — | Path | 临时目录(默认 output_dir/tmp)｜Temporary directory (defaults to output_dir/tmp) |
| `--max-runtime` | `48` | int | 最大运行时间(小时)｜Maximum runtime (hours) |
| `--resume` | — |  | 恢复中断的分析｜Resume interrupted analysis |
| `--hifiasm-path` | `hifiasm` | str | HiFiasm软件路径｜HiFiasm software path |
| `--busco-path` | `busco` | str | BUSCO软件路径｜BUSCO software path |
| `--quast-path` | `quast` | str | QUAST软件路径｜QUAST software path |
| `--python-path` | `python3` | str | Python解释器路径｜Python interpreter path |
| `--samtools-path` | `samtools` | str | Samtools软件路径｜Samtools software path |
| `--busco-db-path` | — | Path | BUSCO数据库路径｜BUSCO database path |
| `--busco-download-path` | — | Path | BUSCO数据集下载路径｜BUSCO dataset download path |
| `--debug` | — |  | 启用调试模式｜Enable debug mode |
| `--verbose, -v` | — |  | 详细输出模式｜Verbose output mode |
| `--log-level` | `INFO` | DEBUG/INFO/WARNING/ERROR | 日志级别｜Log level |
| `--config-file` | — |  | 配置文件路径｜Configuration file path |
| `--dry-run` | — |  | 试运行模式(不执行实际命令)｜Dry run mode (do not execute actual commands) |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input-reads` | 必填 |  | 输入HiFi测序数据文件｜Input HiFi sequencing data file |
| `-o, --output-dir` | `./hifiasm_output` |  | 输出目录｜Output directory |
| `-p, --prefix` | `sample` |  | 输出文件前缀｜Output file prefix |
| `-t, --threads` | `32` | int | 线程数｜Number of threads |
| `--hg-size` | `auto` |  | 基因组大小估计(如1.4g, 2.1g)｜Genome size estimation (e.g., 1.4g, 2.1g) |
| `-l, --purge-level` | `3` | int | purge级别(0-3)｜Purge level (0-3) |
| `--purge-max` | `65` | int | 最大purge覆盖度｜Maximum purge coverage |
| `-s, --similarity-threshold` | `0.75` | float | 相似性阈值｜Similarity threshold |
| `--ont-reads` | — |  | ONT长读长数据文件｜ONT long-read data file |
| `--hi-c-1` | — |  | Hi-C第一端数据文件｜Hi-C first-end data file |
| `--hi-c-2` | — |  | Hi-C第二端数据文件｜Hi-C second-end data file |
| `--extra-hifiasm-args` | `` |  | 额外的HiFiasm参数｜Additional HiFiasm arguments |
| `--skip-busco` | — | store_true | 跳过BUSCO质量评估｜Skip BUSCO quality assessment |
| `--busco-lineage` | `auto` |  | BUSCO谱系数据集(如embryophyta_odb10)｜BUSCO lineage dataset |
| `--busco-mode` | `genome` | genome/proteins/transcriptome | BUSCO评估模式｜BUSCO assessment mode |
| `--skip-quast` | — | store_true | 跳过QUAST质量评估｜Skip QUAST quality assessment |
| `--reference-genome` | — |  | 参考基因组文件(用于QUAST)｜Reference genome file (for QUAST) |
| `--analyze-haplotypes` | — | store_true | 分析单倍型差异｜Analyze haplotype differences |
| `--min-contig-length` | `1000` | int | 最小contig长度过滤｜Minimum contig length filter |
| `--generate-plots` | — | store_true | 生成可视化图表｜Generate visualization plots |
| `--assembly-type` | `auto` | auto/diploid/triploid/polyploid | 组装类型｜Assembly type |
| `--keep-intermediate` | — | store_true | 保留中间文件｜Keep intermediate files |
| `--compress-output` | — | store_true | 压缩输出文件｜Compress output files |
| `--output-formats` | `['both']` | fasta/gfa/both | 输出格式选择｜Output format selection |
| `--memory` | `64` | int | 内存大小(GB)｜Memory size (GB) |
| `--tmp-dir` | — |  | 临时目录(默认 output_dir/tmp)｜Temporary directory (defaults to output_dir/tmp) |
| `--max-runtime` | `48` | int | 最大运行时间(小时)｜Maximum runtime (hours) |
| `--resume` | — | store_true | 恢复中断的分析｜Resume interrupted analysis |
| `--hifiasm-path` | `hifiasm` |  | HiFiasm软件路径｜HiFiasm software path |
| `--busco-path` | `busco` |  | BUSCO软件路径｜BUSCO software path |
| `--quast-path` | `quast` |  | QUAST软件路径｜QUAST software path |
| `--python-path` | `python3` |  | Python解释器路径｜Python interpreter path |
| `--samtools-path` | `samtools` |  | Samtools软件路径｜Samtools software path |
| `--busco-db-path` | — |  | BUSCO数据库路径｜BUSCO database path |
| `--busco-download-path` | — |  | BUSCO数据集下载路径｜BUSCO dataset download path |
| `--debug` | — | store_true | 启用调试模式｜Enable debug mode |
| `--verbose, -v` | `0` | count | 详细输出模式(-v, -vv, -vvv)｜Verbose output mode |
| `--log-level` | `INFO` | DEBUG/INFO/WARNING/ERROR | 日志级别｜Log level |
| `--config-file` | — |  | 配置文件路径｜Configuration file path |
| `--dry-run` | — | store_true | 试运行模式(不执行实际命令)｜Dry run mode (do not execute actual commands) |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- hifiasm（默认在 PATH 中，可通过 `--hifiasm-path` 指定）
- samtools
- BUSCO（可选，`--busco-path`）
- QUAST（可选，`--quast-path`）
- purge_dups（用于冗余序列清理）

## 引用 | Citation

- Cheng, H. et al. HiFiasm: haplotype-resolved de novo assembly using PacBio HiFi reads. *Nature Methods* 18, 170-175 (2021).
- Cheng, H. et al. Haplotype-resolved de novo assembly using phased assembly graphs with hifiasm. *Nature Methods* 19, 632-635 (2022).

## 相关链接 | References

- [项目主页](https://github.com/lixiang117423/biopytools)
- [hifiasm 官方仓库](https://github.com/chhylp123/hifiasm)
