# QIIME2微生物组分析 | QIIME2 Microbiome Analysis

**用QIIME2完成扩增子(16S/ITS)微生物组多样性分析全流程, 支持ASV/OTU | Full amplicon (16S/ITS) microbiome diversity pipeline via QIIME2, ASV or OTU**

## 功能概述 | Overview

qiime2 模块封装了 [QIIME2](https://qiime2.org/), 从双端 FASTQ 出发完成扩增子微生物组分析全流程: 引物切除(cutadapt)→去重/聚类(DADA2 ASV 或 vsearch OTU)→建系统发育树→多样性分析(α/β)→物种注释(classify-sklearn)。

- 支持 16S 与 ITS 扩增子(ITS 自动跳过建树)
- 分类器可传入预训练 `.qza`, 否则用 SILVA/UNITE 参考库自动训练并缓存
- 支持断点续传与 `--force`

## 快速开始 | Quick Start

```bash
# 16S 默认流程(ASV)
biopytools qiime2 -i raw_reads/ -o qiime2_output

# ITS, 跳过引物切除(数据已去引物)
biopytools qiime2 -i raw_reads/ --amplicon its --skip-cutadapt

# 指定截断长度(2x250 V3-V4 常用)
biopytools qiime2 -i raw_reads/ --trunc-len-f 220 --trunc-len-r 200

# 传入预训练分类器
biopytools qiime2 -i raw_reads/ --classifier ~/db/silva_classifier.qza
```

## 参数说明 | Parameters

### 必需参数 | Required

| 参数 | 描述 |
|------|------|
| `-i, --input-dir` | 双端 FASTQ 输入目录 |

### 常用可选参数 | Common Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-o, --output-dir` | `./qiime2_output` | 输出目录 |
| `--amplicon` | `16s` | 扩增子类型(16s/its) |
| `--method` | `asv` | 聚类方法: ASV(DADA2) 或 OTU(vsearch) |
| `--fwd-primer` | `341F` 序列 | 正向引物(IUPAC) |
| `--rev-primer` | `806R` 序列 | 反向引物(IUPAC) |
| `--trunc-len-f` / `--trunc-len-r` | `0` | R1/R2 截断长度(0=不截断) |
| `--sampling-depth` | `0` | 抽平深度(0=自动取第 10 百分位) |
| `--perc-identity` | `0.97` | OTU 聚类相似度 |
| `--confidence` | `0.7` | classify-sklearn 置信度 |
| `--classifier` | 自动训练 | 预训练分类器(.qza) |
| `--database-dir` | `~/database/qiime2` | 原始参考库目录(SILVA/UNITE) |
| `--qiime-path` | `~/miniforge3/envs/qiime_v.2024.10.1/bin/qiime` | qiime 可执行文件路径 |
| `-t, --threads` | `12` | 线程数 |
| `--skip-cutadapt` | `False` | 跳过引物切除(数据已去引物) |
| `--skip-phylogeny` | `False` | 跳过系统发育建树(ITS 自动跳过) |

(运行 `biopytools qiime2 -h` 查看完整参数列表)

## 输出 | Output

- 特征表(ASV/OTU)、代表性序列、系统发育树(`.qza`)
- α/β 多样性结果
- 物种注释结果
- 运行日志(99_logs)

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input-dir` | 必填 |  | 双端FASTQ输入目录｜Input directory of paired-end FASTQ |
| `-o, --output-dir` | `./qiime2_output` | Path | 输出目录｜Output directory |
| `--amplicon` | `16s` | 16s/its | 扩增子类型｜Amplicon type |
| `--method` | `asv` | asv/otu | 聚类方法ASV(DADA2)或OTU(vsearch)｜Method: ASV or OTU |
| `--fwd-primer` | `CCTACGGGNGGCWGCAG` |  | 正向引物序列(IUPAC)｜Forward primer |
| `--rev-primer` | `GACTACHVGGGTATCTAATCC` |  | 反向引物序列(IUPAC)｜Reverse primer |
| `--trunc-len-f` | `0` | int | R1截断长度(0=不截断)｜R1 truncation length (0=none) |
| `--trunc-len-r` | `0` | int | R2截断长度(0=不截断)｜R2 truncation length (0=none) |
| `--trim-left-f` | `0` | int | R1左侧裁剪｜R1 trim left |
| `--trim-left-r` | `0` | int | R2左侧裁剪｜R2 trim left |
| `--sampling-depth` | `0` | int | 抽平深度(0=自动)｜Rarefaction depth (0=auto) |
| `--perc-identity` | `0.97` | float | OTU聚类相似度｜OTU identity |
| `--confidence` | `0.7` | float | 分类置信度｜Classification confidence |
| `--min-length` | `50` | int | extract-reads最小长度｜extract-reads min length |
| `--max-length` | `0` | int | extract-reads最大长度(0=不限)｜extract-reads max length (0=none) |
| `-t, --threads` | `12` | int | 线程数｜Number of threads |
| `--validate-level` | `min` | min/max | tools import校验级别｜import validate level |
| `--classifier` | — |  | 预训练分类器(.qza),省略则自动训练｜Pre-trained classifier (.qza) |
| `--database-dir` | — |  | 原始参考库目录(SILVA/UNITE)｜Raw reference DB directory |
| `--qiime-path` | — |  | qiime可执行文件路径｜qiime executable path |
| `--classifier-cache-dir` | — |  | 分类器缓存目录｜Classifier cache directory |
| `--r1-suffix` | `_1.clean.fq.gz` |  | R1文件后缀｜R1 file suffix |
| `--r2-suffix` | `_2.clean.fq.gz` |  | R2文件后缀｜R2 file suffix |
| `--skip-cutadapt` | — |  | 跳过引物切除｜Skip primer trimming |
| `--skip-phylogeny` | — |  | 跳过系统发育建树(ITS自动跳过)｜Skip phylogeny (auto for ITS) |
| `--metadata-file` | — |  | 样品元数据TSV(可选)｜Sample metadata TSV (optional) |
| `-f, --force` | — |  | 覆盖已有输出｜Overwrite existing outputs |
| `-v, --verbose` | — |  | 详细输出｜Verbose output |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input-dir` | 必填 |  | 双端FASTQ输入目录｜Input directory of paired-end FASTQ |
| `-o, --output-dir` | `./qiime2_output` |  | 输出目录｜Output directory (default: ./qiime2_output) |
| `--amplicon` | `16s` | 16s/its | 扩增子类型｜Amplicon type (default: 16s) |
| `--method` | `asv` | asv/otu | 聚类方法｜Method: ASV(DADA2) or OTU(vsearch) (default: asv) |
| `--fwd-primer` | `CCTACGGGNGGCWGCAG` |  | 正向引物序列(IUPAC)｜Forward primer (default: 341F) |
| `--rev-primer` | `GACTACHVGGGTATCTAATCC` |  | 反向引物序列(IUPAC)｜Reverse primer (default: 806R) |
| `--trunc-len-f` | `0` | int | R1截断长度(0=不截断)｜R1 truncation length (0=none) |
| `--trunc-len-r` | `0` | int | R2截断长度(0=不截断)｜R2 truncation length (0=none) |
| `--trim-left-f` | `0` | int | R1左侧裁剪｜R1 trim left |
| `--trim-left-r` | `0` | int | R2左侧裁剪｜R2 trim left |
| `--sampling-depth` | `0` | int | 抽平深度(0=自动取第10百分位)｜Rarefaction depth (0=auto) |
| `--perc-identity` | `0.97` | float | OTU聚类相似度｜OTU identity (default: 0.97) |
| `--confidence` | `0.7` | float | classify-sklearn置信度｜classification confidence (default: 0.7) |
| `--min-length` | `50` | int | extract-reads最小长度｜extract-reads min length (default: 50) |
| `--max-length` | `0` | int | extract-reads最大长度(0=不限)｜extract-reads max length (0=none) |
| `-t, --threads` | `12` | int | 线程数｜Number of threads (default: 12) |
| `--validate-level` | `min` | min/max | tools import校验级别｜import validate level (default: min) |
| `--classifier` | — |  | 预训练分类器(.qza),省略则自动训练｜Pre-trained classifier (.qza), auto-train if omitted |
| `--database-dir` | — |  | 原始参考库目录(SILVA/UNITE)｜Raw reference DB directory (default: ~/database/qiime2) |
| `--qiime-path` | — |  | qiime可执行文件路径｜qiime executable path (default: ~/miniforge3/envs/qiime_v.2024.10.1/bin/qiime) |
| `--classifier-cache-dir` | — |  | 分类器缓存目录｜Classifier cache directory (default: <db>/classifier_cache) |
| `--r1-suffix` | `_1.clean.fq.gz` |  | R1文件后缀｜R1 file suffix (default: _1.clean.fq.gz) |
| `--r2-suffix` | `_2.clean.fq.gz` |  | R2文件后缀｜R2 file suffix (default: _2.clean.fq.gz) |
| `--skip-cutadapt` | — | store_true | 跳过引物切除(数据已去引物)｜Skip primer trimming |
| `--skip-phylogeny` | — | store_true | 跳过系统发育建树(ITS自动跳过)｜Skip phylogeny (auto for ITS) |
| `--metadata-file` | — |  | 样品元数据TSV(可选)｜Sample metadata TSV (optional) |
| `-f, --force` | — | store_true | 覆盖已有输出｜Overwrite existing outputs |
| `-v, --verbose` | — | store_true | 详细输出｜Verbose output |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- **QIIME2**: 微生物组分析框架(含 cutadapt/DADA2/vsearch/classify-sklearn 等)
- **SILVA / UNITE**: 分类参考库(用于训练分类器)

## 相关链接 | References

- [项目主页](https://github.com/lixiang117423/biopytools)
