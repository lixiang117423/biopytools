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

## 依赖 | Dependencies

- **QIIME2**: 微生物组分析框架(含 cutadapt/DADA2/vsearch/classify-sklearn 等)
- **SILVA / UNITE**: 分类参考库(用于训练分类器)

## 相关链接 | References

- [项目主页](https://github.com/lixiang117423/biopytools)
