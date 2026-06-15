# NLR-Annotator NLR基因预测模块

**从CDS序列预测NLR（Nucleotide-binding Leucine-rich Repeat）抗病基因 | Predict NLR genes from CDS sequences**

## 功能概述 | Overview

NLR-Annotator模块封装了 [NLR-Annotator](https://github.com/stevenpbaxter/NLR-Annotator) Java工具，基于motif组合方法从CDS序列中识别NLR类抗病基因。模块支持单文件和目录批处理两种模式，自动清洗输出（添加表头、去重排序motif），并在批处理模式下生成多样本汇总文件。运行过程中自动跳过已完成的样本，支持断点续传。

## 快速开始 | Quick Start

```bash
# 单文件预测
biopytools nlr-annotator -i genome.cds.fa -o output_dir/

# 目录批处理（匹配 *.cds.fa 文件）
biopytools nlr-annotator -i cds_dir/ -o output_dir/ -t 16

# 同时输出GFF和BED
biopytools nlr-annotator -i genome.cds.fa -o output_dir/ --output-gff --output-bed
```

## 参数说明 | Parameters

### 必需参数 | Required

| 参数 | 描述 |
|------|------|
| `-i, --input` | 输入CDS FASTA文件或目录（目录模式匹配 `--sample-suffix`） |
| `-o, --output-dir` | 输出目录，默认 `./output` |

### 常用可选参数 | Common Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-t, --threads` | `12` | 线程数 |
| `--sample-suffix` | `*.cds.fa` | 目录模式下文件匹配模式 |
| `--output-gff` | `False` | 输出GFF注释文件 |
| `--output-bed` | `False` | 输出BED文件 |
| `--output-motifs` | `False` | 输出motifs BED文件 |
| `--output-alignment` | `False` | 输出motif比对FASTA |

（运行 `biopytools nlr-annotator -h` 查看完整参数列表，包括距离阈值等高级参数）

### 高级参数 | Advanced Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--jar-path` | 自动查找 | NLR-Annotator JAR文件路径 |
| `--mot-file` | 自动查找 | mot.txt配置文件路径 |
| `--store-file` | 自动查找 | store.txt配置文件路径 |
| `--num-seqs-per-thread` | `1000` | 每线程处理序列数 |
| `--distance-within-motif-combination` | `500` | motif组合内最大距离 |
| `--distance-for-elongating` | `2500` | motif延伸距离 |
| `--distance-between-motif-combinations` | `50000` | motif组合间最大距离 |

可通过环境变量 `NLR_ANNOTATOR_PATH` 自定义JAR路径。

## 输出 | Output

```
output_dir/
├── {sample}.nlr_annotator.tsv      # NLR预测结果（带表头）
├── {sample}.nlr_annotator.gff      # GFF格式（--output-gff）
├── {sample}.nlr_annotator.bed      # BED格式（--output-bed）
├── {sample}.nlr_annotator_motifs.bed    # motif BED（--output-motifs）
├── {sample}.nlr_annotator_alignment.fa  # motif比对（--output-alignment）
├── nlr_annotator_summary.tsv       # 多样本汇总（仅批处理模式）
└── 99_logs/                        # 日志目录
```

TSV文件列：`gene_id  nlr_id  type  start  end  strand  motifs`（motif已去重排序，去除`motif_`前缀）。

## 依赖 | Dependencies

- Java运行环境（JRE >= 1.8）
- NLR-Annotator JAR文件，默认查找路径：`~/software/NLR-Annotator/NLR-Annotator-v2.1b.jar`
- 配套文件 `mot.txt` 和 `store.txt`，默认在同目录下

## 引用 | Citation

如果您在研究中使用本工具，请引用：
- Steuernagel, B. et al. NLR-Annotator: Mining MLO and NLR Genes in Plant Genomes. (https://github.com/stevenpbaxter/NLR-Annotator)

## 相关链接 | References

- [NLR-Annotator GitHub](https://github.com/stevenpbaxter/NLR-Annotator)
- [项目主页](https://github.com/lixiang117423/biopytools)
