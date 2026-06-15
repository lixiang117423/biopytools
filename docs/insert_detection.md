# insert-detection | 插入序列位点检测

**基于 Bowtie2 + soft-clip read 信号，从重测序数据中定位外源插入序列（如 T-DNA、转座子）在宿主基因组中的插入位点 | Detect insertion sites of an exogenous insert sequence in a host genome from paired-end resequencing reads**

## 功能概述 | Overview

本模块面向功能基因组学中常见的"插入序列定位"问题：已知一段外源序列（例如 T-DNA 载体、转座子标签、报告基因），需要从携带该插入的个体重测序数据中找到它在参考基因组上的精确插入位置。

核心策略是：用 Bowtie2 将 read 分别比对到"参考基因组 + 插入序列"组合索引上，从比对结果中寻找一边唯一比对到插入序列、另一端 soft-clip 到参考基因组的嵌合 read（split-read），根据这些 read 的剪接点对候选插入位点进行打分聚类，最终输出候选插入位点、支持 read 数与得分。流程自动识别 FASTQ 目录中的配对样本（默认匹配 fastp 输出的 `_1.clean.fq.gz` / `_2.clean.fq.gz`）。

## 快速开始 | Quick Start

```bash
biopytools insert-detection \
    -i host_genome.fa \
    --insert tdna_vector.fa \
    --fastq-dir ./fastq_output/ \
    -o ./insert_detection_output/
```

## 参数说明 | Parameters

### 必需参数 | Required

| 参数 | 描述 |
|------|------|
| `-i, --genome` | 参考基因组 FASTA 文件 |
| `--insert` | 插入序列 FASTA 文件 |
| `--fastq-dir` | FASTQ 文件目录（自动识别配对样本） |

### 常用可选参数 | Common Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-o, --output-dir` | `./insert_detection_output` | 输出目录 |
| `-t, --threads` | `12` | 线程数 |
| `--min-clip` | `20` | 最小 soft-clip 长度 |
| `--min-support` | `5` | 最小支持 read 数 |
| `--score-threshold` | `1000` | 得分阈值 |
| `--bowtie2-path` | `bowtie2` | Bowtie2 可执行文件路径 |
| `--samtools-path` | `samtools` | samtools 可执行文件路径 |
| `--read1-suffix` | `_1.clean.fq.gz` | R1 文件后缀（含扩展名） |
| `--read2-suffix` | `_2.clean.fq.gz` | R2 文件后缀（含扩展名） |
| `--force` | 关闭 | 强制重新运行所有步骤 |
| `--verbose` | 关闭 | 显示详细日志 |
| `--quiet` | 关闭 | 仅显示错误 |

（运行 `biopytools insert-detection -h` 查看完整参数列表）

## 输出 | Output

```
insert_detection_output/
├── 04_results/
│   ├── insert_sites.tsv          # 候选插入位点（染色体、位置、链、支持 read 数、得分）
│   └── insert_summary.txt        # 运行总结
└── 99_logs/                      # 每个样本/步骤的日志
```

## 依赖 | Dependencies

- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)（必需，比对）
- [samtools](http://www.htslib.org/)（必需，BAM 排序与索引）
- Python：标准库

## 引用 | Citation

- Langmead B, Salzberg SL. Fast gapped-read alignment with Bowtie 2. *Nature Methods*, 2012, 9(4): 357-359. doi:10.1038/nmeth.1923
- Li H et al. The Sequence Alignment/Map format and SAMtools. *Bioinformatics*, 2009, 25(16): 2078-2079. doi:10.1093/bioinformatics/btp352

## 相关链接 | References

- [项目主页](https://github.com/lixiang117423/biopytools)
- [Bowtie2 手册](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)
- [samtools 文档](http://www.htslib.org/doc/samtools.html)
