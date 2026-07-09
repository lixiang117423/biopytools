# BRAKER3 基因组注释分析模块

**专业的真核生物基因组基因结构注释工具 | Professional Eukaryotic Genome Gene Structure Annotation Tool**

## 功能概述 | Overview

BRAKER3 基因组注释分析模块是一个基于 BRAKER3、GeneMark-ETP 和 AUGUSTUS 的真核生物基因组基因结构注释工具，提供从重复序列屏蔽、转录本比对到最终基因注释的完整流程。支持多种证据类型（蛋白质同源性、三代全长转录本、二代RNA-seq），特别适合卵菌等真菌界物种的基因组注释。

## 主要特性 | Key Features

- **完整注释流程**: 重复序列屏蔽→转录本比对→BRAKER3注释→TSEBRA整合五步骤自动化
- **多证据整合**: 支持蛋白质同源性、三代全长转录本、二代RNA-seq多种证据类型
- **Singularity容器**: 使用Singularity镜像运行BRAKER3，包含完整的GeneMark-ETP环境
- **灵活步骤控制**: 支持单独运行任意步骤或完整流程
- **卵菌优化**: 针对卵菌优化的真菌模式和参数配置
- **断点续传**: 自动检测已完成步骤，支持中断后继续运行
- **详细日志**: 完整的处理过程日志和错误追踪
- **重复序列误伤修复**: 智能过滤 repeat 库(剔除误判的基因家族) + 证据驱动还原(蛋白/转录组)，解决效应子等多拷贝基因家族被 RepeatModeler 误判为重复序列的问题

## 快速开始 | Quick Start

### 基本用法 | Basic Usage

```bash
# 使用蛋白质数据进行注释
biopytools braker \
    --genome genome.fa \
    --species my_oomycete \
    --prot_seq oomycete_proteins.fa \
    --threads 12 \
    --fungus

# 结合多种证据类型
biopytools braker \
    --genome genome.fa \
    --species my_oomycete \
    --prot_seq proteins.fa \
    --isoseq isoseq.fa \
    --rnaseq_dirs /path/to/rnaseq1,/path/to/rnaseq2 \
    --rnaseq_sets set1,set2 \
    --threads 12 \
    --fungus

# 使用已屏蔽的基因组跳过重复序列步骤
biopytools braker \
    --genome masked_genome.fa \
    --species my_oomycete \
    --bam rnaseq.bam \
    --skip_repeat \
    --fungus
```

### 高级用法 | Advanced Usage

```bash
# 自定义输出目录和Singularity镜像
biopytools braker \
    --genome genome.fa \
    --species my_oomycete \
    --prot_seq proteins.fa \
    --output_dir /path/to/output \
    --singularity_image /path/to/braker3.sif \
    --threads 24 \
    --fungus

# 启用UTR预测和BUSCO评估
biopytools braker \
    --genome genome.fa \
    --species my_oomycete \
    --prot_seq proteins.fa \
    --busco_lineage stramenopiles_odb10 \
    --utr \
    --fungus

# 使用已有训练参数
biopytools braker \
    --genome genome.fa \
    --species my_oomycete \
    --bam rnaseq.bam \
    --training_genes existing_genes.gtf \
    --use_existing \
    --fungus
```

## 参数说明 | Parameters

### 必需参数 | Required Parameters

| 参数 | 描述 | 示例 |
|------|------|------|
| `-g, --genome` | 基因组FASTA文件 | `-g genome.fa` |
| `-s, --species` | 物种名称(用于BRAKER输出命名) | `-s my_oomycete` |

### 输入数据 | Input Data

| 参数 | 描述 | 默认值 |
|------|------|--------|
| `-p, --prot_seq` | 近缘物种蛋白质序列文件 | 无 |
| `-l, --isoseq` | 三代全长转录本文件 | 无 |
| `--rnaseq_dirs` | 二代RNA-seq目录列表(逗号分隔) | 无 |
| `--rnaseq_sets` | RNA-seq集合ID列表(逗号分隔) | 无 |

### 输出配置 | Output Configuration

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-o, --output_dir` | `./braker_output` | 输出目录路径 |

### 流程参数 | Pipeline Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-t, --threads` | `12` | 线程数 |
| `--fungus` | `False` | 启用真菌模式(适合卵菌) |

### Singularity配置 | Singularity Configuration

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--singularity_image` | `/share/.../braker3_devel.sif` | Singularity镜像路径 |
| `--no_singularity` | `False` | 不使用Singularity镜像 |

### 步骤控制 | Step Control

| 参数 | 描述 |
|------|------|
| `--skip_repeat` | 跳过重复序列屏蔽步骤 |
| `--skip_long_reads` | 跳过三代转录本处理步骤 |
| `--skip_short_reads` | 跳过二代RNA-seq处理步骤 |
| `--skip_tsebra` | 跳过TSEBRA整合步骤 |

### BRAKER3特定参数 | BRAKER3 Specific Parameters

| 参数 | 描述 |
|------|------|
| `--busco_lineage` | BUSCO谱系(用于质量评估) |
| `--utr` | 预测UTR区域 |
| `--training_genes` | 训练基因集文件 |
| `--use_existing` | 使用已有物种参数 |

### 重复序列误伤修复参数 | Repeat Mis-mask Repair Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--skip_repeat_filter` | `False`(开启) | 跳过 repeat 库过滤(方案1) |
| `--skip_rescue` | `False`(开启) | 跳过证据驱动还原(方案2) |
| `--pfam_db` | `~/database/eggnog/pfam/Pfam-A.hmm` | Pfam-A HMM 库(TE domain 扫描) |
| `--te_domain_evalue` | `1e-5` | TE domain hmmscan E-value 阈值 |
| `--filter_min_orf_len` | `30` | 过滤用最小 ORF 长度(aa) |
| `--rescue_min_cds_len` | `100` | rescue 蛋白证据最小覆盖长度(bp) |
| `--rescue_min_identity` | `70` | rescue 蛋白最小 identity(%) |
| `--rescue_min_depth` | `5` | rescue RNA-seq 最小覆盖度 |

## 流程说明 | Pipeline Details

### 步骤1: 重复序列屏蔽 | Step 1: Repeat Masking

使用RepeatModeler构建重复序列库，**智能过滤误判的基因家族**(方案1)，然后用RepeatMasker进行软屏蔽。

```bash
RepeatModeler → [过滤: 六框翻译+hmmscan Pfam(TE域排除)+RxLR motif] → RepeatMasker (soft masking) → genome.masked.fa
```

**repeat 库过滤(方案1)**: RepeatModeler 的 de novo 建库会把多拷贝基因家族(如效应子)误判为重复序列。过滤步骤对每个 consensus 六框翻译后用 hmmscan 扫 Pfam: 命中 TE 特征域(逆转录酶/整合酶/转座酶等)→保留为真 TE; 命中非 TE 蛋白域或 RxLR/EER 效应子 motif →剔除(基因家族)。详见 [重复序列误伤修复](#重复序列误伤修复-repeat-mis-mask-repair)。

### 步骤2: 三代全长转录本处理 | Step 2: Long-read Processing

使用minimap2比对三代全长转录本(Iso-Seq/Nanopore)。

```bash
minimap2 -ax splice -uf -C5 genome.masked.fa isoseq.fa > isoseq.sam
```

### 步骤3: 二代RNA-seq处理 | Step 3: Short-read Processing

使用HISAT2比对二代RNA-seq数据。

```bash
hisat2-build genome.masked.fa index
hisat2 -x index -1 R1.fq -2 R2.fq --rna-strandness RF | samtools view -bS - > rnaseq.bam
```

### 步骤3.5: 证据驱动还原 | Step 3.5: Evidence-based Rescue

在所有证据(蛋白/转录本)就绪后、BRAKER 运行前，用证据把**仍被误 mask 的真基因区还原**(方案2)。对 masked genome 的每段连续小写区，检查是否有蛋白比对(miniprot)或 RNA-seq 覆盖支持，满足则 unmask(小写转大写)，确保 BRAKER 单轮即可注释这些基因。

```bash
miniprot(蛋白→masked genome) + samtools depth(RNA-seq) → 落在 masked 区且达阈值的段 → unmask → refined_masked.fa
```

### 步骤4: BRAKER3注释 | Step 4: BRAKER3 Annotation

在Singularity容器中运行BRAKER3。

```bash
singularity exec braker3_devel.sif perl /opt/BRAKER/scripts/braker.pl \
    --genome=genome.masked.fa \
    --species=my_oomycete \
    --bam=rnaseq.bam \
    --prot_seq=proteins.fa \
    --fungus \
    --softmasking
```

### 步骤5: TSEBRA整合 | Step 5: TSEBRA Merging

使用TSEBRA整合多个BRAKER运行结果。

```bash
singularity exec braker3_devel.sif python /opt/TSEBRA/bin/tsebra.py \
    -g braker_rnaseq.gtf,braker_protein.gtf \
    -o tsebra_combined.gtf
```

## 重复序列误伤修复 | Repeat Mis-mask Repair

### 问题背景 | Background

RepeatModeler(de novo 建库)会把**多拷贝基因家族**(如疫霉 RXLR/CRN 效应子、抗病基因)误判为重复序列家族，RepeatMasker 据此将其 soft-mask，导致 BRAKER 注释完全丢失这些基因。这是效应子注释的经典坑(效应子常成簇扩增 + 邻近转座子，触发 RepeatModeler 误判)。

### 两层修复 | Two-layer Repair

**方案1 · repeat 库过滤(默认开启，Step1 内)**: RepeatModeler 建库后，自动识别并剔除"实为基因家族"的 consensus:
- 六框翻译 consensus → 蛋白 ORF
- hmmscan 扫 Pfam-A: 命中 TE 特征域(逆转录酶/整合酶/转座酶/GAG/RNaseH 等)→ **保留**为真 TE; 命中非 TE 蛋白域 → **剔除**(基因家族)
- RxLR/EER 效应子 motif 命中 → **剔除**(效应子家族)
- 无明确证据 → 保守保留(避免误剔真 TE)
- 输出 `01_repeat_masking/filtered_consensi.fa.classified` + `filtered_families.tsv`(剔除报告)

**方案2 · 证据驱动还原(默认开启，Step3.5)**: RepeatMasker 后、BRAKER 前，用用户传入的证据(蛋白和/或转录组)还原仍被误 mask 的真基因区:
- 蛋白证据(prot_seq): miniprot 蛋白→masked genome，落在 masked 区且 identity≥阈值 → unmask
- 转录组证据(RNA-seq/Iso-Seq BAM): samtools depth 覆盖≥阈值 → unmask
- 输出 `01_repeat_masking/refined_masked.fa` + `rescued_regions.bed`

两层互补: 方案1 从源头不让基因家族被 mask; 方案2 用实证兜底补漏。均在 BRAKER 之前完成，**单轮 BRAKER**。

### 案例 | Example

大豆疫霉菌 6497 注释中，效应子 avr1a 的三个拷贝(Chr04)缺失。排查发现 RepeatModeler 把 avr1a 家族误学成 `rnd-1_family-619`，三段全被 soft-mask。启用 repeat 库过滤后，filter 通过 RxLR motif 判定剔除 family-619，RepeatMasker 不再 mask avr1a，注释恢复。

### 关闭修复 | Disabling

若需复现旧行为(不过滤、不还原):

```bash
biopytools braker -g genome.fa -s species -p prot.fa --skip_repeat_filter --skip_rescue
```

## 输出文件 | Output Files

```
braker_output/
├── 01.repeat_masking/
│   ├── genome.masked.fa        # 屏蔽后的基因组
│   └── consensi.fa.classified   # 重复序列库
├── 02.long_reads/
│   └── isoseq.sam               # 三代转录本比对结果
├── 03.short_reads/
│   ├── rnaseq.bam               # RNA-seq比对结果
│   └── rnaseq.bam.bai           # BAM索引
├── 04.braker_annotation/
│   ├── braker.gtf               # 最终基因注释(GTF格式)
│   ├── braker.aa                # 蛋白质序列
│   └── hintsfile.gff            # Hints文件
├── 05.tsebra_merge/
│   └── tsebra_combined.gtf      # TSEBRA整合结果
└── logs/
    └── braker_pipeline.log      # 完整运行日志
```

## 技术要求 | Requirements

### 软件依赖 | Software Dependencies

#### Singularity镜像中包含 | Included in Singularity Image
- BRAKER3 v3.0.8
- GeneMark-ETP
- AUGUSTUS
- TSEBRA
- SAMtools

#### 宿主机工具 | Host Machine Tools
- RepeatModeler v2.0.7
- RepeatMasker
- minimap2
- HISAT2
- Singularity v3.8.7

### 系统要求 | System Requirements

- **操作系统**: Linux (CentOS 7+/RHEL 7+)
- **内存**: 建议 256GB+ (取决于基因组大小)
- **存储**: 建议 500GB+ (取决于数据量和基因组大小)
- **CPU**: 建议 32核心+

## 常见问题 | FAQ

### Q1: 为什么必须提供至少一种证据数据？

A: BRAKER3需要外部证据来训练基因预测模型。可以选择蛋白质序列、三代转录本或二代RNA-seq中的至少一种。

### Q2: 真菌模式适合哪些物种？

A: 真菌模式适合具有类似剪接位点模式的物种，包括真菌和卵菌（Oomycetes）。

### Q3: 如何选择证据数据类型？

A: 优先级推荐：
1. **蛋白质序列**: 近缘物种的高质量蛋白质注释
2. **三代全长转录本**: 提供最准确的转录本证据
3. **二代RNA-seq**: 提供表达量信息

### Q4: Singularity镜像在哪里？

A: 默认路径为 `/share/org/YZWL/yzwl_lixg/software/singularity/braker3_devel.sif`，包含完整的BRAKER3和GeneMark-ETP环境。

### Q5: 如何中断后继续运行？

A: 流程支持断点续传，重新运行相同命令会自动跳过已完成的步骤。

### Q6: 如何评估注释质量？

A: 可以使用BUSCO评估注释完整性：

```bash
busco -i braker.aa -l stramenopiles_odb10 -m protein
```

### Q7: 效应子/多基因家族的基因注释缺失怎么办？

A: 很可能被 RepeatModeler 误判为重复序列并 mask。检查 `01_repeat_masking/*.out` 中目的基因区域是否被某个 `rnd-X_family-Y` 覆盖。若是，确认 repeat 库过滤(方案1)和证据还原(方案2)已开启(默认开启)，重新运行即可自动剔除误判家族并还原被 mask 的真基因区。详见 [重复序列误伤修复](#重复序列误伤修复-repeat-mis-mask-repair)。

## 参考资料 | References

1. Bruna T, et al. (2021) BRAKER2: automatic eukaryotic genome annotation with GeneMark-EP+ and AUGUSTUS. Bioinformatics, 37(12): 1706-1707.

2. Hoff KJ, et al. (2016) BRAKER1: unsupervised RNA-Seq-based genome annotation with GeneMark-ET and AUGUSTUS. Bioinformatics, 32(5), 767-769.

3. Stanke M, et al. (2006) Gene prediction with a hidden Markov model and a dual intron model. Bioinformatics, 22(13), i47-i55.

## 版本历史 | Version History

| 版本 | 日期 | 主要变更 |
|------|------|----------|
| 1.1.0 | 2026-07-09 | 新增重复序列误伤修复: repeat 库智能过滤(方案1，TE域排除+RxLR motif) + 证据驱动还原(方案2，蛋白/转录组 unmask)，解决多拷贝基因家族被误 mask 问题 |
| 1.0.0 | 2026-02-02 | 初始版本 |