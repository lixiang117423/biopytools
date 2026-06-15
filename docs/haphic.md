# 🧬 HapHiC基因组Scaffolding模块

**基于Hi-C数据的快速、参考基因组独立的等位基因感知scaffolding工具 | Fast, Reference-Independent, Allele-Aware Genome Scaffolding Tool Based on Hi-C Data**

## 📖 功能概述 | Overview

HapHiC是一个专业的基因组scaffolding工具，利用Hi-C数据将contigs或scaffolds组装成染色体级别的伪分子。它具有等位基因感知能力，无需参考基因组即可处理单倍型分相、二倍体和多倍体基因组组装。相比其他Hi-C scaffolding工具，HapHiC在处理嵌合contigs、折叠contigs和交换错误方面具有更高的容错性，并且执行速度极快，大多数基因组可在1小时内完成scaffolding。

### 🆕 Pipeline模式更新 | Pipeline Mode Update

现采用HapHiC原生Pipeline模式，实现更高效、更可靠的scaffolding流程：

- **🔄 断点续传**: 自动检测已完成步骤，支持中断恢复
- **⚡ 一键执行**: 单步完成所有scaffolding流程
- **📊 自动可视化**: 默认生成Hi-C接触图
- **🥤 Juicebox支持**: 自动生成Juicebox兼容文件

## ✨ 主要特性 | Key Features

- **🧬 等位基因感知**: 支持单倍型分相、二倍体和多倍体基因组组装
- **🔍 无需参考基因组**: 完全基于Hi-C数据，不依赖任何参考序列
- **⚡ 超快速度**: 比其他工具快10-100倍，内存使用效率高
- **🛡️ 高容错性**: 对嵌合contigs、折叠contigs和交换错误容错率高
- **🔧 智能校正**: 可选的contig校正功能，自动修复组装错误
- **🔄 Pipeline模式**: 原生四步自动化流程：聚类→重新分配→排序定向→构建scaffolds
- **📊 自动可视化**: 默认生成Hi-C接触图（PDF/PNG格式）
- **🥤 Juicebox集成**: 自动生成Juicebox兼容的.hic和.assembly文件
- **🔄 断点续传**: 智能检测进度，支持中断恢复
- **⚙️ 高度可配置**: 丰富的参数选项，适应不同数据质量

## 🚀 快速开始 | Quick Start

### 基本用法 | Basic Usage

```bash
# 使用BAM文件进行scaffolding
biopytools haphic -a assembly.fa -b hic.bam -c 24

# 使用FASTQ文件（自动执行BWA比对）
biopytools haphic -a assembly.fa -1 hic_R1.fastq.gz -2 hic_R2.fastq.gz -c 24

# 高性能配置示例
biopytools haphic -a assembly.fa -b hic.bam -c 24 \
    --threads 32 --processes 16 --correct-nrounds 2
```

### 断点续传使用 | Resume Mode Usage

```bash
# 断点续传模式（默认启用）
biopytools haphic -a assembly.fa -b hic.bam -c 24
# 如果中断后再次运行，自动跳过已完成步骤

# 强制重新运行所有步骤
biopytools haphic -a assembly.fa -b hic.bam -c 24 --force-rerun
```

### 高级用法 | Advanced Usage

```python
# Python模块方式使用
from biopytools.haphic.main import HapHiCProcessor

processor = HapHiCProcessor(
    asm_file="assembly.fa",
    hic_file="hic.bam",
    nchrs=24,
    correct_nrounds=2,
    remove_allelic_links=2,
    generate_plots=True,
    verbose=True
)

success = processor.run_pipeline()
```

## 📋 输入要求 | Input Requirements

### Hi-C数据格式 | Hi-C Data Formats

1. **BAM文件**（推荐）
   - 按read name排序
   - 包含proper pair信息
   - MAPQ ≥ 1（默认）

2. **FASTQ文件**（自动执行BWA比对）
   - Paired-end FASTQ格式
   - 支持gzip压缩（.fastq.gz）
   - 自动执行：BWA mem + samblaster + HapHiC过滤

### 基因组组装文件 | Assembly File

- FASTA格式
- Contig或scaffold序列
- 无大小限制
- 支持分相组装（hifiasm输出）

## 📂 输出文件 | Output Files

### 主要输出 | Primary Outputs

```
04.build/
├── {prefix}.fa              # 最终scaffold序列
├── {prefix}.agp             # SALSA格式AGP文件
├── {prefix}.raw.agp         # YaHS格式AGP文件
└── juicebox.sh              # Juicebox脚本

05.plots/                    # 可视化图表目录
└── *.pdf/*.png              # Hi-C接触图

06.juicebox/                 # Juicebox兼容文件
├── {prefix}.hic            # Juicebox格式Hi-C文件
├── {prefix}.assembly       # 3D-DNA assembly文件
├── out.links.mnd           # MND格式文件
└── out.sorted.links.mnd    # 排序后的MND文件
```

### 日志文件 | Log Files

- `{prefix}_haphic.log`: 完整的运行日志
- 各步骤子目录中的详细日志

## ⚙️ 参数详解 | Parameter Details

### 核心参数 | Core Parameters

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `-a, --assembly` | - | 基因组组装文件（FASTA）**必需** |
| `-b, --bam` | - | Hi-C BAM文件 |
| `-1, --hic1` | - | Hi-C Read1文件（FASTQ） |
| `-2, --hic2` | - | Hi-C Read2文件（FASTQ） |
| `-c, --chr-number` | 12 | 染色体数量 |
| `-o, --output-dir` | - | 输出目录 |
| `--prefix` | - | 输出文件前缀 |

### 断点续传参数 | Resume Parameters

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--force-rerun` | False | 强制重新运行所有步骤，禁用断点续传 |

### 性能参数 | Performance Parameters

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--threads` | 8 | 线程数 |
| `--processes` | 8 | 并行进程数 |
| `--memory-limit` | - | 内存限制（如64G） |

### 组装校正参数 | Assembly Correction

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--correct-nrounds` | 2 | 组装校正轮数 |
| `--correct-min-coverage` | 10.0 | 校正最小覆盖度 |

### 聚类参数 | Clustering Parameters

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--min-inflation` | 1.0 | 最小膨胀参数 |
| `--max-inflation` | 3.0 | 最大膨胀参数 |
| `--inflation-step` | 0.2 | 膨胀参数步长 |
| `--Nx` | 80 | Nx参数 |

### 可视化参数 | Visualization Parameters

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--bin-size` | 500 | 接触图装箱大小 |
| `--min-len` | 1.0 | 最小scaffold长度 |
| `--no-generate-plots` | - | 禁用可视化图表生成 |

### Juicebox参数 | Juicebox Parameters

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--no-generate-juicebox` | - | 禁用Juicebox文件生成 |
| `--matlock-bin` | matlock | Matlock可执行文件路径 |
| `--three-d-dna-dir` | /share/... | 3D-DNA目录路径 |

## 🔬 Pipeline流程详解 | Pipeline Workflow

### 步骤1: BWA比对（如果输入为FASTQ）
```bash
# 自动执行
bwa mem -5SP assembly.fa read1.fastq.gz read2.fastq.gz |
samblaster |
samtools view -F 3340 - > HiC.bam

# HapHiC过滤
filter_bam -i HiC.bam -o HiC.filtered.bam -q 1 -e 3
```

### 步骤2: HapHiC Pipeline（一步完成）
```bash
haphic pipeline assembly.fa HiC.filtered.bam 12 \
    --aln_format bam \
    --outdir output_dir \
    --steps 1,2,3,4 \
    --threads 64 \
    --processes 64
```

子步骤包括：
- **Cluster**: 预处理和Markov聚类
- **Reassign**: 重新分配和拯救contigs
- **Sort**: 组内contig排序和定向
- **Build**: 构建最终scaffolds

### 步骤3: 可视化生成（默认执行）
```bash
haphic plot assembly.fa 04.build/prefix.agp \
    --bin_size 500
```

### 步骤4: Juicebox文件生成（可选）
```bash
# 生成MND文件
matlock bam2 juicer HiC.filtered.bam out.links.mnd

# 生成assembly文件
agp2assembly.py 04.build/prefix.agp prefix.assembly

# 生成.hic文件
run-assembly-visualizer.sh prefix.assembly out.sorted.links.mnd
```

## 📊 质量控制 | Quality Control

### 输入文件验证

- 自动检查文件存在性
- 验证BAM文件格式
- 检查FASTQ文件配对
- 统计基因组序列数量

### 系统资源检查

- CPU核心数检测
- 可用内存评估
- 磁盘空间检查

### 输出文件验证

- 关键文件存在性检查
- 文件大小合理性验证
- 统计信息生成

## 🔧 高级配置 | Advanced Configuration

### 单倍型分相 | Haplotype Phasing

```bash
# 移除等位基因连接
biopytools haphic -a assembly.fa -b hic.bam -c 24 \
    --remove-allelic-links 2

# 使用GFA文件
biopytools haphic -a assembly.fa -b hic.bam -c 24 \
    --gfa-files phased.gfa
```

### 调整聚类参数

```bash
# 密集聚类（适用于小基因组）
biopytools haphic -a assembly.fa -b hic.bam -c 24 \
    --min-inflation 1.0 --max-inflation 2.0

# 稀疏聚类（适用于大基因组）
biopytools haphic -a assembly.fa -b hic.bam -c 24 \
    --min-inflation 1.5 --max-inflation 4.0
```

### 内存优化

```bash
# 低内存配置
biopytools haphic -a assembly.fa -b hic.bam -c 24 \
    --threads 8 --processes 4 --memory-limit 32G

# 高性能配置
biopytools haphic -a assembly.fa -b hic.bam -c 24 \
    --threads 64 --processes 32
```

## 🚨 故障排除 | Troubleshooting

### 常见问题 | Common Issues

1. **内存不足**
   - 减少threads/processes数量
   - 使用--memory-limit参数

2. **BWA比对失败**
   - 检查FASTQ文件完整性
   - 确认文件权限

3. **HapHiC Pipeline失败**
   - 检查输入文件格式
   - 调整聚类参数
   - 使用--force-rerun重新运行

4. **可视化生成失败**
   - 确认AGP文件存在
   - 检查bin_size参数

### 日志分析 | Log Analysis

查看详细日志：
```bash
tail -f {prefix}_haphic.log
```

关键日志位置：
- `01.cluster/HapHiC_cluster.log`
- `02.reassign/HapHiC_reassign.log`
- `03.sort/HapHiC_sort.log`
- `04.build/HapHiC_build.log`

## 📈 性能基准 | Performance Benchmarks

### 典型运行时间 | Typical Runtime

| 基因组大小 | Contig数量 | 运行时间 | 内存使用 |
|------------|------------|----------|----------|
| 100 Mb | 5000 | 10分钟 | 4 GB |
| 1 Gb | 50000 | 30分钟 | 16 GB |
| 3 Gb | 200000 | 1小时 | 64 GB |

### 硬件建议 | Hardware Recommendations

- **最小配置**: 8核CPU, 16GB内存, 100GB磁盘
- **推荐配置**: 32核CPU, 64GB内存, 500GB磁盘
- **高性能配置**: 64核CPU, 128GB内存, 1TB磁盘

## 📚 参考文献 | References

1. Zeng X, Liu J, Li S, et al. HapHiC: a fast, reference-independent, allele-aware scaffolding tool based on Hi-C data. Bioinformatics. 2022.
2. Dudchenko O, Batra SS, Omer AD, et al. De novo assembly of the Aedes aegypti genome using Hi-C yields chromosome-level scaffolds. Science. 2017.
3. Burton JN, Adey A, Patwardhan RP, et al. Chromosome-scale scaffolding of de novo genome assemblies based on chromatin interactions. Nat Biotechnol. 2013.

## 🔗 相关链接 | Related Links

- [HapHiC GitHub](https://github.com/zengxiaofei/HapHiC)
- [HapHiC Documentation](https://github.com/zengxiaofei/HapHiC/wiki)
- [Juicebox](https://aidenlab.org/juicebox/)
- [3D-DNA](https://github.com/theaidenlab/3d-dna)

---

**版本信息 | Version**: v0.13.0+
**更新时间 | Last Updated**: 2024-12-20