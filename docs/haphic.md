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
04_build/
├── {prefix}.fa              # 最终scaffold序列
├── {prefix}.agp             # SALSA格式AGP文件
├── {prefix}.raw.agp         # YaHS格式AGP文件
└── juicebox.sh              # Juicebox脚本

05_plots/                    # 可视化图表目录
└── *.pdf/*.png              # Hi-C接触图

06_juicebox/                 # Juicebox兼容文件
├── {prefix}.hic            # Juicebox格式Hi-C文件
├── {prefix}.assembly       # 3D-DNA assembly文件
├── out.links.mnd           # MND格式文件
└── out.sorted.links.mnd    # 排序后的MND文件
```

### 日志文件 | Log Files

- `99_logs/{prefix}_haphic.log`: 完整的运行日志
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
haphic plot assembly.fa 04_build/prefix.agp \
    --bin_size 500
```

### 步骤4: Juicebox文件生成（可选）
```bash
# 生成MND文件
matlock bam2 juicer HiC.filtered.bam out.links.mnd

# 生成assembly文件
agp2assembly.py 04_build/prefix.agp prefix.assembly

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
tail -f 99_logs/{prefix}_haphic.log
```

关键日志位置：
- `01_cluster/HapHiC_cluster.log`
- `02_reassign/HapHiC_reassign.log`
- `03_sort/HapHiC_sort.log`
- `04_build/HapHiC_build.log`

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

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input, -i` | 必填 | Path | 基因组组装文件(FASTA)｜Genome assembly file path (FASTA format) |
| `--bam, -b` | — | Path | Hi-C BAM文件(按read名排序)｜Hi-C BAM file path (sorted by read name) |
| `--hic1, -1` | — | Path | Hi-C Read1文件(FASTQ)｜Hi-C Read1 file path (FASTQ format) |
| `--hic2, -2` | — | Path | Hi-C Read2文件(FASTQ)｜Hi-C Read2 file path (FASTQ format) |
| `--chr-number, -c` | `12` | int | 染色体数量｜Number of chromosomes |
| `--output-dir, -o` | `.` | Path | 输出目录路径｜Output directory path |
| `--prefix` | — |  | 输出文件前缀｜Output file prefix |
| `--force-rerun` | — |  | 强制重新运行所有步骤｜Force rerun all steps (disable resume mode) |
| `--mapq-threshold` | `1` | int | MAPQ阈值｜MAPQ threshold |
| `--edit-distance` | `3` | int | 编辑距离阈值｜Edit distance threshold |
| `--re-site-cutoff` | `5` | int | Step1 RE位点过滤阈值｜Step1 RE site filtering threshold |
| `--min-RE-sites` | `25` | int | Step2重分配最小RE位点数｜Step2 reassignment min RE sites |
| `--aln-format` | `auto` | auto/bam/pairs | 比对文件格式｜Alignment file format |
| `--min-inflation` | `1.1` | float | 最小膨胀值｜Min inflation |
| `--max-inflation` | `3.0` | float | 最大膨胀值｜Max inflation |
| `--inflation-step` | `0.1` | float | 膨胀值步长｜Inflation step |
| `--Nx` | `80` | int | Nx参数｜Nx parameter |
| `--min-group-len` | `5.0` | float | 最小组长度(Mbp)｜Min group length (Mbp) |
| `--flank` | `500` | int | 邻接矩阵侧翼区域(kbp)｜Adjacency matrix flank region (kbp) |
| `--bin-size-kbp` | `-1` | int | 聚类分箱大小(kbp),-1=自动｜Clustering bin size (kbp), -1=auto |
| `--processes` | `8` | int | 并行进程数｜Number of parallel processes |
| `--threads, -t` | `12` | int | 线程数｜Number of threads |
| `--memory-limit` | `100G` |  | 内存限制｜Memory limit (e.g., 64G, 300G) |
| `--correct-nrounds` | `2` | int | 组装修正轮数(0=禁用)｜Assembly correction rounds (0=disabled) |
| `--correct-min-coverage` | `10.0` | float | 修正最小覆盖度｜Correction min coverage |
| `--median-cov-ratio` | `0.2` | float | 覆盖率截断乘数｜Coverage cutoff multiplier |
| `--region-len-ratio` | `0.1` | float | 高覆盖区域长度比｜High-coverage region length ratio |
| `--min-region-cutoff` | `5000` | int | 高覆盖区域最小长度(bp)｜Min high-coverage region length (bp) |
| `--skip-fast-sort` | — |  | 跳过快速排序｜Skip fast sorting |
| `--skip-allhic` | — |  | 跳过ALLHiC优化｜Skip ALLHiC optimization |
| `--skip-ga` | — |  | 跳过ALLHiC遗传算法｜Skip ALLHiC genetic algorithm |
| `--sort-by-input` | — |  | 按输入顺序排序｜Sort output by input order |
| `--no-additional-rescue` | — |  | 跳过额外救援轮｜Skip additional rescue round |
| `--remove-concentrated-links` | — |  | 移除高密度集中链接｜Remove concentrated links |
| `--normalize-by-nlinks` | — |  | 按链接数归一化｜Normalize by number of links |
| `--dense-matrix` | — |  | 使用稠密矩阵｜Use dense matrix |
| `--remove-allelic-links` | — | int | 移除等位基因连锁数｜Remove allelic links count |
| `--phasing-weight` | `1.0` | float | 相位权重｜Phasing weight |
| `--gfa-files` | — |  | GFA文件路径(逗号分隔)｜GFA files path (comma-separated) |
| `--generate-plots` | — |  | 生成可视化图表｜Generate visualization plots |
| `--bin-size` | `500` | int | 接触图bin大小｜Contact map bin size |
| `--min-len` | `1.0` | float | 最小scaffold长度｜Min scaffold length |
| `--separate-plots` | — |  | 生成单独的图表｜Generate separate plots |
| `--haphic-bin` | `~/miniforge3/envs/hic/bin/haphic` |  | HapHiC可执行文件路径｜HapHiC executable path |
| `--bwa-bin` | `~/miniforge3/envs/align/bin/bwa` |  | BWA可执行文件路径｜BWA executable path |
| `--samtools-bin` | `~/miniforge3/envs/align/bin/samtools` |  | Samtools可执行文件路径｜Samtools executable path |
| `--samblaster-bin` | `~/miniforge3/envs/hic/bin/samblaster` |  | Samblaster可执行文件路径｜Samblaster executable path |
| `--haphic-filter-bam-bin` | `~/miniforge3/envs/hic/bin/filter_bam` |  | HapHiC filter_bam工具路径｜HapHiC filter_bam tool path |
| `--use-samblaster/--no-use-samblaster` | `True` |  | 使用samblaster去重｜Use samblaster deduplication |
| `--use-haphic-filter/--no-use-haphic-filter` | `True` |  | 使用HapHiC过滤｜Use HapHiC filtering |
| `--generate-juicebox/--no-generate-juicebox` | `True` |  | 生成Juicebox兼容文件｜Generate Juicebox compatible files |
| `--matlock-bin` | `~/miniforge3/envs/hic/bin/matlock` |  | Matlock可执行文件路径｜Matlock executable path |
| `--three-d-dna-dir` | `~/software/3d-dna` |  | 3D-DNA目录路径｜3D-DNA directory path |
| `--agp2assembly-script` | `~/software/3d-dna/utils/agp2assembly.py` |  | agp2assembly脚本路径｜agp2assembly script path |
| `--asm-visualizer-script` | `~/software/3d-dna/visualize/run-assembly-visualizer.sh` |  | asm-visualizer脚本路径｜asm-visualizer script path |
| `--RE` | `GATC` |  | 限制性酶切位点｜Restriction enzyme sites |
| `--quick-view` | — |  | 快速查看模式｜Quick view mode |
| `--no-agp` | — |  | 不输出AGP文件｜Don't output AGP file |
| `--no-fasta` | — |  | 不输出FASTA文件｜Don't output FASTA file |
| `--no-generate-plots` | — |  | 不生成可视化图表｜Don't generate visualization plots |
| `--no-juicebox` | — |  | 不生成Juicebox脚本｜Don't generate Juicebox script |
| `--verbose, -v` | — |  | 详细输出模式｜Verbose output mode |
| `--log-file` | — | Path | 日志文件路径｜Log file path |
| `--dry-run` | — |  | 测试模式,不执行实际命令｜Test mode, do not execute actual commands |

<!-- END PARAMS:auto -->
