# HiCanu 基因组组装模块

**专业的基因组组装工具 | Professional Genome Assembly Tool**

## 功能概述 | Overview

HiCanu 基因组组装模块是一个专业的基因组组装工具，基于 HiCanu 软件构建，专门用于 PacBio HiFi reads 的基因组组装。HiCanu 是一个广泛应用于三代测序数据的组装软件，特别适合处理 PacBio HiFi 和 Nanopore 长读长数据。

## 主要特性 | Key Features

- **自动化组装流程**: 一键式 HiCanu 组装，从 reads 到 contigs
- **HiFi 优化**: 针对 PacBio HiFi 数据优化的参数配置
- **灵活参数配置**: 支持自定义基因组大小、最小 reads 长度等关键参数
- **多种组装阶段**: 支持 correct、trim、assemble 等不同阶段的单独运行
- **详细统计报告**: 自动解析组装结果，提供 N50、L50 等统计指标
- **质量控制**: 内置 reads 和组装结果的质量检查
- **完整日志记录**: 详细记录组装过程，便于问题排查
- **资源可控**: 支持线程数、内存等计算资源的配置

## 快速开始 | Quick Start

### 基本用法 | Basic Usage

```bash
# 基本组装
biopytools hicanu \
    -i reads.fastq \
    -g 120m \
    -p sample1 \
    -o hicanu_output

# 指定更多线程
biopytools hicanu \
    -i reads.fastq \
    -g 120m \
    -p sample1 \
    -o hicanu_output \
    -t 24
```

### 高级用法 | Advanced Usage

```bash
# 自定义组装参数
biopytools hicanu \
    -i reads.fastq \
    -g 120m \
    -p sample1 \
    -o hicanu_output \
    --min-read-length 2000 \
    --min-overlap-length 1000 \
    --corrected-error-rate 0.045 \
    -t 24

# 只运行纠错阶段
biopytools hicanu \
    -i reads.fastq \
    -g 120m \
    -p sample1 \
    -o correction_output \
    --stage correct

# 模拟运行（查看将要执行的命令）
biopytools hicanu \
    -i reads.fastq \
    -g 120m \
    -p sample1 \
    -o hicanu_output \
    --dry-run
```

## 参数说明 | Parameters

### 必需参数 | Required Parameters

| 参数 | 描述 | 示例 |
|------|------|------|
| `-i, --reads` | 输入 reads 文件路径（FASTA/FASTQ 格式）| `-i reads.fastq` |
| `-g, --genome-size` | 基因组大小（如：120m, 1g）| `-g 120m` |
| `-p, --prefix` | 输出文件前缀 | `-p sample1` |

### 软件配置 | Software Configuration

| 参数 | 描述 | 默认值 | 示例 |
|------|------|--------|------|
| `--canu-path` | HiCanu 可执行文件路径 | `/share/org/YZWL/yzwl_lixg/miniforge3/envs/canu_v.2.3/bin/canu` | `--canu-path /path/to/canu` |

### 输出配置 | Output Configuration

| 参数 | 描述 | 默认值 | 示例 |
|------|------|--------|------|
| `-o, --output-dir` | 输出目录路径 | `./hihicanu_output` | `-o results` |

### 组装参数 | Assembly Parameters

| 参数 | 描述 | 默认值 | 示例 |
|------|------|--------|------|
| `--min-read-length` | 最小 reads 长度 | `1000` | `--min-read-length 2000` |
| `--min-overlap-length` | 最小重叠长度 | `500` | `--min-overlap-length 1000` |
| `--corrected-error-rate` | 纠正错误率 | `自动` | `--corrected-error-rate 0.045` |
| `--raw-error-rate` | 原始错误率 | `自动` | `--raw-error-rate 0.300` |
| `--max-input-coverage` | 最大输入覆盖度 | `自动` | `--max-input-coverage 200` |
| `--stage` | 组装阶段 | `assemble` | `--stage trim-assemble` |

**组装阶段选项 | Assembly Stage Options:**
- `haplotype`: 生成单倍型特异性 reads
- `correct`: 生成纠正后的 reads
- `trim`: 生成修剪后的 reads
- `assemble`: 生成组装（默认）
- `trim-assemble`: 修剪后组装

### 计算资源 | Computing Resources

| 参数 | 描述 | 默认值 | 示例 |
|------|------|--------|------|
| `-t, --threads` | 线程数 | `12` | `-t 24` |
| `-m, --memory` | 内存限制 | `80G` | `-m 120G` |
| `--use-grid` | 使用网格调度 | `False` | `--use-grid` |
| `--grid-options` | 网格调度选项 | `None` | `--grid-options "-P project_q"` |

### 执行控制 | Execution Control

| 参数 | 描述 | 默认值 | 示例 |
|------|------|--------|------|
| `--dry-run` | 模拟运行（不实际执行）| `False` | `--dry-run` |
| `--keep-intermediate` | 保留中间文件 | `False` | `--keep-intermediate` |

## 输出文件 | Output Files

组装完成后，将在输出目录中生成以下文件：

```
hihicanu_output/
├── sample1.contigs.fasta          # 组装的 contigs 序列
├── sample1.unitigs.fasta          # 组装的 unitigs 序列
├── sample1.report                 # HiCanu 组装报告
└── hicanu_YYYYMMDD_HHMMSS.log  # 运行日志
```

### 主要输出文件说明 | Main Output Files Description

| 文件 | 描述 |
|------|------|
| `*.contigs.fasta` | 组装的 contigs 序列，是最终的组装结果 |
| `*.unitigs.fasta` | 组装的 unitigs 序列，是 contigs 的前体 |
| `*.report` | HiCanu 组装过程的详细报告 |
| `*.log` | 本模块的运行日志 |

## 组装统计 | Assembly Statistics

程序会自动分析组装结果并输出以下统计信息：

- **序列数量 | Number of sequences**: contigs/unitigs 的总数
- **总长度 | Total bp**: 组装的总碱基数
- **N50**: 组装的 N50 值（按长度排序的 contigs 累计达到总长度 50% 时的 contig 长度）
- **L50**: 达到 N50 所需的 contigs 数量
- **最小长度 | Min length**: 最短 contig 的长度
- **最大长度 | Max length**: 最长 contig 的长度
- **平均长度 | Average length**: contigs 的平均长度

## 使用示例 | Usage Examples

### 示例 1: 基本组装 | Example 1: Basic Assembly

```bash
biopytools hicanu \
    -i pacbio_hifi_reads.fastq \
    -g 120m \
    -p strain_A \
    -o assembly_results
```

**预期结果 | Expected Results:**
- 生成 `strain_A.contigs.fasta`
- 输出组装统计信息
- 记录详细日志

### 示例 2: 高覆盖度数据 | Example 2: High Coverage Data

对于高覆盖度数据，可以限制输入覆盖度以加快处理：

```bash
biopytools hicanu \
    -i high_coverage_reads.fastq \
    -g 120m \
    -p sample1 \
    -o results \
    --max-input-coverage 150 \
    -t 32
```

### 示例 3: 分阶段运行 | Example 3: Run by Stages

```bash
# 第一步：纠正 reads
biopytools hicanu \
    -i reads.fastq \
    -g 120m \
    -p sample1 \
    -o corrected \
    --stage correct

# 第二步：修剪并组装
biopytools hicanu \
    -i corrected/sample1.correctedReads.fasta.gz \
    -g 120m \
    -p sample1 \
    -o final_assembly \
    --stage trim-assemble
```

### 示例 4: 质量控制组装 | Example 4: Quality Controlled Assembly

对于质量较低的数据，可以调整错误率参数：

```bash
biopytools hicanu \
    -i lower_quality_reads.fastq \
    -g 120m \
    -p sample1 \
    -o results \
    --raw-error-rate 0.400 \
    --corrected-error-rate 0.060 \
    --min-read-length 2000
```

## Python API 使用 | Python API Usage

除了命令行接口，也可以直接在 Python 代码中使用：

```python
from biopytools.hicanu import HiCanuAssemblyPipeline

# 创建组装流程
pipeline = HiCanuAssemblyPipeline(
    reads_file="reads.fastq",
    genome_size="120m",
    prefix="sample1",
    output_dir="./hicanu_output",
    threads=24
)

# 运行组装
pipeline.run()
```

## 常见问题 | FAQ

### Q1: 如何确定 genomeSize 参数？

**A:** genomeSize 应该是你对基因组大小的最佳估计值，用于估计覆盖度，而不是限制组装大小。可以通过以下方式确定：
- 流式细胞仪检测
- Kmer 分析（如 GenomeScope）
- 近缘物种的基因组大小

### Q2: 组装需要多长时间？

**A:** 组装时间取决于：
- 数据量（reads 数量和总碱基数）
- 基因组复杂度
- 计算资源（线程数、内存）

典型情况下，一个 100Mb 基因组的组装可能需要数小时到数天。

### Q3: 如何评估组装质量？

**A:** 可以通过以下指标评估：
- **总长度**: 是否接近预期基因组大小
- **N50**: 越大表示 contigs 越长，质量越好
- **Contigs 数量**: 越少通常表示组装越完整
- **BUSCO**: 使用 BUSCO 评估基因组完整性

### Q4: 内存不足怎么办？

**A:** HiCanu 会根据 genomeSize 自动调整内存使用。如果仍然不足：
1. 减少 `--max-input-coverage`
2. 分阶段运行（`--stage correct`，然后 `--stage trim-assemble`）
3. 增加可用内存或使用网格调度

### Q5: 为什么 contigs.fasta 和 unitigs.fasta 都存在？

**A:**
- `unitigs.fasta`: 未进行 bubble popping 的初级组装结果
- `contigs.fasta`: 经过 bubble popping 优化的最终组装结果

通常使用 `contigs.fasta` 作为最终组装结果。

## 技术支持 | Technical Support

如有问题或建议，请联系：
- 作者 | Author: Xiang LI
- 版本 | Version: 1.0.0
- 日期 | Date: 2025-12-30

## 参考资源 | References

- [HiCanu 官方文档](http://canu.readthedocs.org/en/latest/)
- [PacBio HiFi 测序技术](https://www.pacb.com/technologies/hifi-sequencing/)
- [基因组组装最佳实践](https://github.com/marbl/canu)
