# 🚀 Parabricks WGS分析模块

**基于GPU加速的全基因组测序数据分析工具 | GPU-Accelerated Whole Genome Sequencing Data Analysis Tool**

## 📖 功能概述 | Overview

Parabricks WGS分析模块是一个基于NVIDIA Parabricks的GPU加速全基因组测序数据分析工具，提供从FASTQ原始数据到最终变异检测结果的完整流程。支持快速序列比对、GVCF生成、Joint Calling联合变异检测，以及多种质量控制和分析选项，专为大规模基因组数据分析而优化。

## ✨ 主要特性 | Key Features

- **🚀 GPU加速处理**: 基于NVIDIA Parabricks，相比传统CPU工具提升数十倍处理速度
- **🔄 完整分析流程**: FASTQ→BAM→GVCF→Joint Calling→变异检测五步骤自动化
- **🎯 灵活输出格式**: 支持BAM、VCF、GVCF多种输出格式
- **🧬 联合变异检测**: 支持多样本Joint Calling联合分析，提高检测准确性
- **⚙️ 高度可配置**: 自定义软件路径、线程数、质量阈值等参数
- **🛡️ 智能质量控制**: 可配置的质量过滤和置信度阈值
- **📊 详细日志记录**: 完整的处理过程日志和性能监控
- **🐳 容器化支持**: 支持Singularity容器环境部署
- **🔍 批处理优化**: 支持多样本批量处理和资源管理

## 🚀 快速开始 | Quick Start

### 基本用法 | Basic Usage

```bash
# GPU加速的WGS分析（默认输出GVCF）
biopytools parabricks \
    -i /path/to/clean_fastq \
    -o /path/to/output \
    -r /path/to/reference.fa

# 使用自定义线程数和质量阈值
biopytools parabricks \
    -i /path/to/clean_fastq \
    -o /path/to/output \
    -r /path/to/reference.fa \
    --threads 64 \
    --min-confidence 20
```

### 高级用法 | Advanced Usage

```bash
# 完整的Joint Calling分析流程
biopytools parabricks \
    -i /path/to/clean_fastq \
    -o /path/to/output \
    -r /path/to/reference.fa \
    --threads 88 \
    --gvcf \
    --joint-calling \
    --combined-output-name cohort.g.vcf \
    --min-base-quality 30 \
    --ploidy 2 \
    --tmp-dir /scratch/tmp

# 自定义容器路径和文件模式
biopytools parabricks \
    -i /path/to/fastq \
    -o /path/to/output \
    -r /path/to/reference.fa \
    --parabricks-path /custom/path/parabricks.sif \
    --read1-pattern "*_R1.clean.fq.gz" \
    --read2-pattern "*_R2.clean.fq.gz"
```

## 📋 命令行参数 | Command Line Parameters

### 必需参数 | Required Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-i, --input-dir` | `None` | 📁 输入目录路径（包含清洁的FASTQ文件） |
| `-o, --output-dir` | `None` | 📂 输出目录路径 |
| `-r, --reference` | `None` | 🧬 参考基因组文件路径 |

### 性能配置参数 | Performance Configuration

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--threads` | `88` | 🔧 并行线程数 |
| `--parabricks-path` | `/share/apps/containers/parabricks.sif` | 🛠️ Parabricks容器路径 |
| `--tmp-dir` | `None` | 🗂️ 临时目录路径 |

### 质量控制参数 | Quality Control Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--min-confidence-threshold` | `30` | 🎯 最小置信度阈值 |
| `--min-base-quality` | `20` | 🎯 最小碱基质量 |
| `--pcr-indel-model` | `CONSERVATIVE` | 🧬 PCR InDel模型 |
| `--ploidy` | `2` | 🔢 样本倍性 |

### 输出格式参数 | Output Format Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--gvcf` | `True` | 📄 输出GVCF格式文件 |
| `--joint-calling` | `True` | 🧬 执行Joint Calling联合变异检测 |
| `--combined-output-name` | `combined.g.vcf` | 📝 合并输出文件名 |

### 文件模式参数 | File Pattern Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--read1-pattern` | `*_1.clean.fq.gz` | 📄 R1文件匹配模式 |
| `--read2-pattern` | `*_2.clean.fq.gz` | 📄 R2文件匹配模式 |

### 步骤说明 | Step Descriptions

| 步骤 | 名称 | 描述 |
|------|------|------|
| **1** | 🧬 FASTQ→BAM | GPU加速序列比对和排序 |
| **2** | 📄 BAM→GVCF | 生成基因组变异格式文件 |
| **3** | 🧬 Joint Calling | 多样本联合变异检测 |
| **4** | 🧪 变异过滤 | 质量控制和变异过滤 |
| **5** | 📊 结果汇总 | 生成分析报告和统计信息 |

## 📁 输入文件格式 | Input File Formats

### 清洁的FASTQ文件 | Clean FASTQ Files

经过质量控制的paired-end FASTQ文件：

```bash
# 文件命名示例
sample1_1.clean.fq.gz
sample1_2.clean.fq.gz
sample2_1.clean.fq.gz
sample2_2.clean.fq.gz
```

**文件要求**:
- gzip压缩的FASTQ格式
- 配对的R1和R2文件
- 质量控制在Q30以上
- 片段长度150bp或以上

### 参考基因组文件 | Reference Genome File

标准FASTA格式的参考基因组：

```fasta
>chromosome1
ATCGATCGATCGATCGATCGATCGATCG...
>chromosome2
GCTAGCTAGCTAGCTAGCTAGCTAGCTA...
```

**文件要求**:
- 标准FASTA格式
- 包含所有染色体的完整序列
- 经过indexing（可通过samtools faidx）
- 建议使用hg38、hg19等标准参考基因组

## 📂 输出文件说明 | Output Files Description

### BAM文件 | BAM Files

排序后的比对文件：
```
output_dir/
├── bam/
│   ├── sample1.sorted.bam
│   ├── sample1.sorted.bam.bai
│   ├── sample2.sorted.bam
│   └── sample2.sorted.bam.bai
```

### GVCF文件 | GVCF Files

基因组变异格式文件：
```
output_dir/
├── gvcf/
│   ├── sample1.g.vcf.gz
│   ├── sample1.g.vcf.gz.tbi
│   ├── sample2.g.vcf.gz
│   └── sample2.g.vcf.gz.tbi
```

### Joint Calling结果 | Joint Calling Results

联合变异检测结果：
```
output_dir/
├── joint_calling/
│   ├── cohort.g.vcf.gz
│   ├── cohort.g.vcf.gz.tbi
│   └── joint_calling.log
```

### 分析报告 | Analysis Reports

详细的分析报告和统计信息：
```
output_dir/
├── reports/
│   ├── summary_report.html
│   ├── quality_metrics.txt
│   └── performance_stats.txt
```

## 🔧 使用示例 | Usage Examples

### 单样本快速分析

```bash
# 基本单样本分析
biopytools parabricks \
    -i /data/clean_fastq \
    -o /results/sample_analysis \
    -r /reference/hg38.fa \
    --threads 32

# 输出GVCF格式
biopytools parabricks \
    -i /data/clean_fastq \
    -o /results/gvcf_output \
    -r /reference/hg38.fa \
    --gvcf \
    --min-confidence 30
```

### 多样本Joint Calling分析

```bash
# 批量多样本分析
biopytools parabricks \
    -i /data/cohort_fastq \
    -o /results/cohort_analysis \
    -r /reference/hg38.fa \
    --gvcf \
    --joint-calling \
    --combined-output-name cohort_combined.g.vcf \
    --threads 88 \
    --tmp-dir /scratch/temp

# 高质量变异检测
biopytools parabricks \
    -i /data/high_quality_fastq \
    -o /results/stringent_analysis \
    -r /reference/hg38.fa \
    --min-confidence 50 \
    --min-base-quality 30 \
    --pcr-indel-model CONSERVATIVE
```

### 自定义容器和路径

```bash
# 使用自定义Parabricks容器
biopytools parabricks \
    -i /data/fastq \
    -o /results/custom_analysis \
    -r /reference/genome.fa \
    --parabricks-path /custom/containers/parabricks_v4.0.sif \
    --tmp-dir /high_speed_ssd/tmp

# 特定文件命名模式
biopytools parabricks \
    -i /data/project_fastq \
    -o /results/project_analysis \
    -r /reference/genome.fa \
    --read1-pattern "*_R1.clean.fq.gz" \
    --read2-pattern "*_R2.clean.fq.gz"
```

## ⚡ 性能优化 | Performance Optimization

### GPU资源配置建议

| 样本数量 | GPU内存 | 线程数 | 建议配置 |
|----------|---------|--------|----------|
| 1-10 | 16GB | 32 | 基本分析 |
| 10-50 | 32GB | 64 | 批量处理 |
| 50+ | 64GB | 88+ | 大规模项目 |

### 存储和内存优化

```bash
# 使用高速SSD作为临时目录
biopytools parabricks \
    -i /data/fastq \
    -o /results/output \
    -r /reference/genome.fa \
    --tmp-dir /nvme/tmp \
    --threads 88

# 并行处理多个样本
for sample in $(cat sample_list.txt); do
    biopytools parabricks \
        -i "/data/${sample}/clean_fastq" \
        -o "/results/${sample}" \
        -r /reference/hg38.fa \
        --threads 32 &
done
wait
```

## 🐛 故障排除 | Troubleshooting

### 常见问题 | Common Issues

**问题1: GPU内存不足**
```bash
# 解决方案：减少线程数或使用更小的批次
biopytools parabricks \
    -i /data/fastq \
    -o /results/output \
    -r /reference/genome.fa \
    --threads 16  # 减少线程数
```

**问题2: 临时目录空间不足**
```bash
# 解决方案：指定更大的临时目录
biopytools parabricks \
    -i /data/fastq \
    -o /results/output \
    -r /reference/genome.fa \
    --tmp-dir /large_disk/tmp
```

**问题3: 容器权限问题**
```bash
# 确保容器有读写权限
chmod +x /share/apps/containers/parabricks.sif
export SINGULARITY_BINDPATH="/data,/results"
```

### 性能监控

```bash
# 监控GPU使用情况
nvidia-smi -l 1

# 监控内存和CPU使用
htop

# 检查临时目录空间
df -h /tmp
```

## 📊 输出质量评估 | Output Quality Assessment

### BAM文件质量检查

```bash
# 使用samtools统计比对质量
samtools flagstat sample.sorted.bam
samtools stats sample.sorted.bam > bam_stats.txt

# 使用qualimap评估覆盖度
qualimap bamqc -bam sample.sorted.bam -outdir coverage_analysis
```

### VCF文件质量检查

```bash
# 使用bcftools统计变异信息
bcftools stats cohort.g.vcf.gz > vcf_stats.txt
plot-vcfstats vcf_stats.txt -p vcf_plots

# 检查变异质量分布
bcftools query -f '%CHROM\t%POS\t%QUAL\t%DP\n' cohort.g.vcf.gz | \
    awk '$3>30 && $4>10' > high_quality_variants.txt
```

## 🔗 相关文档 | Related Documentation

- [NVIDIA Parabricks官方文档](https://docs.nvidia.com/clara/parabricks/)
- [GATK Best Practices](https://gatk.broadinstitute.org/hc/en-us/articles/360035531132)
- [biopytools快速开始指南](../README.md)
- [Annovar变异注释模块](annovar.md)
- [Fastq2VCF Parabricks流程](fastq2vcf_parabricks.md)

## 📄 许可证 | License

本模块遵循MIT许可证。详细信息请参见LICENSE文件。

## 🤝 贡献指南 | Contributing

欢迎提交Issue和Pull Request来改进本模块。

## 📞 技术支持 | Support

如有技术问题，请联系：
- 邮箱: yzwl_lixg@outlook.com
- 项目地址: https://github.com/your-org/biopytools

---

**最后更新**: 2024年12月17日
**版本**: 1.0.0
**作者**: biopytools开发团队