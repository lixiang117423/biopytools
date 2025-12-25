# 🧬 Fastq到VCF (Parabricks) 全流程分析模块

**基于GPU加速的FASTQ到VCF全流程自动化分析工具 | GPU-Accelerated End-to-End FASTQ to VCF Automation Tool**

## 📖 功能概述 | Overview

Fastq到VCF (Parabricks) 模块是一个完整的端到端生物信息学分析流程，专为高通量测序数据分析而设计。该模块集成了质量控制、基因组索引构建、GPU加速序列比对、联合变异检测和变异过滤等关键步骤，基于NVIDIA Parabricks技术提供卓越的处理性能。

## ✨ 主要特性 | Key Features

- **🚀 GPU加速处理**: 基于NVIDIA Parabricks，比传统CPU工具快10-50倍
- **🔄 完整自动化流程**: FASTQ→质控→索引构建→比对→变异检测→变异过滤六步骤一体化
- **🧬 智能流程选择**: 根据样本数量自动选择GATK、GTX单机或GTX集群模式
- **⚙️ 高度可配置**: 支持自定义线程数、质量阈值、文件模式等参数
- **🛡️ 断点续传**: 支持检查点机制，可从中断步骤继续执行
- **📊 详细日志**: 完整的处理过程记录和性能监控
- **🎯 灵活执行**: 支持完整流程或单步骤执行
- **🔍 质量控制**: 集成biopytools fastp进行数据质控
- **💾 存储优化**: 智能的文件管理和空间规划

## 🚀 快速开始 | Quick Start

### 基本用法 | Basic Usage

```bash
# 完整的端到端分析（推荐）
biopytools fastq2vcf-parabricks \
    -i /path/to/raw_fastq \
    -r /path/to/reference.fa \
    -p /path/to/project

# 跳过质量控制步骤
biopytools fastq2vcf-parabricks \
    -i /path/to/raw_fastq \
    -r /path/to/reference.fa \
    -p /path/to/project \
    --skip-qc
```

### 高级用法 | Advanced Usage

```bash
# 自定义线程和质量参数
biopytools fastq2vcf-parabricks \
    -i /path/to/raw_fastq \
    -r /path/to/reference.fa \
    -p /path/to/project \
    --threads-mapping 64 \
    --threads-filter 32 \
    --snp-min-dp 10 \
    --snp-min-qual 30

# 分步骤执行
biopytools fastq2vcf-parabricks \
    -i /path/to/raw_fastq \
    -r /path/to/reference.fa \
    -p /path/to/project \
    --step 3  # 只执行第3步：序列比对
```

## 📋 分析流程 | Analysis Pipeline

### 🔄 完整流程步骤 | Complete Pipeline Steps

1. **🧹 Step 1: 质量控制 | Quality Control**
   - 使用biopytools fastp进行数据质控
   - 生成清洁的FASTQ文件和质控报告
   - 自动检测和修复测序数据质量问题

2. **📊 Step 2: 构建基因组索引 | Build Genome Index**
   - BWA索引构建（用于序列比对）
   - SAMtools faidx索引（用于快速访问）
   - GATK Sequence Dictionary（用于变异检测）
   - Parabricks专用索引（用于GPU加速）

3. **🗺️ Step 3: 序列比对 | Sequence Mapping**
   - GPU加速序列比对（Parabricks fq2bam）
   - 生成排序和去重的BAM文件
   - 生成gVCF文件用于联合分析
   - 支持多GPU并行处理

4. **🧬 Step 4: 联合变异检测 | Joint Variant Calling**
   - **GATK模式**: 样本数 < 4，使用GATK GenotypeGVCFs
   - **GTX单机模式**: 4 ≤ 样本数 < 200，使用GTX单机处理
   - **GTX集群模式**: 样本数 ≥ 200，生成集群作业脚本
   - 支持大规模队列分析

5. **🧹 Step 5: 变异过滤 | Variant Filtering**
   - SNP和InDel分别过滤
   - 基于深度和质量阈值过滤
   - 生成最终的过滤后VCF文件

## ⚙️ 参数配置 | Parameters

### 🔧 必需参数 | Required Parameters

```bash
-i, --input              # 原始FASTQ文件目录路径
-r, --ref-genome         # 参考基因组文件路径
-p, --project-base       # 项目根目录路径
```

### 📂 目录配置 | Directory Configuration

```bash
--clean-fastq-dir        # 清洁FASTQ文件目录路径
--mapping-dir           # 比对结果目录路径
--gvcf-dir              # gVCF文件目录路径
--bam-dir               # BAM文件目录路径
--joint-dir             # 联合检测结果目录路径
--filter-dir            # 过滤结果目录路径
```

### 🧵 性能配置 | Performance Configuration

```bash
--threads-mapping       # 比对线程数 (默认: 88)
--threads-filter        # 过滤线程数 (默认: 88)
--gatk-threshold        # GATK模式样本数阈值 (默认: 4)
--gtx-single-threshold  # GTX单机模式样本数阈值 (默认: 200)
--gtx-window-size       # GTX分块窗口大小 bp (默认: 20000000)
```

### 🎯 过滤参数 | Filtering Parameters

```bash
--snp-min-dp            # SNP最小深度 (默认: 5)
--snp-min-qual          # SNP最小质量 (默认: 30)
--indel-min-dp          # InDel最小深度 (默认: 5)
--indel-min-qual        # InDel最小质量 (默认: 30)
```

### 🔧 高级选项 | Advanced Options

```bash
--step [1-5]            # 只运行指定步骤
--dry-run              # 测试模式，不执行实际命令
--verbose              # 详细输出模式
--no-checkpoint        # 禁用断点续传
--skip-qc              # 跳过质控步骤
--skip-mapping         # 跳过比对步骤
```

## 📊 输出结构 | Output Structure

```
project_base/
├── 01.data/
│   ├── raw/           # 原始FASTQ文件
│   ├── clean/         # 质控后FASTQ文件
│   └── genome/        # 参考基因组和索引文件
├── 02.mapping/        # 比对结果文件
│   ├── bam/          # BAM文件
│   └── tmp/          # 临时文件
├── 03.gvcf/          # gVCF文件
├── 04.joint/         # 联合检测结果
├── 05.filter/        # 过滤后结果
├── 99.logs/          # 日志文件
├── .checkpoints/     # 检查点文件
└── ANALYSIS_REPORT.txt  # 分析报告
```

## 🎯 使用场景 | Use Cases

### 场景1: 小规模研究项目 (< 10样本)
```bash
# 使用GATK模式进行高精度分析
biopytools fastq2vcf-parabricks \
    -i ./raw_data \
    -r ./reference/genome.fa \
    -p ./project \
    --threads-mapping 32
```

### 场景2: 中等规模队列分析 (10-100样本)
```bash
# 使用GTX单机模式平衡速度和精度
biopytools fastq2vcf-parabricks \
    -i ./raw_data \
    -r ./reference/genome.fa \
    -p ./project \
    --threads-mapping 88 \
    --gtx-single-threshold 50
```

### 场景3: 大规模研究项目 (> 100样本)
```bash
# 使用GTX集群模式处理超大规模数据
biopytools fastq2vcf-parabricks \
    -i ./raw_data \
    -r ./reference/genome.fa \
    -p ./project \
    --gtx-single-threshold 100  # 提前切换到集群模式
```

### 场景4: 分步骤执行和调试
```bash
# 第1步: 质量控制
biopytools fastq2vcf-parabricks -i ./raw -r ./ref.fa -p ./proj --step 1

# 第2步: 构建索引
biopytools fastq2vcf-parabricks -i ./raw -r ./ref.fa -p ./proj --step 2

# 第3步: 序列比对
biopytools fastq2vcf-parabricks -i ./raw -r ./ref.fa -p ./proj --step 3

# 第4步: 联合变异检测
biopytools fastq2vcf-parabricks -i ./raw -r ./ref.fa -p ./proj --step 4

# 第5步: 变异过滤
biopytools fastq2vcf-parabricks -i ./raw -r ./ref.fa -p ./proj --step 5
```

## ⚡ 性能优化 | Performance Optimization

### 🎮 GPU资源配置
- **单GPU项目**: 88线程，适合小规模项目
- **多GPU项目**: 128线程，充分利用并行性能
- **集群项目**: 256线程，超大规模数据处理

### 💾 存储空间规划
- **输入数据**: 原始FASTQ文件大小
- **中间数据**: FASTQ大小×3-5倍（BAM + 临时文件）
- **输出数据**: FASTQ大小×0.5-1倍（VCF + 索引）
- **总空间建议**: 输入数据×5-7倍

### 🔄 批处理策略
1. **连续流程**: 适合<50样本的标准项目
2. **分步并行**: 适合50-200样本的中等规模项目
3. **集群处理**: 适合>200样本的大规模项目

## 🛠️ 故障排除 | Troubleshooting

### 常见问题及解决方案

**1️⃣ 质量控制失败**
```bash
# 检查原始数据格式和命名
ls -la /path/to/raw_fastq/*_1.fq.gz
ls -la /path/to/raw_fastq/*_2.fq.gz

# 重新执行质控步骤
biopytools fastq2vcf-parabricks -i ./raw -r ./ref.fa -p ./proj --step 1
```

**2️⃣ 索引构建失败**
```bash
# 检查参考基因组文件
samtools faidx reference.fa
head -5 reference.fa

# 清除检查点重新构建
rm -rf ./project/.checkpoints/
biopytools fastq2vcf-parabricks -i ./raw -r ./ref.fa -p ./proj --step 2
```

**3️⃣ 比对步骤失败**
```bash
# 检查质控结果
ls -la ./project/01.data/clean/

# 检查GPU可用性
nvidia-smi

# 重新执行比对
biopytools fastq2vcf-parabricks -i ./raw -r ./ref.fa -p ./proj --step 3
```

**4️⃣ 联合变异检测失败**
```bash
# 检查gVCF文件
ls -la ./project/03.gvcf/*.g.vcf.gz

# 检查样本数量和阈值
echo "样本数: $(ls ./project/03.gvcf/*.g.vcf.gz | wc -l)"

# 重新执行联合检测
biopytools fastq2vcf-parabricks -i ./raw -r ./ref.fa -p ./proj --step 4
```

**5️⃣ 内存不足问题**
```bash
# 调整GTX窗口大小
biopytools fastq2vcf-parabricks \
    -i ./raw -r ./ref.fa -p ./proj \
    --gtx-window-size 10000000  # 减小窗口大小
```

## 📚 相关资源 | Resources

- **NVIDIA Parabricks**: https://developer.nvidia.com/clara-parabricks
- **GATK Best Practices**: https://gatk.broadinstitute.org/hc/en-us
- **GTX Toolkit**: https://github.com/BGI-Genomics/GTX
- **BWA-MEM**: http://bio-bwa.sourceforge.net/bwa.shtml
- **SAMtools**: http://www.htslib.org/doc/samtools.html

## ⚠️ 注意事项 | Important Notes

- **GPU要求**: 需要NVIDIA GPU和CUDA支持
- **内存建议**: 大规模项目建议256GB+内存
- **存储规划**: 预留充足的磁盘空间（输入数据×5-7倍）
- **参考基因组**: 确保参考基因组文件完整且格式正确
- **文件命名**: 原始FASTQ文件需符合*_1.fq.gz和*_2.fq.gz命名模式
- **权限要求**: 需要对项目目录有读写权限
- **网络连接**: 首次运行可能需要下载依赖

## 🏆 最佳实践 | Best Practices

1. **测试先行**: 使用小数据集测试完整流程
2. **参数调优**: 根据数据特点调整质量阈值
3. **监控日志**: 关注处理进度和错误信息
4. **备份重要**: 定期备份中间结果和配置文件
5. **资源监控**: 监控GPU、内存和磁盘使用情况
6. **文档记录**: 记录分析参数和结果解读

---

🧬 **让复杂的基因组分析变得简单高效 | Making complex genomic analysis simple and efficient**