# 🔬 Fastq到VCF (GTX) 全流程分析模块

**基于GTX的CPU优化FASTQ到VCF全流程自动化分析工具 | CPU-Optimized End-to-End FASTQ to VCF Automation Tool**

## 📖 功能概述 | Overview

Fastq到VCF (GTX) 模块是一个专为CPU环境优化的全流程基因组分析工具，集成了质量控制、基因组索引构建、高效序列比对、联合变异检测和变异过滤等完整流程。基于GTX（Genome Toolkit）技术，在CPU环境下提供卓越的处理性能和内存效率。

## ✨ 主要特性 | Key Features

- **🚀 CPU优化设计**: 专门针对CPU环境优化，无需GPU支持
- **🔄 完整自动化流程**: FASTQ→质控→索引构建→比对→变异检测→变异过滤六步骤一体化
- **🧬 智能流程选择**: 根据样本数量自动选择GATK、GTX单机或GTX集群模式
- **⚙️ 高度可配置**: 支持自定义线程数、质量阈值、窗口大小等参数
- **🛡️ 断点续传**: 支持检查点机制，可从中断步骤继续执行
- **📊 详细日志**: 完整的处理过程记录和性能监控
- **🎯 灵活执行**: 支持完整流程或单步骤执行
- **🔍 质量控制**: 集成biopytools fastp进行数据质控
- **💾 内存优化**: 智能的内存管理和分块处理策略

## 🚀 快速开始 | Quick Start

### 基本用法 | Basic Usage

```bash
# 完整的端到端分析（推荐）
biopytools fastq2vcf-gtx \
    -i /path/to/raw_fastq \
    -r /path/to/reference.fa \
    -p /path/to/project

# 跳过质量控制步骤
biopytools fastq2vcf-gtx \
    -i /path/to/raw_fastq \
    -r /path/to/reference.fa \
    -p /path/to/project \
    --skip-qc
```

### 高级用法 | Advanced Usage

```bash
# 自定义GTX二进制路径和线程参数
biopytools fastq2vcf-gtx \
    -i /path/to/raw_fastq \
    -r /path/to/reference.fa \
    -p /path/to/project \
    --gtx-bin /custom/path/to/gtx \
    --threads-gtx 64 \
    --gtx-window-size 10000000

# 分步骤执行
biopytools fastq2vcf-gtx \
    -i /path/to/raw_fastq \
    -r /path/to/reference.fa \
    -p /path/to/project \
    --step 4  # 只执行第4步：联合变异检测
```

## 📋 分析流程 | Analysis Pipeline

### 🔄 完整流程步骤 | Complete Pipeline Steps

1. **🧹 Step 1: 质量控制 | Quality Control**
   - 使用biopytools fastp进行数据质控
   - 生成清洁的FASTQ文件和质控报告
   - **立即执行GTX索引构建**: 在质控后直接构建GTX索引

2. **📊 Step 2: 构建其他索引 | Build Other Indexes**
   - BWA索引构建（用于序列比对）
   - SAMtools faidx索引（用于快速访问）
   - GATK Sequence Dictionary（用于变异检测）
   - GTX索引已在Step 1中构建完成

3. **🗺️ Step 3: 序列比对 | Sequence Mapping**
   - CPU优化序列比对（GTX WGS）
   - 一体化比对和变异检测流程
   - 生成排序的BAM文件和gVCF文件
   - 支持分块处理大数据集

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
--threads-gtx          # GTX处理线程数 (默认: 88)
--threads-filter       # 过滤线程数 (默认: 88)
--gatk-threshold       # GATK模式样本数阈值 (默认: 4)
--gtx-single-threshold # GTX单机模式样本数阈值 (默认: 200)
--gtx-window-size      # GTX分块窗口大小 bp (默认: 20000000)
```

### 🛠️ 工具路径 | Tool Paths

```bash
--gtx-bin              # GTX可执行文件路径
                       # (默认: /share/org/YZWL/yzwl_lixg/software/gtx/bin/gtx)
```

### 🎯 过滤参数 | Filtering Parameters

```bash
--snp-min-dp           # SNP最小深度 (默认: 5)
--snp-min-qual         # SNP最小质量 (默认: 30)
--indel-min-dp         # InDel最小深度 (默认: 5)
--indel-min-qual       # InDel最小质量 (默认: 30)
```

### 🔧 高级选项 | Advanced Options

```bash
--step [1-5]           # 只运行指定步骤
--dry-run             # 测试模式，不执行实际命令
--verbose             # 详细输出模式
--no-checkpoint       # 禁用断点续传
--skip-qc             # 跳过质控步骤
--skip-mapping        # 跳过比对步骤
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
├── .tmp/             # GTX临时目录
└── ANALYSIS_REPORT.txt  # 分析报告
```

## 🎯 使用场景 | Use Cases

### 场景1: CPU环境下的标准项目
```bash
# 标准CPU环境下的完整分析
biopytools fastq2vcf-gtx \
    -i ./raw_data \
    -r ./reference/genome.fa \
    -p ./project \
    --threads-gtx 64
```

### 场景2: 大规模CPU集群分析
```bash
# 大规模队列分析，使用集群模式
biopytools fastq2vcf-gtx \
    -i ./raw_data \
    -r ./reference/genome.fa \
    -p ./project \
    --threads-gtx 128 \
    --gtx-single-threshold 100  # 提前切换到集群模式
    --gtx-window-size 10000000  # 优化内存使用
```

### 场景3: 资源受限环境
```bash
# 低资源配置环境
biopytools fastq2vcf-gtx \
    -i ./raw_data \
    -r ./reference/genome.fa \
    -p ./project \
    --threads-gtx 32 \
    --gtx-window-size 5000000  # 更小的窗口大小
```

### 场景4: 分步骤执行和调试
```bash
# 第1步: 质量控制 + GTX索引构建
biopytools fastq2vcf-gtx -i ./raw -r ./ref.fa -p ./proj --step 1

# 第2步: 构建其他索引
biopytools fastq2vcf-gtx -i ./raw -r ./ref.fa -p ./proj --step 2

# 第3步: GTX WGS比对
biopytools fastq2vcf-gtx -i ./raw -r ./ref.fa -p ./proj --step 3

# 第4步: 联合变异检测
biopytools fastq2vcf-gtx -i ./raw -r ./ref.fa -p ./proj --step 4

# 第5步: 变异过滤
biopytools fastq2vcf-gtx -i ./raw -r ./ref.fa -p ./proj --step 5
```

## ⚡ 性能优化 | Performance Optimization

### 🖥️ CPU资源配置
- **小规模项目**: 32-64线程，适合<20样本
- **中等规模项目**: 64-128线程，适合20-100样本
- **大规模项目**: 128+线程，适合>100样本

### 💾 内存管理策略
- **窗口大小调整**: 根据可用内存调整`--gtx-window-size`
- **分块处理**: GTX自动将大数据集分块处理
- **临时文件**: 使用项目目录下的`.tmp/`目录

### 🔄 批处理优化
1. **连续流程**: 适合<50样本
2. **分块处理**: 使用窗口大小优化内存使用
3. **集群处理**: 超大规模数据分布式处理

## 🛠️ 故障排除 | Troubleshooting

### 常见问题及解决方案

**1️⃣ GTX索引构建失败**
```bash
# 检查GTX二进制文件
ls -la /share/org/YZWL/yzwl_lixg/software/gtx/bin/gtx
/gtx/bin/gtx --version

# 检查参考基因组
samtools faidx reference.fa
head -5 reference.fa

# 强制重新构建索引
rm -f reference.fa.gtx*
biopytools fastq2vcf-gtx -i ./raw -r ./ref.fa -p ./proj --step 1
```

**2️⃣ GTX WGS比对失败**
```bash
# 检查GTX索引文件
ls -la reference.fa.gtx*

# 检查质控结果
ls -la ./project/01.data/clean/

# 重新执行比对步骤
biopytools fastq2vcf-gtx -i ./raw -r ./ref.fa -p ./proj --step 3
```

**3️⃣ 内存不足问题**
```bash
# 调整窗口大小减少内存使用
biopytools fastq2vcf-gtx \
    -i ./raw -r ./ref.fa -p ./proj \
    --gtx-window-size 5000000  # 5MB窗口

# 减少线程数
biopytools fastq2vcf-gtx \
    -i ./raw -r ./ref.fa -p ./proj \
    --threads-gtx 32
```

**4️⃣ 临时文件空间不足**
```bash
# 检查磁盘空间
df -h ./project/

# 清理临时文件
rm -rf ./project/.tmp/

# 使用自定义临时目录
mkdir -p /tmp/gtx_tmp
# GTX会自动使用项目目录下的.tmp/
```

**5️⃣ faketime相关问题**
```bash
# 检查faketime是否安装
which faketime
faketime '2020-10-20 00:00:00' date

# 手动测试GTX索引构建
faketime '2020-10-20 00:00:00' /path/to/gtx index reference.fa
```

## 🔧 GTX特性详解 | GTX Features Details

### GTX索引系统
- **快速访问**: 专为快速序列访问设计
- **内存映射**: 支持大基因组的内存映射访问
- **faketime兼容**: 使用固定时间戳确保可重复性

### GTX WGS流程
- **一体化处理**: 比对和变异检测在单一步骤中完成
- **内存优化**: 智能的内存管理和垃圾回收
- **分块处理**: 自动将大数据集分成可处理的块

### GTX集群模式
- **作业分割**: 自动将大数据集分割成集群作业
- **负载均衡**: 智能的任务分配和负载均衡
- **结果合并**: 自动合并集群作业结果

## 📚 相关资源 | Resources

- **GTX Toolkit**: https://github.com/BGI-Genomics/GTX
- **GATK Best Practices**: https://gatk.broadinstitute.org/hc/en-us
- **BWA-MEM**: http://bio-bwa.sourceforge.net/bwa.shtml
- **SAMtools**: http://www.htslib.org/doc/samtools.html
- **faketime**: https://github.com/wolfcw/libfaketime

## ⚠️ 注意事项 | Important Notes

- **GTX依赖**: 需要正确安装GTX工具包
- **faketime要求**: GTX索引构建需要faketime工具
- **内存规划**: 根据窗口大小规划内存使用
- **存储空间**: 预留充足的磁盘空间（输入数据×5-7倍）
- **文件权限**: 确保对GTX二进制文件有执行权限
- **时间兼容**: GTX使用固定时间戳确保可重复性

## 🏆 最佳实践 | Best Practices

1. **环境准备**: 确保GTX和faketime正确安装
2. **参数调优**: 根据硬件配置调整线程和窗口大小
3. **内存监控**: 监控内存使用情况，必要时调整窗口大小
4. **索引检查**: 定期检查GTX索引文件完整性
5. **日志分析**: 关注详细日志中的性能指标
6. **备份策略**: 定期备份GTX索引和中间结果

---

🔬 **CPU环境下的高性能基因组分析 | High-performance genomic analysis in CPU environments**