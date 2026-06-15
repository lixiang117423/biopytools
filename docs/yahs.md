# YaHS Hi-C Scaffolding 流程模块

**专业的Hi-C scaffolding分析工具 | Professional Hi-C Scaffolding Analysis Tool**

## 功能概述 | Overview

YaHS模块是基于YaHS软件构建的Hi-C scaffolding分析流程，提供从基因组索引构建、Hi-C数据比对、scaffolding到质量评估的完整流程。支持断点续传、单步执行和灵活的参数配置，适用于各种基因组scaffolding分析研究。

## 主要特性 | Key Features

- **完整流程支持**: 索引→比对→挂载→热图→JBAT→评估六步自动化流程
- **断点续传**: 自动检测已完成步骤，支持中断后继续执行
- **单步执行**: 支持单独运行任意步骤或完整流程
- **工具自动检测**: BWA、samtools等常用工具自动检测可用性
- **灵活参数配置**: 支持YaHS所有核心参数和可选参数
- **JBAT支持**: 生成Juicebox JBAT手动校正所需文件
- **质量评估**: 自动计算N50、N90、L50、L90等统计指标
- **详细日志**: 完整的运行日志和错误追踪

## 快速开始 | Quick Start

### 基本用法 | Basic Usage

```bash
# 运行完整流程
biopytools yahs \
    -r reference.fa \
    -1 hic_R1.fq.gz \
    -2 hic_R2.fq.gz

# 自定义酶切位点和线程数
biopytools yahs \
    -r reference.fa \
    -1 hic_R1.fq.gz \
    -2 hic_R2.fq.gz \
    -e GATC \
    -t 24

# 只运行YaHS挂载步骤
biopytools yahs \
    -r reference.fa \
    -1 hic_R1.fq.gz \
    -2 hic_R2.fq.gz \
    -s 3
```

## 参数说明 | Parameters

### 必需参数 | Required Parameters

| 参数 | 描述 | 示例 |
|------|------|------|
| `-r, --ref` | 参考基因组FASTA文件 | `-r genome.fa` |
| `-1, --hic-r1` | Hi-C R1测序文件 | `-1 hic_R1.fq.gz` |
| `-2, --hic-r2` | Hi-C R2测序文件 | `-2 hic_R2.fq.gz` |

### 输出配置 | Output Configuration

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-o, --output-dir` | `./yahs_output` | 输出目录路径 |

### 资源配置 | Resource Configuration

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-t, --threads` | `12` | 线程数 |
| `--java-ram` | `32G` | Java内存 |
| `--sam-ram` | `4G` | Samtools排序内存 |

### YaHS 核心参数 | YaHS Core Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-e, --enzyme` | `GATC` | 限制性酶切位点序列 |
| `--min-len` | `10000` | 最小contig长度 |
| `--min-mapq` | `30` | 最小MAPQ值 |
| `--no-contig-ec` | `False` | 跳过contig错误校正 |
| `--no-scaffold-ec` | `False` | 跳过scaffold错误校正 |
| `--resolutions` | `None` | 分辨率列表(逗号分隔) |
| `--rounds` | `1` | 每分辨率运行轮数 |
| `--telo-motif` | `None` | 端粒序列模体 |

### 工具路径 | Tool Paths

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--yahs-bin` | `yahs` | YaHS可执行文件路径 |
| `--juicer-bin` | `juicer` | juicer可执行文件路径 |
| `--juicer-jar` | `None` | juicer_tools.jar文件路径 |
| `--bwa-bin` | `bwa` | BWA可执行文件路径 |
| `--samtools-bin` | `samtools` | samtools可执行文件路径 |
| `--java-cmd` | `java` | Java可执行文件路径 |

### 执行控制 | Execution Control

| 参数 | 描述 |
|------|------|
| `-s, --step` | 运行指定步骤 (1-6) |
| `--force-rerun` | 强制重新运行所有步骤 |
| `--keep-temp` | 保留临时文件 |

## 流程步骤 | Pipeline Steps

| 步骤 | 名称 | 描述 |
|------|------|------|
| **1** | 构建基因组索引 | 为参考基因组构建BWA和SAMtools索引 |
| **2** | Hi-C数据比对 | 使用BWA MEM进行比对，排序并标记重复 |
| **3** | YaHS染色体挂载 | 使用YaHS进行Hi-C scaffolding |
| **4** | 生成标准Hi-C热图 | �认.hic文件用于可视化 |
| **5** | 生成JBAT文件 | 生成Juicebox JBAT手动校正文件 |
| **6** | 组装质量评估 | 计算N50、N90等统计指标 |

## 输出目录结构 | Output Directory Structure

```
yahs_output/
├── 00_pipeline_info/          # 流程元数据
│   └── software_versions.yml
├── 01_indexing/               # 步骤1：索引
│   └── genome.fa.*
├── 02_mapping/                # 步骤2：比对
│   ├── aligned_sorted_dedup.bam
│   └── aligned_sorted_dedup.bam.bai
├── 03_scaffolding/            # 步骤3：YaHS挂载
│   ├── yahs_out_scaffolds_final.fa
│   ├── yahs_out_scaffolds_final.agp
│   └── yahs_out.bin
├── 04_hic_standard/           # 步骤4：标准Hi-C热图
│   └── yahs_out_final.hic
├── 05_jbat/                   # 步骤5：JBAT文件
│   ├── out_JBAT.hic
│   ├── out_JBAT.assembly
│   └── out_JBAT.liftover.agp
├── 06_assessment/             # 步骤6：质量评估
│   └── assembly_metrics.txt
└── 99_logs/                   # 日志文件
    └── yahs_pipeline.log
```

## 使用示例 | Usage Examples

### 示例1：基本流程 | Example 1: Basic Pipeline

```bash
# 使用默认参数运行完整流程
biopytools yahs \
    -r reference.fa \
    -1 hic_R1.fq.gz \
    -2 hic_R2.fq.gz
```

### 示例2：自定义参数 | Example 2: Custom Parameters

```bash
# 使用自定义酶切位点、最小长度和线程数
biopytools yahs \
    -r reference.fa \
    -1 hic_R1.fq.gz \
    -2 hic_R2.fq.gz \
    -e GANTC \
    --min-len 50000 \
    -t 24
```

### 示例3：单步执行 | Example 3: Single Step Execution

```bash
# 只运行Hi-C比对步骤（需要已完成索引）
biopytools yahs \
    -r reference.fa \
    -1 hic_R1.fq.gz \
    -2 hic_R2.fq.gz \
    -s 2

# 只运行YaHS挂载步骤（需要已完成比对）
biopytools yahs \
    -r reference.fa \
    -1 hic_R1.fq.gz \
    -2 hic_R2.fq.gz \
    -s 3
```

### 示例4：指定juicer_tools | Example 4: Specify juicer_tools

```bash
# 指定juicer_tools.jar路径以生成热图
biopytools yahs \
    -r reference.fa \
    -1 hic_R1.fq.gz \
    -2 hic_R2.fq.gz \
    --juicer-jar /path/to/juicer_tools.jar
```

### 示例5：高级YaHS参数 | Example 5: Advanced YaHS Parameters

```bash
# 使用高级YaHS参数
biopytools yahs \
    -r reference.fa \
    -1 hic_R1.fq.gz \
    -2 hic_R2.fq.gz \
    --resolutions 10000,50000,100000,500000,1000000 \
    --rounds 2 \
    --no-contig-ec \
    --telo-motif TTAGGG
```

## 系统要求 | System Requirements

### 依赖软件 | Dependencies

- **YaHS** (v1.2.2或更新)
  - 下载地址: https://github.com/c-zhou/yahs
- **BWA** (v0.7.17或更新)
- **samtools** (v1.10或更新)
- **juicer** (YaHS自带)
- **juicer_tools.jar** (可选，用于生成热图)
- **Java** (v8或更新，如果使用juicer_tools)

### 安装依赖 | Installing Dependencies

```bash
# 使用conda安装YaHS
conda create -n yahs_v.1.2.2 -c bioconda yahs

# 安装BWA和samtools
conda install -c bioconda bwa samtools

# 或使用系统包管理器
# Ubuntu/Debian
sudo apt-get install bwa samtools

# CentOS/RHEL
sudo yum install bwa samtools
```

### 硬件建议 | Hardware Recommendations

- **CPU**: 多核处理器（推荐8核以上）
- **RAM**: 最少16GB（大基因组推荐64GB以上）
- **存储**: 预留Hi-C数据大小5倍的磁盘空间
- **运行时间**: 根据基因组大小和测序深度，通常2-12小时

## 注意事项 | Important Notes

1. **输入文件格式**:
   - 参考基因组必须是FASTA格式
   - Hi-C数据支持FASTQ/GZIP压缩格式
   - 参考基因组需要预先建立索引（或由流程自动构建）

2. **酶切位点选择**:
   - DpnII: `GATC`
   - MboI: `GATC`
   - Arima 2-酶: `GATC,GANTC`
   - Arima 4-酶: `GATC,GANTC,CTNAG,TTAA`

3. **内存使用**:
   - 大基因组（>1Gb）建议增加sam-ram参数
   - Java内存建议根据系统总内存调整（50-60%）

4. **断点续传**:
   - 默认启用断点续传，已完成步骤会自动跳过
   - 使用`--force-rerun`强制重新运行所有步骤

5. **JBAT文件**:
   - 需要指定`--juicer-jar`参数
   - 生成的文件可用于Juicebox JBAT手动校正

## 故障排除 | Troubleshooting

### 常见问题 | Common Issues

**Q: "BWA not found" 错误**
```bash
# 确保BWA在PATH中或使用--bwa-bin参数
which bwa
biopytools yahs ... --bwa-bin /path/to/bwa
```

**Q: 内存不足错误**
```bash
# 减少线程数或增加内存分配
biopytools yahs ... -t 8 --sam-ram 2G
```

**Q: YaHS运行失败**
```bash
# 检查YaHS版本
yahs --version

# 使用完整路径
biopytools yahs ... --yahs-bin /path/to/yahs
```

**Q: juicer_tools未找到**
```bash
# 指定完整路径
biopytools yahs ... --juicer-jar /share/org/YZWL/yzwl_lixg/software/juicer/scripts/common/juicer_tools.jar
```

## 相关资源 | Related Resources

- [YaHS官方文档](https://github.com/c-zhou/yahs)
- [YaHS论文](https://academic.oup.com/bioinformatics/article/39/1/btac808/6782626)
- [Juicebox文档](https://github.com/aidenlab/Juicebox)
- [Hi-C数据分析最佳实践](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6606163/)

## 引用信息 | Citation

如果在学术研究中使用此工具，请引用YaHS相关文献：

```
Zhou, C., McCarthy, S. A., & Durbin, R. (2023).
YaHS: yet another Hi-C scaffolding tool.
Bioinformatics, 39(1), btac808.
```
