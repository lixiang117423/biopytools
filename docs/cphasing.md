# 🧬 CPhasing 基因组分相和挂载工具

**基于Hi-C数据的多倍体基因组分相和染色体级别挂载 | Polyploid Genome Phasing and Chromosome-level Scaffolding with Hi-C Data**

## 📖 功能概述 | Overview

CPhasing是一个专门用于多倍体基因组分相和挂载的工具，基于Hi-C、Pore-C、HiFi-C等染色体构象捕获数据。本模块提供了CPhasing pipeline的便捷封装，支持Hi-C数据的完整分析流程。

### ✨ 主要特性 | Key Features

- **🔄 完整分析流程**: 从序列比对到分相挂载的一站式流程
- **🎯 多种模式支持**: 支持单倍体(haploid)、多倍体(phasing)、基础模式(basal)
- **🛡️ 灵活参数配置**: 丰富的参数选项，满足不同分析需求
- **📊 自动分组**: 支持自动或指定分组数的多轮分相
- **🚀 高效处理**: 优化的处理流程，支持大规模基因组
- **⚙️ 多种预设**: precision、sensitive、very-sensitive、nofilter四种分析预设

## 🚀 快速开始 | Quick Start

### 基本用法 | Basic Usage

```bash
# 最简单的使用方式
biopytools cphasing \
    -f genome.fa \
    -hic1 hic_R1.fq.gz \
    -hic2 hic_R2.fq.gz \
    -t 10 \
    -n 8:4
```

### 高级用法 | Advanced Usage

```bash
# 使用敏感模式进行分相
biopytools cphasing \
    -f genome.fa \
    -hic1 hic_R1.fq.gz \
    -hic2 hic_R2.fq.gz \
    -t 20 \
    -n 0 \
    --preset sensitive \
    -o results
```

## 📋 参数说明 | Parameters

### 必需参数 | Required Parameters

| 参数 | 描述 | 示例 |
|------|------|------|
| `-f, --fasta` | 基因组FASTA文件 | `-f genome.fa` |
| `--hic1` | Hi-C R1 reads文件 | `--hic1 R1.fq.gz` |
| `--hic2` | Hi-C R2 reads文件 | `--hic2 R2.fq.gz` |

### 核心参数 | Core Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-n, --groups` | `0` | 🎯 分组数 (如: `8:4` 表示第一轮8组，第二轮4组；`0` 表示自动) |
| `-t, --threads` | `10` | 🔧 线程数 |

### 模式和预设 | Mode and Preset

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--mode` | `phasing` | 📊 分相模式：<br>- `phasing`: 多倍体分相 (默认)<br>- `haploid`: 单倍体模式<br>- `basal`: 基础模式<br>- `basal_withprune`: 基础模式+修剪 |
| `--preset` | `precision` | ⚙️ 分析预设：<br>- `precision`: 精确模式 (默认)<br>- `sensitive`: 敏感模式<br>- `very-sensitive`: 超敏感模式<br>- `nofilter`: 无过滤模式 |

### 输出控制 | Output Control

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-o, --output-dir` | `./cphasing_output` | 📁 输出目录 |

### 步骤控制 | Step Control

| 参数 | 描述 |
|------|------|
| `--steps` | 🎯 运行指定步骤 (如: `1,2,3`) |
| `--skip-steps` | ⏭️ 跳过步骤 (如: `1,2`) |

### Hi-C Mapper 参数 | Hi-C Mapper Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--hic-aligner` | `_chromap` | 🔧 Hi-C比对器 (`_chromap`, `chromap`, `minimap2`, `bwa-mem2`) |
| `--hic-mapper-k` | `自动` | kmer大小 |
| `--hic-mapper-w` | `自动` | 窗口大小 |
| `--mapping-quality` | `0` | 最小比对质量 |

### HCR 参数 | HCR Parameters

| 参数 | 描述 |
|------|------|
| `--hcr` | 启用高置信区域过滤 |
| `--pattern` | 酶切位点模式 (如: `AAGCTT`) |

### 高级选项 | Advanced Options

| 参数 | 描述 |
|------|------|
| `--low-memory` | 📉 启用低内存模式 |

## 💡 使用示例 | Usage Examples

### 示例1: 四倍体基因组分相 | Example 1: Tetraploid Genome Phasing

```bash
# 四倍体基因组，8条基础染色体，第一轮分8组，第二轮分4组
biopytools cphasing \
    -f tetraploid_genome.fa \
    -hic1 hic_R1.fq.gz \
    -hic2 hic_R2.fq.gz \
    -t 20 \
    -n 8:4 \
    -o tetraploid_results
```

### 示例2: 单倍体基因组挂载 | Example 2: Haploid Genome Scaffolding

```bash
# 单倍体模式，用于二倍体或单倍体基因组挂载
biopytools cphasing \
    -f diploid_genome.fa \
    -hic1 hic_R1.fq.gz \
    -hic2 hic_R2.fq.gz \
    -t 16 \
    --mode haploid \
    -n 10 \
    -o haploid_results
```

### 示例3: 使用敏感模式 | Example 3: Using Sensitive Preset

```bash
# 使用敏感模式，聚类更多低信号contigs
biopytools cphasing \
    -f genome.fa \
    -hic1 hic_R1.fq.gz \
    -hic2 hic_R2.fq.gz \
    --preset sensitive \
    -n 0 \
    -t 24
```

### 示例4: 启用HCR过滤 | Example 4: Enable HCR Filtering

```bash
# 使用高置信区域过滤提高分相质量
biopytools cphasing \
    -f genome.fa \
    -hic1 hic_R1.fq.gz \
    -hic2 hic_R2.fq.gz \
    --hcr \
    --pattern AAGCTT \
    -t 16 \
    -n 8:4
```

### 示例5: 只运行特定步骤 | Example 5: Run Specific Steps Only

```bash
# 只运行步骤3 (hyperpartition)
biopytools cphasing \
    -f genome.fa \
    -hic1 hic_R1.fq.gz \
    -hic2 hic_R2.fq.gz \
    --steps 3 \
    -n 8:4
```

### 示例6: 跳过某些步骤 | Example 6: Skip Specific Steps

```bash
# 跳过步骤1和2 (alleles和prepare)
biopytools cphasing \
    -f genome.fa \
    -hic1 hic_R1.fq.gz \
    -hic2 hic_R2.fq.gz \
    --skip-steps 1,2 \
    -n 8:4
```

## 📊 输出结果 | Output Results

### 输出目录结构 | Output Directory Structure

```
cphasing_output/
├── 0.mapper/                  # 序列比对结果
│   ├── aligned.bam
│   └── alignment.log
├── 1.alleles/                 # 等位基因识别
│   ├── alleles.txt
│   └── alleles.log
├── 2.prepare/                 # 数据准备
│   ├── contact_matrix.txt
│   └── prepare.log
├── 3.hyperpartition/          # 超图分区
│   ├── first.clusters.txt
│   ├── second.clusters.txt
│   └── partition.log
├── 4.scaffolding/             # 挂载
│   ├── groups.agp
│   ├── groups.asm.fasta
│   └── scaffolding.log
├── 5.plot/                    # 热图绘制
│   ├── heatmap.png
│   └── plot.log
└── cphasing_pipeline.log      # 完整日志
```

### 主要输出文件 | Main Output Files

| 文件 | 描述 |
|------|------|
| `groups.agp` | 挂载结果AGP格式文件 |
| `groups.asm.fasta` | 挂载后的基因组序列 |
| `first.clusters.txt` | 第一轮分组结果 |
| `second.clusters.txt` | 第二轮分组结果 |
| `heatmap.png` | Hi-C接触矩阵热图 |

## ⚙️ 系统要求 | System Requirements

### 依赖软件 | Dependencies

- **CPhasing** (版本 0.2.10 或更新)
  - conda环境: `cphasing_v.0.2.10`
- **Python** (版本 3.7+)
- **Python包**:
  - `click` - 命令行界面
  - `dataclasses` - 数据类 (Python 3.7+)

### 安装CPhasing | Installing CPhasing

```bash
# 下载CPhasing
cd ~/software/cphasing
wget https://github.com/wangyibin/CPhasing/releases/download/v0.2.10/CPhasing_v0.2.10.linux-x86.tar.gz

# 解压
tar -xzf CPhasing_v0.2.10.linux-x86.tar.gz
cd CPhasing_v0.2.10

# 创建conda环境
conda env create -f environment.yml

# 激活环境
conda activate cphasing_v.0.2.10
```

### 硬件建议 | Hardware Recommendations

- **CPU**: 多核处理器 (推荐16核以上)
- **RAM**: 最少32GB (大基因组推荐128GB以上)
- **存储**: 预留基因组文件大小10倍的磁盘空间

## ⚠️ 注意事项 | Important Notes

1. **CPhasing环境**: 确保已安装并配置好CPhasing conda环境 (`cphasing_v.0.2.10`)
2. **Hi-C数据质量**: Hi-C数据质量直接影响分相效果
3. **分组数设置**:
   - 二倍体基因组: `-n 0` (自动) 或 `-n 倍数:0`
   - 四倍体基因组: `-n 8:4` 或 `-n 0`
   - 六倍体基因组: `-n 0` (推荐自动)
4. **内存使用**: 大基因组需要大量内存，建议使用 `--low-memory` 模式
5. **conda自动包装**: 本工具会自动使用 `conda run --no-capture-output` 包装命令

## 🐛 故障排除 | Troubleshooting

### 常见问题 | Common Issues

**Q: "cphasing命令未找到" 错误**
```bash
# 检查conda环境
conda env list | grep cphasing

# 激活环境
conda activate cphasing_v.0.2.10

# 验证cphasing命令
cphasing --help
```

**Q: 内存不足错误**
```bash
# 启用低内存模式
biopytools cphasing ... --low-memory

# 或者减少线程数
biopytools cphasing ... -t 8
```

**Q: Hi-C比对失败**
```bash
# 检查Hi-C reads文件
zcat hic_R1.fq.gz | head -n 2

# 尝试不同的比对器
biopytools cphasing ... --hic-aligner minimap2
```

**Q: 分组效果不好**
```bash
# 使用敏感模式
biopytools cphasing ... --preset sensitive

# 启用HCR过滤
biopytools cphasing ... --hcr --pattern AAGCTT

# 调整分组数
biopytools cphasing ... -n 0  # 自动分组
```

## 📚 相关资源 | Related Resources

- [CPhasing官方文档](https://wangyibin.github.io/CPhasing)
- [CPhasing GitHub仓库](https://github.com/wangyibin/CPhasing)
- [Hi-C数据分析最佳实践](https://www.nature.com/articles/s41588-021-00951-w)

## 📄 许可证 | License

本项目采用MIT许可证 - 详见 [LICENSE](LICENSE) 文件

**注意**: CPhasing软件本身有单独的许可证，请查看官方文档。

## 🔬 引用信息 | Citation

如果在学术研究中使用CPhasing，请引用：

```
Yibin Wang, Ping Zhao, Xiaofei Zeng, et al.
Enhanced Pore-C with C-Phasing Enables Chromosomal-Scale, Haplotype-Resolved Assembly of Ultra-Complex Genomes.
PREPRINT (2025). https://doi.org/10.21203/rs.3.rs-7343323/v1
```

---

**📧 技术支持 | Technical Support**: 如有问题，请提交 [Issue](https://github.com/your-org/biopytools/issues)
