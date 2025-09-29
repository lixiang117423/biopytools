# 🧬 OrthoFinder泛基因组分析工具

[![Python](https://img.shields.io/badge/Python-3.6+-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Platform](https://img.shields.io/badge/Platform-Linux%20%7C%20macOS-lightgrey.svg)]()

一个完整的泛基因组分析流程工具，使用OrthoFinder进行同源基因群聚类，然后根据基因在不同基因组中的分布模式将基因家族分类为核心基因、软核心基因、附属基因和特异基因，适用于比较基因组学和进化分析研究。

## 📋 目录

- [功能特点](#-功能特点)
- [安装](#-安装)
- [快速开始](#-快速开始)
- [使用方法](#-使用方法)
- [参数说明](#-参数说明)
- [输入文件要求](#-输入文件要求)
- [输出文件说明](#-输出文件说明)
- [使用示例](#-使用示例)
- [系统要求](#-系统要求)
- [泛基因组分类标准](#-泛基因组分类标准)
- [应用场景](#-应用场景)
- [故障排除](#-故障排除)
- [最佳实践](#-最佳实践)
- [引用](#-引用)

## ✨ 功能特点

- 🤖 **全自动分析**：完整的OrthoFinder同源基因群分析流程
- 🎯 **智能分类**：基于分布模式的泛基因组基因家族自动分类
- 🔍 **多算法支持**：支持BLAST、Diamond、MMseqs2等多种序列搜索算法
- ⚙️ **灵活配置**：可自定义基因家族阈值和分析参数
- 📊 **详细报告**：生成综合统计报告和详细数据表
- 🚀 **性能优化**：针对大规模基因组数据进行优化
- 🌳 **系统发育分析**：可选的系统发育树构建功能
- 🧬 **双模式**：同时支持蛋白质和DNA序列分析

## 📦 安装

### 依赖软件

在使用本工具前，请确保已安装以下依赖：

1. **必需软件**
   - Python 3.6+
   - OrthoFinder (v2.3.0+)
   - BLAST+ 或 Diamond
   - MCL (Markov Cluster Algorithm)

2. **可选软件**（用于系统发育分析）
   - FastTree / RAxML / IQ-TREE
   - MAFFT / MUSCLE

## 🚀 快速开始

```bash
# 基本用法
biopytools orthofinder -i protein_sequences/ -o pangenome_results

# 查看帮助
biopytools orthofinder --help
```

## 💻 使用方法

### 基本命令格式

```bash
biopytools orthofinder [OPTIONS]
```

### 必需参数

- `-i, --input`: 输入蛋白质序列文件目录
- `-o, --output`: 输出结果目录

## ⚙️ 参数说明

### 核心参数

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `-i, --input` | 路径 | 必需 | 输入蛋白质序列文件目录 |
| `-o, --output` | 路径 | 必需 | 输出目录 |
| `-n, --project-name` | 字符串 | - | 项目名称 |
| `--soft-threshold` | 字符串/数值 | all-1 | Soft基因家族阈值 |

### 性能参数

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `-t, --threads` | 整数 | 88 | 并行计算线程数 |
| `-s, --search` | 选项 | blastp | 序列搜索程序 |

**搜索算法选项**：
- `blast` / `blastp`: 标准BLAST蛋白质搜索
- `diamond`: 快速Diamond搜索
- `diamond_ultra_sens`: Diamond高敏感模式
- `mmseqs`: MMseqs2快速搜索
- `blast_nucl`: 核酸序列BLAST（DNA模式）

### OrthoFinder参数

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `--mcl-inflation` | 浮点数 | 1.2 | MCL聚类膨胀参数 |
| `--basic-only` | 布尔 | True | 仅进行基础分析 |
| `--orthofinder-path` | 路径 | orthofinder | OrthoFinder程序路径 |

### 系统发育参数

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `--generate-trees` | 布尔 | False | 生成系统发育树 |
| `--msa-program` | 选项 | mafft | 多序列比对程序 |
| `--tree-program` | 选项 | fasttree | 系统发育树构建程序 |

### 其他参数

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `-d, --dna` | 布尔 | False | 输入序列为DNA |
| `--resume` | 布尔 | True | 使用已有结果续跑 |
| `--force` | 布尔 | False | 强制重新分析 |
| `--skip-orthofinder` | 布尔 | False | 跳过OrthoFinder步骤 |

### Soft阈值说明

`--soft-threshold` 参数定义了软核心基因的判定标准：

- `all-1`: 存在于 N-1 个基因组（默认）
- `all-2`: 存在于 N-2 个基因组
- `specific`: 根据数据自动确定
- `数值`: 直接指定基因组数量阈值（如 5、8、10）

## 📁 输入文件要求

### 蛋白质序列文件

- **格式**：标准FASTA格式 (`.fa`, `.faa`, `.fasta`)
- **组织**：每个基因组一个独立文件
- **命名**：文件名将作为基因组标识符
- **序列ID**：需要在同一基因组内唯一

### DNA序列文件（使用--dna时）

- **格式**：标准FASTA格式核酸序列
- **内容**：通常为CDS序列或基因序列
- **要求**：需要正确的开放阅读框

### 目录结构示例

```
protein_sequences/
├── genome1.faa
├── genome2.faa
├── genome3.faa
└── genome4.faa
```

### 序列格式示例

```fasta
>gene001 hypothetical protein
MKILVFASLLSLLAAGVQAAPEAQPVIKLEEATGKGQRWLWAGLASRLVDAKMQTIHSASLVS
>gene002 DNA polymerase
MAKTFEKVKLAASQAGEEAATQSNAQLQMKLVLKATTQADQAIKAKTQASLAALNQADEQLAS
```

## 📄 输出文件说明

### 核心结果文件

```
output_directory/
├── pangenome_classification.txt      # 基因家族分类结果
├── gene_families_table.txt           # 详细基因家族数据表
├── comprehensive_report.txt          # 综合分析报告
├── orthofinder_results/              # OrthoFinder原始输出
│   ├── Orthogroups/
│   │   ├── Orthogroups.txt
│   │   └── Orthogroups.GeneCount.txt
│   ├── Phylogenetic_Hierarchical_Orthogroups/
│   └── Species_Tree/
└── analysis.log                       # 分析日志
```

### 输出文件说明

- **pangenome_classification.txt**: 每个基因家族的分类结果（核心/软核心/附属/特异）
- **gene_families_table.txt**: 包含每个基因家族在各基因组中的基因数量
- **comprehensive_report.txt**: 包含统计摘要和分析细节
- **orthofinder_results/**: OrthoFinder的完整输出结果

## 📝 使用示例

### 1. 基本泛基因组分析

```bash
biopytools orthofinder -i protein_sequences/ -o pangenome_results
```

### 2. 指定项目名称和自定义soft阈值

```bash
biopytools orthofinder \
    -i data/proteomes/ \
    -o results/ \
    --project-name "Ecoli_pangenome" \
    --soft-threshold all-2
```

### 3. 使用Diamond加速和高线程数

```bash
biopytools orthofinder \
    -i sequences/ \
    -o output/ \
    -s diamond \
    -t 64 \
    --soft-threshold 5
```

### 4. 数值型soft阈值设置

```bash
biopytools orthofinder \
    -i proteomes/ \
    -o results/ \
    --soft-threshold 8 \
    --search diamond_ultra_sens
```

### 5. 包含系统发育树的完整分析

```bash
biopytools orthofinder \
    -i sequences/ \
    -o comprehensive_analysis/ \
    --generate-trees \
    --msa-program muscle \
    --tree-program iqtree \
    --search diamond \
    -t 96
```

### 6. DNA序列分析模式

```bash
biopytools orthofinder \
    -i dna_sequences/ \
    -o dna_results/ \
    --dna \
    --search blast_nucl \
    --mcl-inflation 1.5
```

### 7. 高性能服务器优化配置

```bash
biopytools orthofinder \
    -i large_dataset/ \
    -o hpc_results/ \
    -t 128 \
    --search mmseqs \
    --mcl-inflation 1.1 \
    --soft-threshold all-3 \
    --project-name "HighThroughput_Analysis"
```

## 🖥️ 系统要求

### 硬件建议

| 数据规模 | 基因组数 | 内存 | 存储 | CPU |
|---------|---------|------|------|-----|
| 小规模 | < 50 | 4GB | 20GB | 4核+ |
| 中规模 | 50-200 | 16GB | 100GB | 16核+ |
| 大规模 | 200-500 | 64GB | 500GB | 32核+ |
| 超大规模 | > 500 | 128GB+ | 1TB+ | 64核+ |

### 软件版本要求

- Python: 3.6+
- OrthoFinder: 2.3.0+
- BLAST+: 2.2.31+ 或 Diamond: 0.8.0+
- MCL: 14-137+

### 性能估算

- **时间复杂度**: O(n²m)（n=基因组数，m=平均基因数）
- **空间复杂度**: O(nm)
- 使用Diamond可显著提升速度（约10-100倍）

## 🎯 泛基因组分类标准

本工具将基因家族分为四类：

| 类别 | 标准 | 特征 |
|------|------|------|
| 🔴 **核心基因** (Core) | 存在于所有基因组 | 基础代谢和维持功能 |
| 🟠 **软核心基因** (Soft-core) | 存在于大多数基因组 | 重要但非必需功能 |
| 🟡 **附属基因** (Accessory) | 存在于部分基因组 | 环境适应相关 |
| 🟢 **特异基因** (Specific) | 存在于少数基因组 | 特殊功能或水平转移 |

## 🔬 应用场景

- 🦠 细菌和古菌泛基因组分析
- 🧬 物种进化和基因组可塑性研究
- 🔴 核心功能基因识别
- 📊 基因组比较和系统发育分析
- 🦠 病原菌毒力因子挖掘
- 🏭 工业微生物功能基因发现

## 🔧 故障排除

### 常见问题

#### 1. OrthoFinder not found
```bash
# 检查OrthoFinder是否在PATH中
which orthofinder

# 如果不在PATH中，使用--orthofinder-path指定路径
biopytools orthofinder -i input/ -o output/ --orthofinder-path /path/to/orthofinder
```

#### 2. Memory error
- 减少线程数：`-t 32`（降低并行度）
- 使用Diamond代替BLAST：`-s diamond`
- 增加系统交换空间

#### 3. Empty orthogroups
- 检查输入序列格式是否正确
- 确认序列文件不为空
- 验证FASTA格式规范性

#### 4. MCL inflation error
- 调整MCL参数：`--mcl-inflation 1.5`
- 尝试不同的inflation值（0.5-5.0）

#### 5. Tree generation failed
- 检查MSA软件是否安装
- 确认树构建软件可用
- 尝试其他树构建方法

### 优化建议

- ⚡ 大数据集优先使用Diamond搜索
- 🧵 根据内存限制调整线程数（一般为CPU核心数）
- ✅ 使用`--resume`参数进行断点续跑
- 💿 使用SSD存储提高I/O性能
- 📊 使用`htop`或`top`监控系统资源

## 🏆 最佳实践

### 数据准备

1. ✅ 使用高质量的基因组注释
2. 🔍 确保蛋白质序列完整性
3. 🏷️ 统一序列命名规范
4. 🧹 移除冗余和低质量序列

### 参数选择

1. **小数据集**（< 50基因组）
   - 使用BLAST获得最佳精度
   - 可以使用较高线程数

2. **大数据集**（> 200基因组）
   - 使用Diamond平衡速度与精度
   - 考虑使用MMseqs2
   - 调整MCL参数优化聚类

3. **Soft阈值**
   - 近缘物种：使用`all-1`
   - 远缘物种：使用`all-2`或更低
   - 多样性高：使用数值型阈值

### 结果解释

- 🔴 **核心基因**：通常为基础代谢功能，高度保守
- 🟡 **附属基因**：可能与环境适应、毒力相关
- 🟢 **特异基因**：可能包含新功能或水平基因转移
- 📚 结合功能注释（GO、KEGG）进行生物学解释

### 质量控制

1. 🔍 检查基因组完整性（BUSCO评估）
2. ✅ 验证同源基因群的合理性
3. 🔄 比较不同参数设置的结果
4. 🧪 使用已知基因验证分类准确性

## 📚 引用

如果在学术研究中使用此工具，请引用相关文献：

### OrthoFinder
```
Emms, D.M., Kelly, S. OrthoFinder: phylogenetic orthology inference for comparative genomics.
Genome Biology 20, 238 (2019). https://doi.org/10.1186/s13059-019-1832-y
```

### Diamond
```
Buchfink, B., Xie, C. & Huson, D. Fast and sensitive protein alignment using DIAMOND.
Nature Methods 12, 59–60 (2015). https://doi.org/10.1038/nmeth.3176
```

### MCL Algorithm
```
Van Dongen, S. Graph Clustering by Flow Simulation.
PhD thesis, University of Utrecht (2000).
```

### 泛基因组概念
```
Tettelin, H. et al. Genome analysis of multiple pathogenic isolates of Streptococcus agalactiae: 
implications for the microbial "pan-genome". PNAS 102, 13950-13955 (2005).
```

## 📧 联系方式

- **Issues**: [GitHub Issues](https://github.com/yourusername/biopytools/issues)
- **Email**: your.email@example.com
- **Documentation**: [完整文档](https://biopytools.readthedocs.io)

## 📜 许可证

本项目采用 MIT 许可证 - 详见 [LICENSE](LICENSE) 文件

---

**⭐ 如果这个工具对您有帮助，请给我们一个Star！**