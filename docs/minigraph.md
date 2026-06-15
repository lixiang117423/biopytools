# Minigraph 泛基因组图构建和分析模块

**泛基因组图构建、序列映射、SV调用和bubble提取 | Pangenome Graph Construction, Sequence Mapping, SV Calling and Bubble Extraction**

## 功能概述 | Overview

Minigraph 模块是一个高效的泛基因组图构建和分析工具包，基于 minigraph 算法实现。它可以从多个基因组组装构建泛基因组图，支持序列到图的映射、结构变异调用和bubble提取，是泛基因组学研究的核心工具。

## 主要特性 | Key Features

- **高效图构建**: 基于minimap2算法，快速构建泛基因组图
- **增量式构建**: 支持逐步添加样本到现有图中
- **序列映射**: 支持多种测序类型的序列到图映射
- **SV调用**: 为每个样本生成基于图的SV调用结果
- **Bubble提取**: 从图中提取所有结构变异区域
- **灵活配置**: 支持多种预设和参数调整

## 快速开始 | Quick Start

### 查看所有子命令 | View All Subcommands

```bash
biopytools minigraph -h
```

### 基本用法 | Basic Usage

```bash
# 1. 构建泛基因组图
biopytools minigraph build \
    --ref reference.fa \
    --samples sample1.fa sample2.fa sample3.fa \
    -o pangenome.gfa \
    -t 16

# 2. 调用结构变异（生成BED文件）
biopytools minigraph call \
    --graph-gfa pangenome.gfa \
    --samples sample1.fa sample2.fa \
    -o call_results/

# 3. 提取SV bubbles
biopytools minigraph bubble \
    --graph-gfa pangenome.gfa \
    -o sv_bubbles.bed

# 4. 序列到图映射
biopytools minigraph map \
    --graph-gfa pangenome.gfa \
    --queries reads.fa \
    -o mapping.gaf
```

## 子命令详解 | Subcommands Details

### 1. build - 构建泛基因组图

构建泛基因组图是 Minigraph 的核心功能，用于从多个基因组创建表示变异的图结构。

#### 基本用法

```bash
biopytools minigraph build \
    --ref reference.fa \
    --samples sample1.fa sample2.fa sample3.fa \
    -o pangenome.gfa
```

#### 高级用法

```bash
# 自定义参数
biopytools minigraph build \
    --ref reference.fa \
    --samples sample1.fa sample2.fa \
    -o pangenome.gfa \
    --preset ggs \
    --min-identity 0.95 \
    --min-aln-len 50000 \
    --max-gap 500000 \
    -t 24 \
    --batch-size 1000

# 追加模式（添加新样本到现有图）
biopytools minigraph build \
    --ref existing_graph.gfa \
    --samples new_sample.fa \
    --append-mode \
    -o updated_graph.gfa
```

#### 预设选项说明

| 预设 | 描述 | 推荐场景 |
|------|------|----------|
| **g** | 基础图构建 | 简单快速构建 |
| **gs** | 包含序列相似度 | 需要质量控制 |
| **ggs** | 包含序列相似度和基础比对 | **推荐默认**，最准确 |

### 2. call - 调用结构变异

为每个样本生成基于泛基因组图的SV调用结果（BED格式）。

#### 基本用法

```bash
biopytools minigraph call \
    --graph-gfa pangenome.gfa \
    --samples sample1.fa sample2.fa \
    -o bed_files/
```

#### 批量处理

```bash
# 处理多个样本
biopytools minigraph call \
    --graph-gfa pangenome.gfa \
    --samples *.fa \
    -o all_samples/ \
    -t 16
```

#### 输出格式

每个样本生成一个BED文件，包含：
- 染色体和位置信息
- 通过bubble的路径
- 基因型信息
- 比对质量分数

### 3. bubble - 提取SV bubbles

从泛基因组图中提取所有结构变异bubble信息。

#### 基本用法

```bash
biopytools minigraph bubble \
    --graph-gfa pangenome.gfa \
    -o sv_bubbles.bed
```

#### BED文件格式说明

输出的BED文件包含以下列：

| 列 | 描述 |
|----|------|
| 1-3 | 染色体、起始、结束位置 |
| 4 | bubble中的GFA片段数量 |
| 5 | 所有可能的路径数量 |
| 6 | 是否包含倒位 (1/0) |
| 7 | 最短路径长度 |
| 8 | 最长路径长度 |
| 9-11 | 保留字段（忽略） |
| 12 | bubble中的片段列表 |
| 13 | 最短路径序列 |
| 14 | 最长路径序列 |

#### 示例

```bash
# 提取bubbles并统计
biopytools minigraph bubble \
    --graph-gfa pangenome.gfa \
    -o bubbles.bed

# 统计bubble数量
wc -l bubbles.bed

# 提取特定类型的SV
awk '$6==1' bubbles.bed > inversions.bed  # 倒位
awk '$7>1000' bubbles.bed > large_sv.bed  # 大片段SV
```

### 4. map - 序列到图映射

将序列（reads、contigs等）映射到泛基因组图。

#### 基本用法

```bash
biopytools minigraph map \
    --graph-gfa pangenome.gfa \
    --queries query.fa \
    -o mapping.gaf
```

#### 不同测序类型

```bash
# 短读长 (Short reads)
biopytools minigraph map \
    --graph-gfa pangenome.gfa \
    --queries illumina_reads.fa \
    --preset sr \
    -o sr_mapping.gaf

# 长读长 (Long reads)
biopytools minigraph map \
    --graph-gfa pangenome.gfa \
    --queries pacbio_reads.fa \
    --preset map-pb \
    -o pb_mapping.gaf

# 基因组组装 (Assembly)
biopytools minigraph map \
    --graph-gfa pangenome.gfa \
    --queries assembly.fa \
    --preset asm \
    --max-intron-len 100000 \
    -o asm_mapping.gaf
```

#### 映射预设选项

| 预设 | 描述 | 适用场景 |
|------|------|----------|
| **sr** | 短读长 | Illumina reads |
| **lr** | 长读长 | PacBio/Nanopore reads |
| **map-pb** | PacBio映射 | PacBio CLR reads |
| **map-ont** | ONT映射 | Oxford Nanopore reads |
| **asm** | 基因组组装 | Contig/Scaffold 映射 |

## 参数说明 | Parameters

### build 子命令参数

#### 必需参数 | Required Parameters

| 参数 | 描述 | 示例 |
|------|------|------|
| `--ref` | 参考基因组FASTA文件 | `--ref reference.fa` |
| `--samples` | 样本基因组FASTA文件列表 | `--samples s1.fa s2.fa` |

#### 输出配置 | Output Configuration

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-o, --output-gfa` | `./pangenome.gfa` | 输出GFA文件路径 |

#### 构建参数 | Build Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--preset` | `ggs` | 图构建预设 (g/gs/ggs) |
| `--min-identity` | `0.9` | 最小序列相似度 |
| `--min-aln-len` | `100000` | 最小比对长度 |
| `--max-gap` | `1000000` | 最大gap大小 |

#### 性能参数 | Performance Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-t, --threads` | `16` | 线程数 |
| `--batch-size` | `None` | 批处理大小(MB) |

#### 外部工具路径 | External Tool Paths

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--minigraph-path` | `minigraph` | minigraph工具路径 |
| `--gfatools-path` | `gfatools` | gfatools工具路径 |

#### 处理选项 | Processing Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--keep-intermediate` | `False` | 保留中间文件 |
| `--append-mode` | `False` | 追加模式 |

### call 子命令参数

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--graph-gfa` | *必需* | 泛基因组图GFA文件 |
| `--samples` | *必需* | 样本基因组FASTA文件列表 |
| `-o, --output-dir` | `./minigraph_call` | 输出目录 |
| `--preset` | `asm` | SV调用预设 |
| `-t, --threads` | `16` | 线程数 |
| `--minigraph-path` | `minigraph` | minigraph工具路径 |

### bubble 子命令参数

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--graph-gfa` | *必需* | 泛基因组图GFA文件 |
| `-o, --output-bed` | `./sv_bubbles.bed` | 输出BED文件路径 |
| `--gfatools-path` | `gfatools` | gfatools工具路径 |

### map 子命令参数

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--graph-gfa` | *必需* | 泛基因组图GFA文件 |
| `--queries` | *必需* | 查询序列FASTA文件列表 |
| `-o, --output-gaf` | `./mapping.gaf` | 输出GAF文件路径 |
| `--preset` | `lr` | 映射预设 |
| `--max-intron-len` | `None` | 最大内含子长度 |
| `-t, --threads` | `16` | 线程数 |
| `--batch-size` | `None` | 批处理大小(MB) |
| `--minigraph-path` | `minigraph` | minigraph工具路径 |

## 使用示例 | Usage Examples

### 示例1：完整流程 - 从基因组到SV

```bash
# Step 1: 构建泛基因组图
biopytools minigraph build \
    --ref GRCh38.fa \
    --samples sample1.fa sample2.fa sample3.fa \
    -o human_pangenome.gfa \
    -t 24

# Step 2: 提取SV bubbles
biopytools minigraph bubble \
    --graph-gfa human_pangenome.gfa \
    -o sv_bubbles.bed

# Step 3: 为样本调用SV
biopytools minigraph call \
    --graph-gfa human_pangenome.gfa \
    --samples sample1.fa sample2.fa sample3.fa \
    -o sv_calls/
```

### 示例2：增量式添加样本

```bash
# 初始图
biopytools minigraph build \
    --ref ref.fa \
    --samples sample1.fa sample2.fa \
    -o graph_v1.gfa

# 添加新样本
biopytools minigraph build \
    --ref graph_v1.gfa \
    --samples sample3.fa sample4.fa \
    --append-mode \
    -o graph_v2.gfa
```

### 示例3：质量控制参数调整

```bash
# 严格质量控制（高相似度、长比对）
biopytools minigraph build \
    --ref ref.fa \
    --samples sample1.fa sample2.fa \
    --preset ggs \
    --min-identity 0.95 \
    --min-aln-len 200000 \
    -o strict_graph.gfa

# 宽松参数（发现更多变异）
biopytools minigraph build \
    --ref ref.fa \
    --samples sample1.fa sample2.fa \
    --preset g \
    --min-identity 0.85 \
    --min-aln-len 50000 \
    -o sensitive_graph.gfa
```

### 示例4：结合Swave使用

```bash
# 使用Minigraph构建图
biopytools minigraph build \
    --ref reference.fa \
    --samples sample1.fa sample2.fa sample3.fa \
    -o pangenome.gfa

# 生成BED文件（用于Swave）
biopytools minigraph call \
    --graph-gfa pangenome.gfa \
    --samples sample1.fa sample2.fa sample3.fa \
    -o bed_files/

# 准备assemblies.tsv
cat > assemblies.tsv << EOF
NAME\thap1
sample1\t$(pwd)/sample1.fa
sample2\t$(pwd)/sample2.fa
sample3\t$(pwd)/sample3.fa
EOF

# 使用Swave检测SV
biopytools swave call \
    -i assemblies.tsv \
    -r reference.fa \
    -g pangenome.gfa \
    -s minigraph \
    -o swave_results
```

### 示例5：大规模数据批处理

```bash
# 使用for循环处理多个样本
ref="reference.fa"
samples=($(ls samples/*.fa))

# 构建图（分批处理）
biopytools minigraph build \
    --ref $ref \
    --samples "${samples[@]}" \
    -o pangenome.gfa \
    -t 32

# 为所有样本调用SV
biopytools minigraph call \
    --graph-gfa pangenome.gfa \
    --samples "${samples[@]}" \
    -o all_calls \
    -t 32
```

## 输出结果 | Output Results

### GFA文件结构

GFA (Graphical Fragment Assembly) 格式包含三种行类型：

```
S       segment1        ACCTT...                # 序列片段
L       segment1+       segment2-       0M      # 连接
P       path1   segment1+,segment2-              # 路径
```

### call 输出 (BED格式)

每个样本生成一个BED文件，包含：
- 基因组坐标
- 通过bubble的路径
- 覆盖度和质量信息
- 样本的基因型

### bubble 输出 (BED格式)

包含所有检测到的结构变异bubble，可用于：
- SV频率分析
- 特定SV提取
- 图可视化

### map 输出 (GAF格式)

Graph Alignment Format，包含：
- 查询序列ID
- 图上的路径
- 比对质量
- CIGAR字符串

## 系统要求 | System Requirements

### 依赖软件 | Dependencies

- **minigraph** (版本 0.5 或更新)
  - GitHub: https://github.com/lh3/minigraph
  - 包含在 swave_v.1.2 conda 环境中

- **gfatools** (可选)
  - 用于 bubble 提取
  - GitHub: https://github.com/lh3/gfatools

### 编译安装 | Compilation

```bash
# 如果需要手动安装minigraph
git clone https://github.com/lh3/minigraph
cd minigraph
make
sudo cp minigraph /usr/local/bin/
```

### Conda 安装

```bash
# 使用已有环境
conda activate swave_v.1.2

# 或创建新环境
conda create -n minigraph -c bioconda minigraph
conda activate minigraph
```

## 注意事项 | Important Notes

1. **内存需求**: 构建大型泛基因组图需要大量内存（建议32GB+）
2. **样本质量**: 输入样本应高质量基因组组装，避免contigs
3. **参数选择**: 根据数据类型选择合适的预设参数
4. **增量构建**: 大样本建议分批构建图，避免一次性处理
5. **输出格式**: GFA文件可用Bandage等工具可视化

## 性能优化 | Performance Optimization

### 线程配置

```bash
# 根据可用CPU核心数调整
biopytools minigraph build \
    --ref ref.fa \
    --samples s1.fa s2.fa \
    -t $(nproc)  # 使用所有可用核心
```

### 内存管理

```bash
# 使用批处理参数控制内存
biopytools minigraph build \
    --ref ref.fa \
    --samples s1.fa s2.fa \
    --batch-size 500  # 500MB批处理
```

### 分批处理

```bash
# 处理大量样本时分批构建
for batch in batch1 batch2 batch3; do
    biopytools minigraph build \
        --ref ref.gfa \
        --samples ${batch}/*.fa \
        --append-mode \
        -o ${batch}_graph.gfa \
        -t 24
done
```

## 故障排除 | Troubleshooting

### 常见问题

**Q: "内存不足" 错误**

```bash
# 解决方案：
# 1. 减少线程数
biopytools minigraph build ... -t 8

# 2. 使用批处理
biopytools minigraph build ... --batch-size 500

# 3. 分批处理样本
```

**Q: "比对率低"**

```bash
# 调整参数提高灵敏度
biopytools minigraph build \
    ... \
    --min-identity 0.85 \
    --min-aln-len 50000
```

**Q: "图构建失败"**

```bash
# 检查输入文件
# 1. 验证FASTA格式
seqkit stats sample.fa

# 2. 检查文件大小
ls -lh sample.fa

# 3. 使用简单预设测试
biopytools minigraph build \
    ... \
    --preset g
```

## 结果解读 | Result Interpretation

### GFA文件解读

- **S行**: 序列片段，代表基因组中的连续区域
- **L行**: 连接，表示片段之间的关系
- **P行**: 路径，表示特定样本在图中的路径

### BED文件解读

- **位置**: 参考基因组上的坐标
- **Bubble**: 包含变异的区域
- **路径**: 样本选择的路径

### GAF文件解读

- 比对位置、路径、CIGAR等信息
- 可用于下游分析（如变异检测）

## 相关工具 | Related Tools

- **Bandage**: GFA图可视化
- **VG**: 另一个泛基因组图工具包
- **Gfatools**: GFA文件处理工具
- **Swave**: 基于图的SV检测（本模块）

## 参考文献 | References

Li, H. (2021). *Minigraph: a fast and efficient tool for pangenome graph construction and mapping*. Bioinformatics, 37(22), 4116-4122.

## 许可证 | License

本项目采用MIT许可证 - 详见 [LICENSE](LICENSE) 文件

Minigraph本身为MIT许可证，详见 https://github.com/lh3/minigraph
