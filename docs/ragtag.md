# RagTag 基因组 Scaffolding 分析模块

**专业的参考基因组scaffolding工具 | Professional Reference-Based Genome Scaffolding Tool**

## 功能概述 | Overview

RagTag 基因组scaffolding分析模块是一个基于参考基因组的序列scaffolding工具，提供完整的序列组装、scaffolding、序列ID重命名和分类输出功能。支持自动化的处理流程、灵活的参数配置和详细的输出结果，适用于各种基因组组装研究。

## 主要特性 | Key Features

- **参考基因组scaffolding**: 基于近缘物种参考基因组进行序列排序和定向
- **自动序列ID重命名**: 移除_RagTag等后缀，自动重命名序列ID
- **智能分类输出**: 自动区分scaffolded和unscaffolded序列并分别输出
- **多种比对器支持**: 支持minimap2、unimap、nucmer三种比对算法
- **灵活参数配置**: 可配置的scaffolding参数和gap推断选项
- **详细日志记录**: 完整的处理过程日志和错误追踪
- **高效处理**: 优化的处理流程，支持多线程加速

## 快速开始 | Quick Start

### 基本用法 | Basic Usage

```bash
# 基本scaffolding分析
biopytools ragtag \
    -r reference.fa \
    -q query.fa \
    -s Sample1 \
    -o output_dir

# 使用序列ID前缀
biopytools ragtag \
    -r reference.fa \
    -q query.fa \
    -s Sample1 \
    -p hifiasm \
    -o output_dir

# 使用gap推断功能
biopytools ragtag \
    -r reference.fa \
    -q query.fa \
    -s Sample1 \
    -R \
    -o output_dir
```

### 高级用法 | Advanced Usage

```bash
# 使用自定义比对器和线程数
biopytools ragtag \
    -r reference.fa \
    -q query.fa \
    -s Sample1 \
    --aligner nucmer \
    -t 24 \
    -o output_dir

# 合并未定位contigs并推断gap
biopytools ragtag \
    -r reference.fa \
    -q query.fa \
    -s Sample1 \
    -C \
    -R \
    -o output_dir
```

## 参数说明 | Parameters

### 必需参数 | Required Parameters

| 参数 | 描述 | 示例 |
|------|------|------|
| `-r, --reference` | 参考基因组FASTA文件路径 | `-r reference.fa` |
| `-q, --query` | 查询基因组FASTA文件路径 | `-q query.fa` |
| `-s, --sample-name` | 样品名称，用于输出文件命名 | `-s Sample1` |

### 处理配置 | Processing Configuration

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-t, --threads` | `12` | 线程数 |
| `-o, --output-dir` | `./ragtag_output` | 输出目录路径 |
| `-p, --prefix` | `None` | 序列ID前缀（添加到所有输出序列ID前面） |
| `--aligner` | `minimap2` | 比对器 (minimap2/unimap/nucmer) |

### RagTag选项 | RagTag Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-C, --concatenate-unplaced` | `False` | 将未定位的contigs合并为chr0 |
| `-R, --infer-gaps` | `False` | 推断gap大小 |

## 输入文件格式 | Input File Formats

### 参考基因组文件 | Reference Genome File

标准FASTA格式的参考基因组序列：

```fasta
>chromosome1
ATCGATCGATCGATCGATCGATCGATCG...
>chromosome2
GCTAGCTAGCTAGCTAGCTAGCTAGCTA...
```

### 查询基因组文件 | Query Genome File

标准FASTA格式的查询基因组序列（待scaffold的序列）：

```fasta
>contig1
ATCGATCGATCGATCGATCGATCGATCG...
>contig2
GCTAGCTAGCTAGCTAGCTAGCTAGCTA...
```

## 使用示例 | Usage Examples

### 示例1：基本scaffolding分析 | Example 1: Basic Scaffolding Analysis

```bash
# 对植物基因组进行scaffolding
biopytools ragtag \
    -r reference_genome.fa \
    -q assembled_contigs.fa \
    -s PlantSample1 \
    -o scaffolding_results
```

### 示例2：使用序列ID前缀 | Example 2: With Sequence ID Prefix

```bash
# 使用前缀标识序列来源
biopytools ragtag \
    -r reference.fa \
    -q query.fa \
    -s Sample2 \
    -p hifiasm \
    -o results_with_prefix

# 输出序列ID示例：
# >hifiasm_Chr1
# >hifiasm_Chr2
# >hifiasm_scaffold_1
# >hifiasm_scaffold_2
```

### 示例3：使用gap推断 | Example 3: With Gap Inference

```bash
# 使用gap推断功能获得更准确的gap大小
biopytools ragtag \
    -r reference.fa \
    -q query.fa \
    -s Sample3 \
    -R \
    -t 16 \
    -o results_with_gaps
```

### 示例4：使用不同比对器 | Example 4: Using Different Aligner

```bash
# 使用nucmer比对器
biopytools ragtag \
    -r reference.fa \
    -q query.fa \
    -s Sample3 \
    --aligner nucmer \
    -t 24 \
    -o nucmer_results
```

### 示例5：合并未定位contigs | Example 5: Concatenate Unplaced Contigs

```bash
# 将未定位的contigs合并为chr0
biopytools ragtag \
    -r reference.fa \
    -q query.fa \
    -s Sample5 \
    -C \
    -R \
    -o complete_results
```

## 输出结果 | Output Results

### 输出目录结构 | Output Directory Structure

```
ragtag_output/
├── Sample1_RagTag_scaffolded.fa           # Scaffolded序列（已重命名）
├── Sample1_RagTag_unscaffolded.fa         # Unscaffolded序列（已重命名为scaffold_1, scaffold_2...）
├── Sample1_ragtag.scaffold.fasta          # RagTag原始输出
├── Sample1_ragtag.scaffold.agp            # AGP格式文件
├── Sample1_ragtag.scaffold.asm.paf        # 比对文件
├── Sample1_ragtag.scaffold.asm.paf.log    # 比对日志
├── Sample1_ragtag.scaffold.confidence.txt # 置信度文件
├── Sample1_ragtag.scaffold.err            # 错误日志
├── Sample1_ragtag.scaffold.stats          # 统计信息
└── Sample1_ragtag.log                     # 运行日志
```

### 关键输出文件说明 | Key Output Files Description

- **{sample}_RagTag_scaffolded.fa**: 已成功scaffold的序列，序列ID已重命名
  - `Chr12_RagTag` → `Chr12`
  - 染色体级别的序列

- **{sample}_RagTag_unscaffolded.fa**: 未成功scaffold的序列
  - `ptg000019l` → `scaffold_1`
  - `ptg000020l` → `scaffold_2`
  - 按顺序编号的未定位序列

- **{sample}_ragtag.scaffold.fasta**: RagTag原始输出文件
- **{sample}_ragtag.scaffold.agp**: AGP格式的scaffolding结果
- **{sample}_ragtag.scaffold.stats**: Scaffolding统计信息
- **{sample}_ragtag.log**: 完整的运行日志

### 序列ID重命名规则 | Sequence ID Renaming Rules

1. **移除_RagTag后缀**
   - `Chr12_RagTag` → `Chr12`
   - `scaffold_123_RagTag` → `scaffold_123`

2. **未定位序列重命名**
   - `ptg000019l` → `scaffold_1`
   - `ptg000020l` → `scaffold_2`
   - 按顺序从1开始编号

3. **保留原始染色体名称**
   - 染色体级别的序列保持原有名称（去除后缀后）

4. **可选前缀添加（使用-p参数）**
   - 如果使用 `-p hifiasm`，所有序列ID会添加前缀：
     - `Chr12` → `hifiasm_Chr12`
     - `scaffold_1` → `hifiasm_scaffold_1`

## 系统要求 | System Requirements

### 依赖软件 | Dependencies

- **RagTag** (版本 1.1.0 或更新)
  - 下载地址: https://github.com/malonglu/RagTag
- **Python** (版本 3.7+)
- **Python包**:
  - `click` - 命令行界面

### 安装依赖软件 | Installing Dependencies

```bash
# 安装RagTag
conda install -c bioconda ragtag
# 或
pip install ragtag

# 安装Python包
pip install click
```

### 硬件建议 | Hardware Recommendations

- **CPU**: 多核处理器（推荐4核以上）
- **RAM**: 最少4GB（大基因组推荐16GB以上）
- **存储**: 预留基因组文件大小3倍的磁盘空间

## 注意事项 | Important Notes

1. **参考基因组选择**: 参考基因组应与查询基因组亲缘关系较近
2. **序列ID格式**: 确保输入FASTA文件的序列ID格式正确
3. **样品名称**: 样品名称将用于输出文件命名，建议使用有意义的名称
4. **比对器选择**: minimap2适用于大多数情况，nucmer适用于差异较大的基因组
5. **线程数配置**: 合理配置线程数可以显著提高处理速度

## 故障排除 | Troubleshooting

### 常见问题 | Common Issues

**Q: "RagTag not found" 错误**
```bash
# 检查RagTag安装
which ragtag.py

# 使用conda安装
conda install -c bioconda ragtag
```

**Q: 没有scaffolded序列输出**
```bash
# 检查参考基因组和查询基因组的亲缘关系
# 尝试降低比对严格度参数（需要在config中调整）
```

**Q: 内存不足错误**
```bash
# 减少线程数
biopytools ragtag ... -t 4

# 或分批处理大基因组
```

**Q: 序列ID没有正确重命名**
```bash
# 检查原始序列ID格式
# 确保样品名称参数正确
```

## 相关资源 | Related Resources

- [RagTag官方文档](https://github.com/malonglu/RagTag)
- [基因组scaffolding最佳实践](https://www.nature.com/articles/s41587-020-00750-5)
- [FASTA格式规范](https://en.wikipedia.org/wiki/FASTA_format)

## 许可证 | License

本项目采用MIT许可证 - 详见 [LICENSE](LICENSE) 文件

## 引用信息 | Citation

如果在学术研究中使用此工具，请引用RagTag相关文献：

```
Alonge, M., Williams, A.L. and Grimwood, J. et al.
RagTag: Reference-guided assembly and scaffolding of draft genomes.
 bioRxiv (2019) doi: 10.1101/667985
```
