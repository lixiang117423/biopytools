# Hi-C数据质量评估模块

**使用pairtools评估Hi-C mapping数据质量 | Assess Hi-C Mapping Data Quality with pairtools**

## 功能概述 | Overview

Hi-C数据质量评估模块用于评估Hi-C数据mapping结果的质量，通过pairtools计算关键质量指标并与预设阈值进行比较，生成详细的质量评估报告。

## 主要特性 | Key Features

- **自动化质量评估**: 自动运行pairtools stats并解析结果
- **多维度指标**: 评估未比对率、单端比对率、双端比对率、PCR重复率和cis/trans比例
- **灵活的阈值**: 可自定义各项指标的评估阈值
- **详细报告**: 生成包含统计数据和评估结果的详细报告
- **标准日志**: 符合BioPyTools开发规范的日志输出格式

## 快速开始 | Quick Start

### 基本用法 | Basic Usage

```bash
# 使用pairs文件进行评估（推荐）
biopytools hic-qc -i sample.pairs.gz

# 使用BAM文件进行评估（需要提供chrom.sizes）
biopytools hic-qc -i sample.bam -c assembly.chrom.sizes

# 指定输出目录
biopytools hic-qc -i sample.pairs.gz -o qc_results
```

### 高级用法 | Advanced Usage

```bash
# 自定义质量阈值
biopytools hic-qc -i sample.bam -c assembly.chrom.sizes \
    --max-unmapped-rate 15 \
    --max-single-sided-rate 8 \
    --min-mapped-rate 85 \
    --max-dup-rate 25 \
    --min-cis-trans-ratio 5

# 指定pairtools路径
biopytools hic-qc -i sample.pairs.gz \
    -p /path/to/pairtools
```

## 参数说明 | Parameters

### 必需参数 | Required Parameters

| 参数 | 描述 | 示例 |
|------|------|------|
| `-i, --input` | 输入的pairs或BAM文件路径 | `-i sample.pairs.gz` 或 `-i sample.bam` |

### 可选参数 | Optional Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-p, --pairtools-path` | `/share/org/YZWL/yzwl_lixg/miniforge3/envs/pairtools_v.1.1.3/bin/pairtools` | Pairtools可执行文件路径 |
| `-o, --output-dir` | `./pairtools_qc_output` | 输出目录 |
| `-c, --chroms-path` | 无 | Chromosome sizes文件路径（BAM输入时必需）|

### 质量阈值参数 | Quality Threshold Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--max-unmapped-rate` | `20.0` | 未比对reads比例阈值(%) |
| `--max-single-sided-rate` | `10.0` | 单端比对比例阈值(%) |
| `--min-mapped-rate` | `80.0` | 双端比对率阈值(%) |
| `--max-dup-rate` | `30.0` | PCR重复率阈值(%) |
| `--min-cis-trans-ratio` | `4.0` | Cis/Trans比例阈值 |

## 质量指标 | Quality Metrics

### 1. 未比对reads比例 | Unmapped Reads Rate

- **含义**: 两条reads都没比对上的比例
- **阈值**: < 20%
- **计算**: `total_unmapped / total_pairs * 100`
- **重要性**: 高比例说明mapping质量差或基因组不匹配

### 2. 单端比对比例 | Single-sided Mapping Rate

- **含义**: 只有一条read比对上的比例
- **阈值**: < 10%
- **计算**: `total_single_sided_mapped / total_pairs * 100`
- **重要性**: 高比例说明数据质量不均或测序深度不足

### 3. 双端比对率 | Paired Mapping Rate

- **含义**: 两条reads都比对上的比例
- **阈值**: > 80%
- **计算**: `total_mapped / total_pairs * 100`
- **重要性**: 主要评估指标，反映整体mapping质量

### 4. PCR重复率 | PCR Duplication Rate

- **含义**: PCR重复序列占比对序列的比例
- **阈值**: < 30%
- **计算**: `total_dups / total_mapped * 100`
- **重要性**: 高比例说明PCR扩增过度，影响数据多样性

### 5. Cis/Trans比例 | Cis/Trans Ratio

- **含义**: 染色体内pairs与染色体间pairs的比值
- **阈值**: > 4
- **计算**: `cis / trans`
- **重要性**: 评估Hi-C实验质量的关键指标，高比值说明数据质量好

## 输出结果 | Output Results

### 报告格式 | Report Format

```
------------------------------------------------------------
Hi-C数据质量评估报告|Hi-C Data Quality Assessment Report
------------------------------------------------------------

统计数据|Statistics:
  总reads数|Total reads: 100,000,000
  未比对reads|Unmapped reads: 15,000,000 (15.00%)
  单端比对|Single-sided mapped: 5,000,000 (5.00%)
  双端比对|Paired mapped: 80,000,000 (80.00%)
  PCR重复|PCR duplicates: 16,000,000 (20.00%)

  染色体内pairs|Cis pairs: 60,000,000
  染色体间pairs|Trans pairs: 4,000,000
  Cis/Trans比例|Cis/Trans ratio: 15.00

质量评估|Quality Assessment:
  [通过|PASS] 未比对reads比例|Unmapped rate: 15.00% (阈值|threshold <= 20%)
  [通过|PASS] 单端比对比例|Single-sided rate: 5.00% (阈值|threshold <= 10%)
  [通过|PASS] 双端比对率|Paired mapping rate: 80.00% (阈值|threshold >= 80%)
  [通过|PASS] PCR重复率|Duplication rate: 20.00% (阈值|threshold <= 30%)
  [通过|PASS] Cis/Trans比例|Cis/Trans ratio: 15.00 (阈值|threshold >= 4)

总体评估|Overall: 通过|PASSED
------------------------------------------------------------
```

### 结果解读 | Result Interpretation

| 结果 | 含义 | 建议 |
|------|------|------|
| **全部通过** | 数据质量良好，可用于下游分析 | 继续分析 |
| **部分未通过** | 某些指标未达标 | 检查对应问题并考虑是否需要过滤数据 |
| **多数未通过** | 数据质量较差 | 建议重新测序或重新mapping |

## 使用示例 | Usage Examples

### 示例1：基本质量评估 | Example 1: Basic Quality Assessment

```bash
# 评估Hi-C数据质量
biopytools hic-qc -i SRR123456.pairs.gz
```

### 示例2：严格质量标准 | Example 2: Strict Quality Standards

```bash
# 使用更严格的质量阈值
biopytools hic-qc -i SRR123456.pairs.gz \
    --max-unmapped-rate 10 \
    --max-single-sided-rate 5 \
    --min-mapped-rate 90 \
    --max-dup-rate 20 \
    --min-cis-trans-ratio 10
```

### 示例3：批量评估多个样本 | Example 3: Batch Assessment of Multiple Samples

```bash
# 评估多个样本
for sample in *.pairs.gz; do
    biopytools hic-qc -i "$sample" -o "qc_$(basename $sample .pairs.gz)"
done
```

### 示例4：在Hi-C分析流程中使用 | Example 4: Use in Hi-C Analysis Pipeline

```bash
# 在Hi-C分析流程中进行质量控制
biopytools hic-qc -i aligned.pairs.gz -o qc_results

# 根据QC结果决定是否继续
if [ $? -eq 0 ]; then
    echo "质量评估通过，继续分析"
    # 继续下游分析...
else
    echo "质量评估未通过，请检查数据"
fi
```

## 输入文件格式 | Input File Format

### 支持的输入格式 | Supported Input Formats

本工具支持两种输入格式：

1. **Pairs/Pairsam格式**（推荐）| Pairs/Pairsam format (recommended)
   - 直接进行质量评估
   - 处理速度更快
   - 示例：`sample.pairs.gz`, `sample.pairsam.gz`

2. **BAM/SAM格式** | BAM/SAM format
   - 需要提供 `--chroms-path` 参数指定chrom.sizes文件
   - 程序会自动调用 `pairtools parse` 转换为pairs格式
   - 转换后的pairs文件会保存在输出目录中
   - 示例：`sample.bam`, `sample.sam`

### Chrom.sizes文件 | Chrom.sizes File

使用BAM/SAM文件作为输入时，必须提供chrom.sizes文件。这是一个简单的两列文本文件：

```text
chr1    248956422
chr2    242193529
chr3    198295559
...
```

**生成chrom.sizes文件**：

```bash
# 从fasta文件生成
samtools faidx assembly.fa
cut -f1,2 assembly.fa.fai > assembly.chrom.sizes

# 或使用seqkit
seqkit faidx -i assembly.fa | cut -f1,2 > assembly.chrom.sizes
```

### Pairs文件格式 | Pairs File Format

输入文件应为pairtools生成的.pairs或.pairsam格式文件，可以是压缩格式(.gz)。

```text
#columns: chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type
chr1    1000    chr1    50000   +       +       LL
chr1    2000    chr2    30000   -       +       LU
chr2    1500    chr2    80000   +       -       UU
...
```

**文件要求**:
- 标准pairs/pairsam格式
- 包含完整的header信息
- 可以是gzip压缩格式

## 常见问题 | FAQ

### Q1: 为什么要评估这些指标？

**A**: 这些指标是评估Hi-C数据质量的关键：
- **未比对率**: 反映基因组匹配程度
- **单端比对率**: 反映测序质量
- **双端比对率**: 反映整体mapping效果
- **PCR重复率**: 反映文库复杂度
- **Cis/Trans比例**: 反映Hi-C实验质量

### Q2: 如何调整质量阈值？

**A**: 根据具体研究需求和数据类型调整：
- 染色体级别组装通常要求更严格
- scaffold级别组装可以适当放宽
- 某些特殊物种可能需要特殊标准

### Q3: QC未通过怎么办？

**A**: 可以尝试：
1. 检查基因组版本是否正确
2. 调整mapping参数重新比对
3. 检查测序数据质量
4. 考虑是否需要更多数据

### Q4: Cis/Trans比例为什么很重要？

**A**: Cis/Trans比例是Hi-C实验质量的关键指标：
- 正常的Hi-C数据应该有更多的染色体内相互作用
- 低比值可能说明实验失败或数据质量问题
- 通常比值在5-20之间是正常的

### Q5: 支持哪些输入格式？

**A**: 支持pairtools生成的所有格式：
- .pairs
- .pairsam
- .gz压缩格式

## 系统要求 | System Requirements

### 依赖软件 | Dependencies

- **pairtools** (版本 1.1.3 或更新)
- **Python** (版本 3.7+)
- **Python包**:
  - `click` - 命令行界面
  - `pathlib` - 路径处理

### 安装pairtools | Installing pairtools

```bash
# 使用conda安装
conda install -c conda-forge -c bioconda pairtools

# 或从源码安装
pip install pairtools
```

## 技术细节 | Technical Details

### 评估流程 | Assessment Pipeline

1. **运行pairtools stats**: 生成统计文件
2. **解析统计结果**: 提取关键指标
3. **计算百分比和比例**: 转换为可比较的格式
4. **与阈值比较**: 评估各项指标
5. **生成报告**: 输出详细评估结果

### 性能说明 | Performance Notes

- 处理速度取决于pairs文件大小
- 大文件(>10GB)可能需要较长时间
- 建议在计算服务器上运行大样本

## 相关资源 | Related Resources

- [pairtools官方文档](https://pairtools.readthedocs.io/)
- [Hi-C数据分析最佳实践](https://www.4dnucleome.org/)
- [BioPyTools开发规范](/share/org/YZWL/yzwl_lixg/software/scripts/develop_python_guides.md)

## 许可证 | License

本项目采用MIT许可证 - 详见 [LICENSE](LICENSE) 文件

## 引用信息 | Citation

如果在学术研究中使用此工具，请引用pairtools：

```
S. V. 4D Nucleome Consortium, et al.
pairtools: ecosystem for pre-processing chromatin conformation capture data.
Bioinformatics, 2024.
```
