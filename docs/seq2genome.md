# Seq2Genome - 序列到基因组比对工具

## 功能介绍|Overview

Seq2Genome是一个自动检测序列类型并选择合适比对工具的分析流程：

- **DNA序列**：使用Minimap2进行比对
- **蛋白质序列**：使用Miniprot进行比对
- **自动检测**：根据序列特征自动判断是DNA还是蛋白质

## 主要特性|Features

1. ✅ **自动序列类型检测** - 智能判断输入是DNA还是蛋白质序列
2. ✅ **双模式支持** - 同时支持DNA序列和蛋白质序列比对
3. ✅ **统一输出格式** - PAF、GFF3、BED格式输出
4. ✅ **向后兼容** - 完全兼容原有pep2genome功能

## 安装依赖|Dependencies

```bash
# DNA序列比对（Minimap2）
conda install -c bioconda minimap2

# 蛋白质序列比对（Miniprot）
conda install -c bioconda miniprot

# Python依赖
pip install biopython
```

## 使用方法|Usage

### 基本用法（自动检测）

```bash
# 自动检测序列类型
biopytools seq2genome --genome genome.fa --query sequences.fa -o results
```

### 指定序列类型

```bash
# DNA序列比对
biopytools seq2genome --genome genome.fa --query dna.fa --query-type dna -o results

# 蛋白质序列比对
biopytools seq2genome --genome genome.fa --query protein.fa --query-type protein -o results
```

### 完整参数

```bash
biopytools seq2genome \
  --genome genome.fa \
  --query sequences.fa \
  --query-type auto \  # auto/dna/protein
  -o results \
  -t 12 \  # 线程数
  --minimap2-path minimap2 \
  --miniprot-path miniprot
```

## 输出文件|Output Files

- `alignment.paf` - 比对结果（PAF格式）
- `alignment_statistics.txt` - 统计报告
- `alignment.gff3` - GFF3格式注释
- `alignment.bed` - BED格式注释
- `alignment.fa` - 提取的基因组序列

## 序列类型检测逻辑|Sequence Type Detection

工具通过分析FASTA文件中的序列字符组成来判断类型：

- **DNA序列特征**：90%以上字符为ATCGN
- **蛋白质序列特征**：30%以上字符为20种氨基酸字母

检测置信度会输出到日志中，如果自动检测失败，可以手动指定`--query-type`参数。

## 向后兼容|Backward Compatibility

原有pep2genome用户可以继续使用旧参数名：

```bash
# 原有命令仍然有效
biopytools seq2genome --genome genome.fa --protein protein.fa -o results
```

## 示例|Examples

### 示例1：DNA序列定位

```bash
# 找出DNA序列在基因组上的位置
biopytools seq2genome \
  --genome reference.fa \
  --query target_sequences.fa \
  -o dna_mapping_results
```

### 示例2：蛋白质基因定位

```bash
# 将蛋白质序列定位到基因组
biopytools seq2genome \
  --genome genome.fa \
  --query proteins.fa \
  --query-type protein \
  -o protein_mapping_results
```

## 技术细节|Technical Details

### Minimap2参数（DNA序列）

- Preset: `asm5` (用于~5%差异的序列)
- 输出格式: PAF
- 只输出primary alignments

### Miniprot参数（蛋白质序列）

- 标准蛋白质比对模式
- 输出格式: PAF

## 版本历史|Version History

- **v2.0.0** - 重命名为seq2genome，添加DNA序列支持和自动检测
- **v1.0.0** - 原pep2genome，仅支持蛋白质序列

## 相关模块|Related Modules

- `assembly_qc` - 基因组组装质量评估
- `bam_cov` - BAM覆盖度分析
- `primers` - 引物设计

## 联系方式|Contact

如有问题或建议，请联系 BioPyTools 开发团队。
