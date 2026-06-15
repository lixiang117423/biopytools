# 基因组ID重命名工具|Genome ID Renamer

## 功能描述|Function Description

按顺序重命名FASTA文件中所有序列为Chr01, Chr02, Chr03...格式，并可选择性提取前N条作为染色体文件。适用于NCBI下载的各种基因组文件。

Rename all sequences in FASTA file in order as Chr01, Chr02, Chr03..., optionally extract first N sequences as chromosomes. Suitable for various genome files from NCBI.

## 主要特性|Main Features

- **简单直接|Simple and direct**: 按顺序重命名所有序列，无需复杂类型检测
- **染色体提取|Chromosome extraction**: 自动生成只包含染色体的文件
- **ID映射保存|ID mapping saving**: 保存原始ID和新ID的对应关系
- **灵活的格式选项|Flexible format options**: 支持自定义前缀和零填充

## 使用方法|Usage

### 基本用法|Basic Usage

```bash
# 重命名所有序列|Rename all sequences
biopytools rename-genome-id -i input.fa -o output.fa
```

### 指定染色体数量|Specify Chromosome Count

```bash
# 重命名并提取前20条作为染色体|Rename and extract first 20 as chromosomes
biopytools rename-genome-id -i input.fa -o output.fa -n 20
```

### 完整参数说明|Full Parameters

| 参数|短参数|描述|默认值|
|------|------|------|------|
| `--input`|`-i`|输入FASTA文件|必需|
| `--output`|`-o`|输出FASTA文件|必需|
| `--prefix`|`-p`|序列前缀|Chr|
| `--no-zero-padding`||不使用零填充(Chr1而非Chr01)|False|
| `--padding-width`|`-w`|填充宽度|2|
| `--chr-count`|`-n`|染色体数量|0(不提取)|
| `--no-mapping`||不保存ID映射文件|False|
| `--mapping-file`||指定映射文件路径|自动生成|
| `--log-level`||日志级别|INFO|

## 输出文件|Output Files

### 1. 完整重命名文件|Renamed File (所有序列)
包含所有序列，全部按顺序重命名：
```
output.fa (或用户指定名称)
```

### 2. 染色体文件|Chromosome-Only File (仅当使用-n参数时)
仅包含前N条序列（染色体）：
```
output_chr_only.fa
```

### 3. ID映射文件|ID Mapping File
原始ID和新ID的对应关系：
```
output_id_mapping.txt
```

格式：
```
# Original_ID	New_ID	Description
CM081539.1	Chr01	CM081539.1 ... chromosome 1 ...
CM081540.1	Chr02	CM081540.1 ... chromosome 2 ...
```

## 使用示例|Examples

### 示例1: 花生基因组|Example 1: Peanut Genome (20条染色体)

```bash
biopytools rename-genome-id \
    -i GCA_039854485.1_S83_genomic.fna \
    -o peanut_renamed.fna \
    -n 20
```

**结果|Results:**
- `peanut_renamed.fna`: 全部20条序列 → Chr01-Chr20
- `peanut_renamed_chr_only.fna`: 仅20条染色体
- `peanut_renamed_id_mapping.txt`: ID映射文件

### 示例2: 小花糖芥基因组|Example 2: Erysimum Genome (8条染色体+101个scaffold)

```bash
biopytools rename-genome-id \
    -i GCA_040802135.1_genomic.fna \
    -o ech_renamed.fna \
    -n 8
```

**结果|Results:**
- `ech_renamed.fna`: 全部109条序列 → Chr01-Chr109
- `ech_renamed_chr_only.fna`: 前8条染色体 → Chr01-Chr08
- `ech_renamed_id_mapping.txt`: ID映射文件

### 示例3: 自定义格式|Example 3: Custom Format

```bash
# 不使用零填充
biopytools rename-genome-id \
    -i input.fa -o output.fa \
    --no-zero-padding

# 三位数字填充
biopytools rename-genome-id \
    -i input.fa -o output.fa \
    -w 3

# 自定义前缀
biopytools rename-genome-id \
    -i input.fa -o output.fa \
    -p chromosome
```

## 应用场景|Use Cases

### 1. 共线性分析|Collinearity Analysis
简化后的序列ID更易读：
```
Before: CM081539.1 vs LNXXXXXX.1
After:  Chr01 vs Chr01
```

### 2. 基因组浏览器|Genome Browser
简短的ID显示更清晰

### 3. 批量更新注释|Batch Update Annotations
使用ID映射文件批量更新GFF/GTF：
```bash
# 使用映射文件更新GFF
awk 'NR==FNR{map[$1]=$2; next} $3 in map{$3=map[$3]}1' id_mapping.txt annotation.gff > annotation_renamed.gff
```

## 注意事项|Notes

1. **文件大小|File Size**: 此工具会复制整个FASTA文件，大文件可能需要较长时间和磁盘空间
2. **序列顺序|Sequence Order**: 工具按文件中的顺序重命名，不改变序列顺序
3. **染色体数量|Chromosome Count**: 使用`-n`参数时确保指定的数量正确
4. **备份原文件|Backup**: 建议保留原始FASTA文件作为备份

## 设计思路|Design Philosophy

**为什么按顺序重命名？|Why rename in order?**

1. **通用性|Universality**: 不需要复杂的序列类型检测，适用于所有基因组文件
2. **可预测性|Predictability**: 输出结果完全可预测，易于调试
3. **简单性|Simplicity**: 逻辑简单，不容易出错

**为什么支持提取染色体？|Why support chromosome extraction?**

1. **常见需求|Common Need**: 很多分析只需要染色体序列
2. **便捷性|Convenience**: 一步完成重命名和提取，无需额外命令
3. **灵活性|Flexibility**: 可以选择提取任意数量的序列

## 故障排除|Troubleshooting

### 问题1: 内存不足
**解决方案|Solution**: 对于超大基因组文件，考虑使用seqkit等工具先提取需要的序列

### 问题2: 染色体数量不对
**解决方案|Solution**: 先用`grep -c "^>" input.fa`检查序列总数，确认染色体数量

### 问题3: 需要特定序列而不是前N条
**解决方案|Solution**:
1. 使用本工具重命名所有序列
2. 使用seqkit根据ID映射文件提取特定序列

## 版本历史|Version History

- **v1.0.0** (2026-01-21): 初始版本，支持按顺序重命名和染色体提取|Initial release with sequential renaming and chromosome extraction
