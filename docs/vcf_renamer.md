# 🔄 VCF样品名称重命名模块

**防止长样本名被软件截断的专业工具 | Professional Tool to Prevent Sample Name Truncation**

## 📖 功能概述 | Overview

VCF样品名称重命名模块是一个专门用于解决VCF文件中长样本名被软件截断问题的工具。通过将复杂的样本名简化为序列化的短名称（如 S1, S2, S3...），确保与各类生物信息学软件的兼容性，同时保留完整的名称映射关系以供追溯。

## ✨ 主要特性 | Key Features

- **🔄 自动化重命名**: 将复杂样本名转换为序列化短名称
- **📋 映射文件保留**: 自动生成新旧名称对照表
- **⚙️ 灵活前缀**: 支持自定义新样本名前缀
- **🔍 详细日志**: 完整的处理过程记录和调试信息
- **📦 自动索引**: 重命名后自动创建VCF索引文件
- **🛡️ 输入验证**: 完善的文件格式和依赖检查
- **🚀 高效处理**: 基于bcftools的快速处理流程
- **📊 处理报告**: 详细的统计信息和运行摘要

## 🚀 快速开始 | Quick Start

### 基本用法 | Basic Usage

```bash
# 基本重命名操作
biopytools vcf-renamer \
    -i variation.filtered.snp.vcf.gz \
    -o variation.renamed.vcf.gz

# 使用自定义前缀
biopytools vcf-renamer \
    -i input.vcf.gz \
    -o output.vcf.gz \
    -p Sample
```

### 高级用法 | Advanced Usage

```bash
# 详细输出模式
biopytools vcf-renamer \
    -i input.vcf.gz \
    -o output.vcf.gz \
    -v

# 保存日志文件
biopytools vcf-renamer \
    -i input.vcf.gz \
    -o output.vcf.gz \
    --log-file rename.log \
    -v
```

## 📋 参数说明 | Parameters

### 必需参数 | Required Parameters

| 参数 | 描述 | 示例 |
|------|------|------|
| `-i, --input` | 输入VCF文件路径 | `-i input.vcf.gz` |
| `-o, --output` | 输出VCF文件路径 | `-o output.vcf.gz` |

### 可选参数 | Optional Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-p, --prefix` | `S` | 新样品名的前缀 |
| `--log-file` | `None` | 日志文件路径 |
| `-v, --verbose` | `0` | 详细输出级别 |
| `--quiet` | `False` | 静默模式 |

### 输出文件 | Output Files

| 文件 | 描述 |
|------|------|
| `<output>.vcf.gz` | 重命名后的VCF文件 |
| `<output>.vcf.gz.tbi` | VCF索引文件 |
| `vcf_name_mapping.txt` | 新旧样品名称映射表 |

## 💡 使用示例 | Usage Examples

### 示例1：植物群体样本重命名 | Example 1: Plant Population Sample Renaming

```bash
# 将复杂的样本名简化为S1, S2, S3...
biopytools vcf-renamer \
    -i Arabidopsis_population.vcf.gz \
    -o Arabidopsis_renamed.vcf.gz
```

**映射文件示例** (`vcf_name_mapping.txt`):
```
Col_0_leaf_rep1	S1
Col_0_leaf_rep2	S2
Ler_leaf_rep1	S3
Ler_leaf_rep2	S4
```

### 示例2：使用描述性前缀 | Example 2: Using Descriptive Prefix

```bash
# 使用"Sample"作为前缀
biopytools vcf-renamer \
    -i human_wes_variants.vcf.gz \
    -o human_renamed.vcf.gz \
    -p Sample
```

**结果样本名**: `Sample1`, `Sample2`, `Sample3`, ...

### 示例3：微生物样本批处理 | Example 3: Microbial Sample Batch Processing

```bash
# 批量处理多个VCF文件
for file in strain_*.vcf.gz; do
    output="renamed_$(basename $file)"
    biopytools vcf-renamer -i "$file" -o "$output" -p Strain
done
```

### 示例4：详细日志记录 | Example 4: Detailed Logging

```bash
# 启用详细日志并保存到文件
biopytools vcf-renamer \
    -i large_dataset.vcf.gz \
    -o large_renamed.vcf.gz \
    --log-file renaming_process.log \
    -v
```

### 示例5：质量控制和验证 | Example 5: Quality Control and Validation

```bash
# 重命名后验证样本数量
biopytools vcf-renamer -i input.vcf.gz -o output.vcf.gz -v

# 检查映射文件
cat vcf_name_mapping.txt

# 验证新VCF的样本名
bcftools query -l output.vcf.gz
```

## 📁 输入文件格式 | Input File Formats

### 输入VCF文件 | Input VCF File

必须为压缩的VCF格式（`.vcf.gz` 或 `.bcf`）：

```vcf
##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	VeryLongSampleName_A	VeryLongSampleName_B	VeryLongSampleName_C
chr1	100	.	A	G	30	PASS	DP=100	GT:DP	0/1:50	1/1:48	0/0:52
```

**文件要求**:
- 必须是压缩格式（`.vcf.gz` 或 `.bcf`）
- 必须有对应的索引文件（`.tbi` 或 `.csi`）
- 标准VCF格式（VCF v4.0+）

## 📊 输出结果 | Output Results

### 重命名后的VCF | Renamed VCF

```vcf
##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3
chr1	100	.	A	G	30	PASS	DP=100	GT:DP	0/1:50	1/1:48	0/0:52
```

### 映射文件格式 | Mapping File Format

文本格式，包含原始样本名和新样本名的对应关系：

```
VeryLongSampleName_A	S1
VeryLongSampleName_B	S2
VeryLongSampleName_C	S3
```

**字段说明**:
- 第1列：原始样本名
- 第2列：新样本名（前缀 + 序号）

## 🔧 系统要求 | System Requirements

### 依赖软件 | Dependencies

- **bcftools** (版本 1.10+)
  - 用于VCF文件操作
  - 安装命令：`conda install bcftools` 或 `sudo apt-get install bcftools`

- **Python** (版本 3.7+)
- **Python包**:
  - `click` - 命令行界面
  - `pathlib` - 路径处理
  - `subprocess` - 系统调用
  - `logging` - 日志记录

### 安装依赖软件 | Installing Dependencies

```bash
# 安装bcftools
# Conda
conda install bcftools

# Ubuntu/Debian
sudo apt-get install bcftools

# CentOS/RHEL
sudo yum install bcftools

# 验证安装
bcftools --version
```

## 🔄 工作流程 | Workflow

```
输入VCF (.vcf.gz)
    ↓
步骤1: 提取原始样品名称
    - 使用 bcftools query -l
    - 保存到临时文件
    ↓
步骤2: 生成新旧ID映射表
    - 格式: 旧名称 新名称
    - 例如: SampleA S1
    ↓
步骤3: 重命名VCF样品
    - 使用 bcftools reheader
    - 生成新的VCF文件
    ↓
步骤4: 建立VCF索引
    - 使用 bcftools index -t
    - 生成 .tbi 索引文件
    ↓
输出文件
    - 重命名后的VCF
    - VCF索引
    - 名称映射表
```

## ⚠️ 注意事项 | Important Notes

1. **样本数量**: 支持任意数量的样本，从1个到数千个
2. **文件格式**: 输入必须是压缩的VCF格式（`.vcf.gz` 或 `.bcf`）
3. **索引要求**: 输入VCF必须有对应的索引文件
4. **前缀选择**: 建议使用简短的字母前缀（如S、Sample、Strain等）
5. **映射保存**: 务必保存映射文件以备后续分析追溯
6. **文件覆盖**: 输出文件已存在时会报错，请先删除旧文件或使用新文件名

## 🐛 故障排除 | Troubleshooting

### 常见问题 | Common Issues

**Q: "bcftools: command not found" 错误**
```bash
# 安装bcftools
conda install bcftools
# 或
sudo apt-get install bcftools

# 验证安装
bcftools --version
```

**Q: 输入文件格式错误**
```bash
# 检查文件格式
file input.vcf.gz

# 如果不是压缩格式，先压缩
bgzip input.vcf -o input.vcf.gz
bcftools index input.vcf.gz
```

**Q: 缺少索引文件**
```bash
# 创建索引
bcftools index -t input.vcf.gz
# 或
bcftools index input.vcf.gz  # 自动选择索引格式
```

**Q: 重命名后样本名验证**
```bash
# 查看重命名后的样本名
bcftools query -l output.vcf.gz

# 查看原始样本名
bcftools query -l input.vcf.gz
```

## 💻 高级用法 | Advanced Usage

### 与其他工具集成 | Integration with Other Tools

```bash
# 1. 重命名VCF文件
biopytools vcf-renamer -i raw.vcf.gz -o renamed.vcf.gz

# 2. 使用重命名后的VCF进行后续分析
bcftools stats renamed.vcf.gz > stats.txt

# 3. 根据映射文件追溯原始样本名
grep "S5" vcf_name_mapping.txt
```

### 批量处理脚本 | Batch Processing Script

```bash
#!/bin/bash
# batch_rename.sh

INPUT_DIR="./raw_vcfs"
OUTPUT_DIR="./renamed_vcfs"
PREFIX="Sample"

mkdir -p "$OUTPUT_DIR"

for vcf in "$INPUT_DIR"/*.vcf.gz; do
    basename=$(basename "$vcf")
    output="$OUTPUT_DIR/renamed_$basename"

    echo "Processing: $basename"
    biopytools vcf-renamer -i "$vcf" -o "$output" -p "$PREFIX" -v
done
```

## 📚 应用场景 | Application Scenarios

1. **群体遗传学分析**: 重命名长样本名以兼容GATK、VCFtools等软件
2. **植物育种研究**: 统一不同批次的样本命名规范
3. **微生物基因组**: 规范化菌株名称
4. **临床检测**: 简化患者样本ID以保护隐私
5. **数据共享**: 在数据共享前匿名化样本信息

## 📄 许可证 | License

本项目采用MIT许可证 - 详见 [LICENSE](LICENSE) 文件

---

## 🔬 引用信息 | Citation

如果在学术研究中使用此工具，请引用相关文献：

```
[您的项目名称] (2025).
VCF Sample Renamer: A tool for preventing sample name truncation in VCF files.
https://github.com/yourusername/biopytools
```
