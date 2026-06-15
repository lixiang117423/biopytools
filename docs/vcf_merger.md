# VCF按染色体合并工具 | VCF Merge by Chromosome Tool

**专业的VCF文件按染色体合并工具 | Professional VCF File Merging Tool by Chromosome**

## 功能概述 | Overview

VCF按染色体合并工具是一个高效的VCF文件合并工具，专门用于将按染色体分割的VCF文件自动识别并合并。基于bcftools concat实现，支持并行处理、自动索引生成和详细的日志记录，适用于各种基因组变异分析流程。

## 主要特性 | Key Features

- **智能染色体识别**: 自动从文件名中提取染色体编号并分组
- **高效合并处理**: 使用bcftools concat进行快速合并
- **并行处理支持**: 多线程处理提高合并效率
- **自动索引生成**: 可选的.csi索引文件自动创建
- **灵活文件匹配**: 支持自定义文件名模式匹配
- **详细日志记录**: 完整的处理过程日志和错误追踪
- **压缩输出支持**: 输出压缩的.gz格式VCF文件
- **标准遵循**: 完全符合VCF格式规范

## 快速开始 | Quick Start

### 基本用法 | Basic Usage

```bash
# 基本VCF文件合并
biopytools vcf-merger \
    -i /path/to/vcf_files \
    -o /path/to/output

# 使用多线程加速
biopytools vcf-merger \
    -i /path/to/vcf_files \
    -o /path/to/output \
    -t 8
```

### 高级用法 | Advanced Usage

```bash
# 自定义文件模式和多线程处理
biopytools vcf-merger \
    -i /path/to/vcf_files \
    -o /path/to/output \
    --pattern "*.vcf.gz" \
    -t 16

# 不生成索引文件（节省时间）
biopytools vcf-merger \
    -i /path/to/vcf_files \
    -o /path/to/output \
    --no-index

# 详细输出模式
biopytools vcf-merger \
    -i /path/to/vcf_files \
    -o /path/to/output \
    -v
```

## 参数说明 | Parameters

### 必需参数 | Required Parameters

| 参数 | 描述 | 示例 |
|------|------|------|
| `-i, --input-dir` | 输入VCF文件目录 | `-i /data/vcf_files` |
| `-o, --output-dir` | 输出目录 | `-o /data/merged_vcf` |

### 可选参数 | Optional Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-t, --threads` | `4` | 使用的线程数 |
| `--pattern` | `*.joint.vcf.gz` | VCF文件名模式 |
| `--no-index` | `False` | 不生成索引文件 |

### 日志参数 | Logging Parameters

| 参数 | 描述 |
|------|------|
| `-v, --verbose` | 详细输出模式 (-v: INFO, -vv: DEBUG) |
| `--quiet` | 静默模式 (只输出ERROR) |

## 输入文件格式 | Input File Format

### VCF文件命名规则 | VCF File Naming Rules

工具通过文件名自动识别染色体编号，支持以下命名格式：

```
# 支持的命名格式 | Supported naming formats:
Chr19_1-20000000.joint.vcf.gz    -> 识别为 Chr19
chr19_1-20000000.joint.vcf.gz    -> 识别为 chr19
19_1-20000000.joint.vcf.gz       -> 识别为 19
Chr19_1-20000000.vcf.gz          -> 识别为 Chr19
```

### 输入目录结构示例 | Input Directory Structure Example

```
input_vcf/
├── Chr1_1-25000000.joint.vcf.gz
├── Chr1_25000001-50000000.joint.vcf.gz
├── Chr2_1-30000000.joint.vcf.gz
├── Chr2_30000001-60000000.joint.vcf.gz
├── Chr3_1-20000000.joint.vcf.gz
└── ...
```

## 输出结果 | Output Results

### 输出目录结构 | Output Directory Structure

```
output_dir/
├── Chr1.joint.merged.vcf.gz       # 染色体1合并后的VCF文件
├── Chr1.joint.merged.vcf.gz.csi   # 染色体1索引文件
├── Chr2.joint.merged.vcf.gz       # 染色体2合并后的VCF文件
├── Chr2.joint.merged.vcf.gz.csi   # 染色体2索引文件
├── Chr3.joint.merged.vcf.gz       # 染色体3合并后的VCF文件
├── Chr3.joint.merged.vcf.gz.csi   # 染色体3索引文件
└── ...
```

## 使用示例 | Usage Examples

### 示例1：基本合并流程 | Example 1: Basic Merge Pipeline

```bash
# 合并按染色体分割的GTX联合基因型VCF文件
biopytools vcf-merger \
    -i /data/gtx_output/joint_vcf \
    -o /data/merged_vcf
```

### 示例2：大规模并行处理 | Example 2: Large-Scale Parallel Processing

```bash
# 使用16个线程处理大规模VCF文件
biopytools vcf-merger \
    -i /data/large_scale_vcf \
    -o /data/merged_output \
    -t 16 \
    -v
```

### 示例3：自定义文件模式 | Example 3: Custom File Pattern

```bash
# 处理不同命名格式的VCF文件
biopytools vcf-merger \
    -i /data/custom_vcf \
    -o /data/merged_custom \
    --pattern "*.vcf.gz"
```

### 示例4：快速合并（无索引）| Example 4: Fast Merge (No Index)

```bash
# 快速合并，不创建索引文件
biopytools vcf-merger \
    -i /data/split_vcf \
    -o /data/quick_merge \
    --no-index

# 后续可以单独创建索引
for file in /data/quick_merge/*.vcf.gz; do
    bcftools index "$file"
done
```

### 示例5：详细调试模式 | Example 5: Verbose Debug Mode

```bash
# 使用详细输出模式查看处理细节
biopytools vcf-merger \
    -i /data/vcf_files \
    -o /data/merged_vcf \
    -vv
```

## 系统要求 | System Requirements

### 依赖软件 | Dependencies

- **bcftools** (版本 1.10+)
  - 用于VCF文件合并和索引创建
  - 安装命令: `conda install -c bioconda bcftools`

- **Python** (版本 3.7+)
- **Python包**:
  - `click` - 命令行界面
  - `pathlib` - 路径处理

### 安装依赖软件 | Installing Dependencies

```bash
# 使用conda安装bcftools
conda install -c bioconda bcftools

# 或使用apt-get (Ubuntu/Debian)
sudo apt-get install bcftools

# 或使用yum (CentOS/RHEL)
sudo yum install bcftools

# 验证安装
bcftools --version
```

## 程序输出说明 | Program Output Explanation

### 日志输出示例 | Log Output Example

```
[2025-12-30 10:30:15] INFO: ============================================================
[2025-12-30 10:30:15] INFO: VCF文件按染色体合并工具 | VCF Merge by Chromosome Tool
[2025-12-30 10:30:15] INFO: Version: 1.0.0
[2025-12-30 10:30:15] INFO: ============================================================
[2025-12-30 10:30:15] INFO: 输入目录 | Input directory: /data/vcf_files
[2025-12-30 10:30:15] INFO: 输出目录 | Output directory: /data/merged_vcf
[2025-12-30 10:30:15] INFO: 文件模式 | File pattern: *.joint.vcf.gz
[2025-12-30 10:30:15] INFO: 线程数 | Threads: 8
[2025-12-30 10:30:15] INFO: ============================================================
[2025-12-30 10:30:15] INFO: bcftools已安装 | bcftools is available: bcftools 1.18
[2025-12-30 10:30:15] INFO: ============================================================
[2025-12-30 10:30:15] INFO: STEP 1: 查找VCF文件 | Finding VCF files
[2025-12-30 10:30:15] INFO: ============================================================
[2025-12-30 10:30:15] INFO: 找到 50 个VCF文件 | Found 50 VCF files
[2025-12-30 10:30:15] INFO: ============================================================
[2025-12-30 10:30:15] INFO: STEP 2: 按染色体分组 | Grouping by chromosome
[2025-12-30 10:30:15] INFO: ============================================================
[2025-12-30 10:30:15] INFO: 识别出 10 个染色体 | Identified 10 chromosomes
[2025-12-30 10:30:15] INFO:   Chr1: 5 个文件 | files
[2025-12-30 10:30:15] INFO:   Chr2: 5 个文件 | files
...
[2025-12-30 10:30:15] INFO: ============================================================
[2025-12-30 10:30:15] INFO: STEP 3: 合并VCF文件 | Merging VCF files
[2025-12-30 10:30:15] INFO: ============================================================
[2025-12-30 10:30:15] INFO: 处理染色体 Chr1 | Processing chromosome Chr1
[2025-12-30 10:30:16] INFO: 合并 5 个文件到 Chr1.joint.merged.vcf.gz | Merging 5 files to Chr1.joint.merged.vcf.gz
[2025-12-30 10:30:16] INFO: 成功创建 Chr1.joint.merged.vcf.gz | Successfully created Chr1.joint.merged.vcf.gz
[2025-12-30 10:30:16] INFO: 正在创建索引 | Creating index for: Chr1.joint.merged.vcf.gz
[2025-12-30 10:30:16] INFO: 索引创建成功 | Index created successfully
...
[2025-12-30 10:30:20] INFO: ============================================================
[2025-12-30 10:30:20] INFO: 合并完成 | Merge completed
[2025-12-30 10:30:20] INFO: ============================================================
[2025-12-30 10:30:20] INFO: 成功合并 | Successfully merged: 10/10 个染色体
[2025-12-30 10:30:20] INFO: 输出目录 | Output directory: /data/merged_vcf
[2025-12-30 10:30:20] INFO: 总耗时 | Total runtime: 5.23 秒 | seconds
[2025-12-30 10:30:20] INFO: 所有操作完成 | All operations completed
```

## 常见问题 | Frequently Asked Questions

### Q1: 如何处理命名不规则的VCF文件？

如果文件名不符合自动识别规则，建议重命名文件：

```bash
# 使用rename命令批量重命名
rename 's/^/Chr/' *.vcf.gz

# 或使用循环重命名
for file in chr_*.vcf.gz; do
    mv "$file" "C${file#chr}"
done
```

### Q2: bcftools合并时出现内存不足怎么办？

可以减少线程数或分批处理：

```bash
# 减少线程数
biopytools vcf-merger -i input -o output -t 2

# 或分染色体手动合并
biopytools vcf-merger -i input/Chr1 -o output/Chr1
biopytools vcf-merger -i input/Chr2 -o output/Chr2
```

### Q3: 如何验证合并后的VCF文件？

```bash
# 使用bcftools验证
bcftools view output/Chr1.joint.merged.vcf.gz | head -100

# 统计变异数量
bcftools view output/Chr1.joint.merged.vcf.gz | wc -l

# 检查索引
bcftools index -s output/Chr1.joint.merged.vcf.gz
```

## 故障排除 | Troubleshooting

### 常见错误 | Common Errors

**错误：未找到bcftools**
```
未找到bcftools，请先安装
```
**解决方案**：
```bash
conda install -c bioconda bcftools
```

**错误：未找到匹配的VCF文件**
```
未找到匹配的VCF文件 | No VCF files found matching pattern: *.joint.vcf.gz
```
**解决方案**：
- 检查输入目录路径是否正确
- 使用 `--pattern` 参数指定正确的文件模式
- 确认文件确实存在于目录中

**错误：无法从文件名提取染色体编号**
```
无法从文件名提取染色体编号 | Cannot extract chromosome from filename
```
**解决方案**：
- 检查文件名格式是否符合命名规则
- 考虑重命名文件以符合标准格式

## API使用 | Python API Usage

```python
from biopytools.vcf_merger import VCFMerger

# 创建合并器
merger = VCFMerger(
    input_dir="/path/to/vcf_files",
    output_dir="/path/to/output",
    pattern="*.joint.vcf.gz",
    threads=8,
    create_index=True
)

# 运行合并
success = merger.run()

# 获取统计信息
stats = merger.get_statistics()
print(f"总染色体数: {stats['total_chromosomes']}")
print(f"成功合并: {stats['success_count']}")
print(f"失败数量: {stats['failed_count']}")

# 获取染色体分组信息
chr_groups = merger.get_chromosome_groups()
for chr_id, files in chr_groups.items():
    print(f"{chr_id}: {len(files)} 个文件")
```

## 相关资源 | Related Resources

- [bcftools官方文档](http://www.htslib.org/doc/bcftools.html)
- [VCF格式规范](https://samtools.github.io/hts-specs/VCFv4.3.pdf)
- [GTX联合基因型分析](../docs/gtx_joint.md)

## 许可证 | License

本项目采用MIT许可证 - 详见 [LICENSE](../LICENSE) 文件
