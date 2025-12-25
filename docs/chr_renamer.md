# 染色体重命名工具 (chr-renamer)

## 功能概述

染色体重命名工具是一个专门用于将FASTA文件序列重命名为标准格式的工具。它可以将输入文件中的前N条序列命名为染色体格式（如Chr01, Chr02...），剩余序列命名为scaffold格式（如HiC_scaffold_01, HiC_scaffold_02...）。

## 主要特性

- ✅ **高效流式处理**: 使用Python生成器实现内存高效的大文件处理
- ✅ **灵活命名规则**: 支持自定义前缀和格式化规则
- ✅ **压缩文件支持**: 自动识别并处理gzip压缩的FASTA文件
- ✅ **全面验证**: 内置FASTA格式验证和输入文件检查
- ✅ **详细统计**: 提供处理过程的详细统计信息
- ✅ **备份功能**: 可选的输出文件备份机制
- ✅ **标准日志**: 遵循生信工具规范的日志输出格式

## 安装要求

该工具是BioPyTools的一部分，无需额外安装依赖。

## 基本用法

### 命令行接口

```bash
# 基本重命名
biopytools chr-renamer -i input.fa -o output.fa -n 20

# 自定义前缀
biopytools chr-renamer \
  -i scaffolds.fa \
  -o renamed.fa \
  -n 10 \
  --chr-prefix Chr \
  --scaffold-prefix Scf

# 处理压缩文件
biopytools chr-renamer \
  -i input.fa.gz \
  -o output.fa.gz \
  -n 30 \
  --backup

# 跳过验证和预览（快速处理）
biopytools chr-renamer \
  -i large_genome.fa \
  -o renamed.fa \
  -n 25 \
  --no-validate \
  --preview-lines 0
```

### Python API

```python
from biopytools.chr_renamer import ChromosomeRenamer

# 创建重命名器
renamer = ChromosomeRenamer(
    input_file="input.fa",
    output_file="output.fa",
    chromosome_count=20
)

# 运行重命名
success = renamer.run()

# 获取统计信息
if success:
    stats = renamer.get_statistics()
    print(f"处理了 {stats['total_sequences']} 条序列")
    print(f"染色体: {stats['chr_sequences']}")
    print(f"Scaffolds: {stats['scaffold_sequences']}")
```

## 参数说明

### 必需参数

| 参数 | 描述 | 示例 |
|------|------|------|
| `-i, --input` | 输入FASTA文件路径 | `input.fa` |
| `-o, --output` | 输出FASTA文件路径 | `output.fa` |
| `-n, --number` | 染色体数量 | `20` |

### 可选参数

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--chr-prefix` | `Chr` | 染色体前缀 |
| `--scaffold-prefix` | `HiC_scaffold` | Scaffold前缀 |
| `--chr-format` | `{:02d}` | 染色体编号格式 |
| `--scaffold-format` | `{:02d}` | Scaffold编号格式 |
| `--no-validate` | `False` | 跳过输入文件验证 |
| `--backup` | `False` | 备份输出文件（如果存在） |
| `--preview-lines` | `5` | 预览行数（0=不预览） |

## 输出格式

### 命名规则

- **染色体格式**: `{chr_prefix}{chr_format}`
  - 默认: `Chr01`, `Chr02`, `Chr03`...
- **Scaffold格式**: `{scaffold_prefix}_{scaffold_format}`
  - 默认: `HiC_scaffold_01`, `HiC_scaffold_02`...

### 输出示例

```
>Chr01
ATCGATCGATCGATCG...
>Chr02
GCTAGCTAGCTAGCTA...
>HiC_scaffold_01
TTTTAAAAATTTTAAA...
>HiC_scaffold_02
CCCCGGGGCCCCGGGG...
```

## 日志输出

### INFO输出（stdout）
```
[2025-12-21 10:30:15] INFO: 🧬 染色体重命名工具启动 | Chromosome Renamer Tool Started
[2025-12-21 10:30:15] INFO: 版本 | Version: 1.0.0
[2025-12-21 10:30:15] INFO: ============================================================
[2025-12-21 10:30:15] INFO: STEP: 输入文件验证 | Input File Validation
[2025-12-21 10:30:15] INFO: ============================================================
[2025-12-21 10:30:15] INFO: 输入文件 | Input file: input.fa
[2025-12-21 10:30:15] INFO: 输出文件 | Output file: output.fa
[2025-12-21 10:30:15] INFO: 染色体数量 | Chromosome count: 20
[2025-12-21 10:30:15] INFO: 输入序列数量 | Input sequence count: 150
[2025-12-21 10:30:15] INFO: FASTA文件验证通过 | FASTA file validation passed: 150 sequences
```

### 统计摘要
```
[2025-12-21 10:30:20] INFO: ============================================================
[2025-12-21 10:30:20] INFO: 处理摘要 | Processing Summary
[2025-12-21 10:30:20] INFO: ============================================================
[2025-12-21 10:30:20] INFO: 总序列数量 | Total sequences: 150
[2025-12-21 10:30:20] INFO: 染色体数量 | Chromosomes: 20
[2025-12-21 10:30:20] INFO: Scaffold数量 | Scaffolds: 130
[2025-12-21 10:30:20] INFO: 染色体比例 | Chromosome percentage: 13.3%
[2025-12-21 10:30:20] INFO: Scaffold比例 | Scaffold percentage: 86.7%
[2025-12-21 10:30:20] INFO: 输出文件大小 | Output file size: 1,234,567,890 bytes
```

## 高级用法

### 自定义命名格式

```python
from biopytools.chr_renamer import ChromosomeRenamer

# 使用罗马数字格式
renamer = ChromosomeRenamer(
    input_file="genome.fa",
    output_file="roman_genome.fa",
    chromosome_count=12,
    chr_prefix="Chromosome_",
    chr_format="{:03d}",  # Chromosome_001, Chromosome_002...
    scaffold_prefix="Scaffold_",
    scaffold_format="{:04d}"  # Scaffold_0001, Scaffold_0002...
)
```

### 处理超大型文件

```python
# 对于超大型文件，可以跳过验证以提高速度
renamer = ChromosomeRenamer(
    input_file="huge_genome.fa",
    output_file="renamed.fa",
    chromosome_count=100,
    validate_input=False  # 跳过FASTA格式验证
)

# 不预览输出，直接处理
renamer.run(preview_output=False)
```

### 错误处理

```python
try:
    renamer = ChromosomeRenamer(
        input_file="nonexistent.fa",
        output_file="output.fa",
        chromosome_count=10
    )
    success = renamer.run()

    if not success:
        print("处理失败，请检查日志")

except FileNotFoundError as e:
    print(f"文件未找到: {e}")
except ValueError as e:
    print(f"参数错误: {e}")
```

## 性能优化建议

1. **大文件处理**: 对于大型基因组文件，使用`--no-validate`选项跳过格式验证
2. **内存优化**: 工具使用流式处理，内存占用与文件大小无关
3. **并行处理**: 如需处理多个文件，可以并行运行多个实例
4. **压缩处理**: 工具自动支持gzip压缩，无需额外解压步骤

## 故障排除

### 常见错误

1. **文件不存在**
   ```
   错误: 输入文件 'input.fa' 不存在
   ```
   解决方案: 检查文件路径是否正确

2. **格式错误**
   ```
   警告: 序列 5 包含无效字符: X*
   ```
   解决方案: 检查FASTA文件中的序列是否只包含标准核苷酸字符

3. **染色体数量过多**
   ```
   警告: 染色体数量 (100) 大于输入序列数量 (50)
   ```
   解决方案: 调整染色体数量或检查输入文件

### 调试技巧

1. 使用`-v`选项启用详细输出
2. 使用`--preview-lines 10`查看更多输出预览
3. 检查日志文件中的WARNING和ERROR信息

## 版本历史

- **v1.0.0** (2025-12-21): 初始版本
  - 基本重命名功能
  - 完整的日志和统计功能
  - CLI和Python API接口

## 相关工具

- [annovar](annovar.md): ANNOVAR变异注释工具
- [fastp](fastp.md): FASTQ数据质量控制
- [bwa](bwa.md): 全基因组比对工具

## 许可证

本工具遵循BioPyTools项目的许可证条款。

## 贡献

欢迎提交Issue和Pull Request来改进这个工具。