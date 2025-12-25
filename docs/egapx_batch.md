# EGAPx批量处理工具 (egapx-batch)

## 功能概述

EGAPx批量处理工具是一个专门用于按染色体拆分基因组并生成批量EGAPx处理配置的工具。它可以将大型基因组文件按染色体序列拆分为独立文件，为每个染色体自动生成对应的YAML配置文件和Shell运行脚本，支持路径动态替换和软链接管理，便于并行处理大量基因预测任务。

## 主要特性

- ✅ **智能文件分割**: 按染色体自动拆分大型FASTA基因组文件
- ✅ **模板配置处理**: 自动提取并替换YAML和Shell脚本模板中的路径参数
- ✅ **动态路径替换**: 智能识别并替换配置文件中的各种路径和参数
- ✅ **软链接管理**: 自动为每个任务创建和管理EGAPx软链接
- ✅ **批量脚本生成**: 生成顺序执行和并行执行的任务脚本
- ✅ **通用工具支持**: 不仅限于EGAPx，支持其他类似的基因预测工具
- ✅ **详细日志记录**: 提供完整的处理过程日志和统计信息
- ✅ **输入验证**: 全面的文件格式和配置验证

## 安装要求

该工具是BioPyTools的一部分，无需额外安装依赖。

## 内置模板

工具内置了默认的YAML和Shell脚本模板，基于原始脚本中的模板文件：

### 默认YAML模板
```yaml
genome: PLACEHOLDER_GENOME_PATH
taxid: 71234
short_reads: /path/to/short_reads.txt
long_reads: /path/to/long_reads.txt
# locus_tag_prefix: PLACEHOLDER_LOCUS_TAG_PREFIX  # 默认注释，只在指定时启用
```

### 默认Shell脚本模板
```bash
#!/bin/bash
# 激活conda环境
source ~/.bashrc
conda activate base

export JAVA_HOME=/share/org/YZWL/yzwl_lixg/miniforge3/envs/EGAPx_v.0.4.0-alpha
export PATH=$JAVA_HOME/bin:$PATH
export PATH="/share/org/YZWL/yzwl_lixg/software:$PATH"

# 运行EGAPx
python3 \
    ui/egapx.py \
    PLACEHOLDER_YAML_PATH \
    -e singularity \
    -w PLACEHOLDER_WORK_PATH \
    -o PLACEHOLDER_OUTPUT_PATH \
    -lc /share/org/YZWL/yzwl_lixg/software/EGAPX_v.0.4.1-alpha/local_cache \
    -r PLACEHOLDER_REPORT_NAME
```

工具会自动替换占位符（如`PLACEHOLDER_GENOME_PATH`）为实际的路径值。

## 基本用法

### 命令行接口

```bash
# 🎯 使用内置默认模板（推荐）
biopytools egapx-batch \
  -g genome.fa \
  -o batch_output

# 🔧 自定义工具名称和前缀（使用内置模板）
biopytools egapx-batch \
  -g genome.fa \
  -o batch_output \
  --tool-name GENE \
  --locus-tag-prefix GeneTarget \
  --report-prefix gene_report

# 📦 使用自定义模板文件
biopytools egapx-batch \
  -g genome.fa \
  -y config.yaml \
  -s run.sh \
  -o batch_output

# 📦 过滤特定染色体（使用内置模板）
biopytools egapx-batch \
  -g genome.fa \
  -o batch_output \
  --prefix Chr

# ⚡ 跳过验证快速处理
biopytools egapx-batch \
  -g large_genome.fa \
  -o batch_output \
  --no-validate
```

### Python API

```python
from biopytools.egapx_batch import EGAPxBatchProcessor

# 创建批量处理器
processor = EGAPxBatchProcessor(
    genome_file="genome.fa",
    yaml_template="config.yaml",
    script_template="run.sh",
    output_dir="batch_output",
    tool_name="GENE",
    locus_tag_prefix="GeneTarget"
)

# 运行批量处理
success = processor.run()

# 获取统计信息
if success:
    stats = processor.get_statistics()
    print(f"处理了 {stats['statistics']['total_sequences']} 个序列")
    print(f"创建了 {stats['statistics']['created_tasks']} 个任务")
```

## 参数说明

### 必需参数

| 参数 | 描述 | 示例 |
|------|------|------|
| `-g, --genome` | 基因组FASTA文件路径 | `genome.fa` |
| `-o, --output` | 输出目录路径 | `batch_output` |

### 模板参数（可选）

| 参数 | 描述 | 默认值 |
|------|------|--------|
| `-y, --yaml` | YAML配置模板文件路径 | 使用内置模板 |
| `-s, --script` | Shell脚本模板文件路径 | 使用内置模板 |

### 可选参数

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-e, --egapx` | `/share/org/YZWL/yzwl_lixg/software/EGAPX_v.0.4.1-alpha/egapx` | EGAPx安装路径 |
| `-p, --prefix` | `""` | 染色体前缀过滤（如只处理Chr开头的序列） |
| `--locus-tag-prefix` | `""` | locus_tag前缀（可选，默认不设置） |
| `--report-prefix` | `EGAPx` | 报告前缀 |
| `--tool-name` | `EGAPx` | 工具名称（用于生成脚本文件名） |
| `--no-validate` | `False` | 跳过输入文件验证 |

## 模板文件要求

### YAML模板文件示例

```yaml
genome: /path/to/genome.fa
taxid: 71234
short_reads: /path/to/short_reads.txt
long_reads: /path/to/long_reads.txt
locus_tag_prefix: Target
```

### Shell脚本模板示例

```bash
#!/bin/bash
# 激活conda环境
source ~/.bashrc
conda activate base

# 运行EGAPx
python3 \
    ui/egapx.py \
    config.yaml \
    -e singularity \
    -w work \
    -o output \
    -lc local_cache \
    -r EGAPx
```

## 输出结构

工具会创建以下目录结构：

```
batch_output/
├── all_jobs_submit.list.sh      # 任务列表脚本
├── run_all_parallel.sh          # 并行执行脚本
├── Chr01/                       # 染色体1任务目录
│   ├── Chr01.fa                # 染色体1序列文件
│   ├── Chr01.yaml              # 染色体1配置文件
│   ├── egapx_Chr01.sh          # 染色体1运行脚本
│   ├── work/                   # 工作目录
│   ├── output/                 # 输出目录
│   └── [EGAPx软链接...]        # EGAPx工具软链接
├── Chr02/                       # 染色体2任务目录
│   └── ...
└── ...
```

## 输出文件说明

### 生成的配置文件

每个染色体任务会生成：

1. **FASTA文件**: `ChrXX.fa` - 染色体序列
2. **YAML配置**: `ChrXX.yaml` - 自动替换路径和前缀的配置
3. **Shell脚本**: `egapx_ChrXX.sh` - 可独立执行的运行脚本
4. **软链接**: EGAPx工具的软链接，便于独立运行

### 批量执行脚本

1. **任务列表**: `all_jobs_submit.list.sh` - 包含所有任务的执行命令
2. **并行脚本**: `run_all_parallel.sh` - 支持并行执行的脚本

## 执行方式

### 1. 顺序执行

```bash
bash batch_output/all_jobs_submit.list.sh
```

### 2. 并行执行（推荐）

```bash
# 使用生成的并行脚本
bash batch_output/run_all_parallel.sh 4  # 并行4个任务

# 或使用GNU parallel
cat batch_output/all_jobs_submit.list.sh | parallel -j 4

# 或使用xargs
cat batch_output/all_jobs_submit.list.sh | xargs -P 4 -I {} bash -c {}
```

## 日志输出

### INFO输出（stdout）
```
[2025-12-21 10:30:15] INFO: 🧬 EGAPx批量处理工具启动 | EGAPx Batch Processing Tool Started
[2025-12-21 10:30:15] INFO: 版本 | Version: 1.0.0
[2025-12-21 10:30:15] INFO: ============================================================
[2025-12-21 10:30:15] INFO: STEP: 输入文件验证 | Input File Validation
[2025-12-21 10:30:15] INFO: ============================================================
[2025-12-21 10:30:15] INFO: 基因组文件 | Genome file: genome.fa
[2025-12-21 10:30:15] INFO: YAML模板 | YAML template: config.yaml
[2025-12-21 10:30:15] INFO: 待处理序列数量 | Sequences to process: 12
[2025-12-21 10:30:15] INFO: FASTA文件验证通过 | FASTA file validation passed: 12 sequences
```

### 处理摘要
```
[2025-12-21 10:35:20] INFO: ============================================================
[2025-12-21 10:35:20] INFO: 处理摘要 | Processing Summary
[2025-12-21 10:35:20] INFO: ============================================================
[2025-12-21 10:35:20] INFO: 总序列数量 | Total sequences: 12
[2025-12-21 10:35:20] INFO: 处理序列数量 | Processed sequences: 12
[2025-12-21 10:35:20] INFO: 创建任务数量 | Created tasks: 12
[2025-12-21 10:35:20] INFO: 任务创建成功率 | Task creation success rate: 100.0%
```

## 高级用法

### 使用不同工具

```python
from biopytools.egapx_batch import EGAPxBatchProcessor

# 为其他基因预测工具（如MAKER）创建批量任务
processor = EGAPxBatchProcessor(
    genome_file="genome.fa",
    yaml_template="maker_config.yaml",
    script_template="maker_run.sh",
    output_dir="maker_batch",
    tool_name="MAKER",
    locus_tag_prefix="MAKER",
    report_prefix="maker_output"
)
processor.run()
```

### 过滤特定序列

```bash
# 只处理以Chr开头的染色体
biopytools egapx-batch \
  -g genome.fa \
  -y config.yaml \
  -s run.sh \
  -o batch_output \
  --prefix Chr

# 只处理以scaffold_开头的序列
biopytools egapx-batch \
  -g genome.fa \
  -y config.yaml \
  -s run.sh \
  -o batch_output \
  --prefix scaffold_
```

### 自定义配置参数

```python
processor = EGAPxBatchProcessor(
    genome_file="genome.fa",
    yaml_template="config.yaml",
    script_template="run.sh",
    output_dir="custom_batch",
    locus_tag_prefix="MyProject",  # 自定义locus_tag前缀
    report_prefix="my_analysis",   # 自定义报告名称
    tool_name="MYTOOL"             # 自定义工具名称
)
```

## 性能优化建议

1. **大文件处理**: 对于超大型基因组，使用`--no-validate`选项跳过格式验证
2. **并行执行**: 使用并行脚本或GNU parallel提高处理效率
3. **存储优化**: 确保输出目录有足够的磁盘空间
4. **内存管理**: 工具使用流式处理，内存占用与文件大小无关

## 故障排除

### 常见错误

1. **模板文件格式错误**
   ```
   警告: YAML文件中未找到genome路径配置 | genome path not found in YAML file
   ```
   解决方案: 检查YAML模板文件格式，确保包含`genome:`字段

2. **序列过滤无结果**
   ```
   错误: 没有找到符合条件的序列 | No sequences found matching criteria
   ```
   解决方案: 检查前缀过滤条件或去除过滤器

3. **路径替换失败**
   ```
   警告: Shell脚本中未找到YAML路径配置 | YAML path not found in Shell script
   ```
   解决方案: 检查Shell脚本模板中的YAML路径引用

### 调试技巧

1. 使用`-v`选项启用详细输出
2. 检查生成的配置文件是否正确
3. 先用小文件测试模板格式
4. 验证EGAPx安装路径是否正确

## 版本历史

- **v1.0.0** (2025-12-21): 初始版本
  - 基本批量处理功能
  - 支持YAML和Shell脚本模板处理
  - 动态路径替换功能
  - 完整的日志和统计功能
  - 通用工具支持

## 相关工具

- [annovar](annovar.md): ANNOVAR变异注释工具
- [chr_renamer](chr_renamer.md): 染色体重命名工具
- [fastp](fastp.md): FASTQ数据质量控制

## 许可证

本工具遵循BioPyTools项目的许可证条款。

## 贡献

欢迎提交Issue和Pull Request来改进这个工具。