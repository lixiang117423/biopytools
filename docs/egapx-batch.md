# EGAPx批量运行配置生成工具

## 功能概述

`biopytools egapx-batch` 是一个专门用于批量生成EGAPx运行配置和脚本的工具。它可以将大型基因组按染色体或scaffold拆分，为每个序列生成独立的配置文件和执行脚本，支持并行化批量基因预测和注释。

## 核心功能

- 📂 **按序列拆分基因组** - 自动按染色体/scaffold拆分
- 📝 **生成独立配置** - 为每个序列生成YAML配置文件
- 📜 **生成运行脚本** - 为每个序列生成执行脚本
- 🔗 **自动软链接** - 自动创建EGAPx程序软链接
- 🚀 **并行执行支持** - 生成任务列表和并行执行脚本
- 🏷️ **自定义命名** - 支持自定义locus标签前缀和报告名称
- 🔍 **智能过滤** - 支持按前缀过滤特定染色体

## 使用方法

### 基本用法

```bash
biopytools egapx-batch \
    -g genome.fa \
    -y template.yaml \
    -s template.sh \
    -o output_dir
```

### 自定义locus标签前缀和报告名

```bash
biopytools egapx-batch \
    -g genome.fa \
    -y template.yaml \
    -s template.sh \
    -o output_dir \
    --locus-prefix Gene \
    --report-name MyEGAPx
```

### 只处理主染色体（Chr开头）

```bash
biopytools egapx-batch \
    -g genome.fa \
    -y template.yaml \
    -s template.sh \
    -o output_dir \
    --chr-prefix Chr
```

### 自定义EGAPx安装路径

```bash
biopytools egapx-batch \
    -g genome.fa \
    -y template.yaml \
    -s template.sh \
    -o output_dir \
    --egapx /custom/path/to/egapx
```

## 命令行参数

### 必需参数

| 参数 | 简写 | 说明 | 示例 |
|------|------|------|------|
| `--genome` | `-g` | 基因组FASTA文件路径 | `genome.fa` |
| `--yaml` | `-y` | YAML模板文件路径 | `template.yaml` |
| `--script` | `-s` | Shell脚本模板路径 | `template.sh` |
| `--output` | `-o` | 输出目录路径 | `output_dir` |

### 可选参数

| 参数 | 简写 | 默认值 | 说明 |
|------|------|--------|------|
| `--egapx` | `-e` | `/share/org/YZWL/yzwl_lixg/software/EGAPX_v.0.4.1-alpha/egapx` | EGAPx安装路径 |
| `--chr-prefix` | `-p` | None | 染色体前缀过滤 |
| `--locus-prefix` | | `Target` | locus标签前缀 |
| `--report-name` | | `EGAPx` | 报告名称 |

## 输出目录结构

```
output_dir/
├── Chr01/
│   ├── Chr01.fa              # 染色体序列文件
│   ├── Chr01.yaml            # YAML配置文件
│   ├── egapx_Chr01.sh        # 运行脚本
│   ├── work/                 # 工作目录
│   ├── output/               # 输出目录
│   └── [EGAPx软链接文件]     # EGAPx程序软链接
├── Chr02/
│   └── ... (同上结构)
├── Chr03/
│   └── ... (同上结构)
├── all_jobs_submit.list.sh   # 所有任务列表
├── run_all_parallel.sh       # 并行执行脚本
└── egapx_batch.log           # 运行日志
```

## 模板文件要求

### YAML模板 (template.yaml)

```yaml
genome: /path/to/genome.fa
taxid: 71234
short_reads: /path/to/short_reads.txt
long_reads: /path/to/long_reads.txt
locus_tag_prefix: Target
```

**字段说明：**
- `genome`: 基因组文件路径（会被自动替换）
- `taxid`: 物种分类ID
- `short_reads`: 短读长测序数据列表文件
- `long_reads`: 长读长测序数据列表文件
- `locus_tag_prefix`: locus标签前缀（会被自动修改为前缀_序列名）

### Shell脚本模板 (template.sh)

```bash
#!/bin/bash
# 激活conda环境
source ~/.bashrc
conda activate base

export JAVA_HOME=/path/to/java
export PATH=$JAVA_HOME/bin:$PATH

# 运行EGAPx
python3 \
    ui/egapx.py \
    /path/to/config.yaml \
    -e singularity \
    -w /path/to/work \
    -o /path/to/output \
    -r EGAPx
```

**参数说明：**
- `ui/egapx.py`: EGAPx主程序
- YAML配置文件路径（会被自动替换）
- `-e`: 执行环境（singularity/docker）
- `-w`: 工作目录（会被自动替换）
- `-o`: 输出目录（会被自动替换）
- `-r`: 报告名称（会被自动修改为报告名_序列名）

## 执行方式

### 1. 顺序执行

```bash
bash output_dir/all_jobs_submit.list.sh
```

适合资源有限或需要按顺序处理的情况。

### 2. 并行执行（推荐）

```bash
# 并行4个任务
bash output_dir/run_all_parallel.sh 4

# 并行8个任务
bash output_dir/run_all_parallel.sh 8
```

### 3. 使用GNU Parallel

```bash
# 并行4个任务
cat output_dir/all_jobs_submit.list.sh | parallel -j 4

# 并行8个任务
cat output_dir/all_jobs_submit.list.sh | parallel -j 8
```

### 4. 使用xargs

```bash
# 并行4个任务
cat output_dir/all_jobs_submit.list.sh | xargs -P 4 -I {} bash -c {}
```

## 使用流程示例

### 步骤1：准备模板文件

复制默认模板并修改：

```bash
# 复制YAML模板
cp ~/software/scripts/62.EGAPx_temple.yaml my_template.yaml

# 复制Shell脚本模板
cp ~/software/scripts/63.EGAPx_temple.sh my_template.sh

# 编辑模板文件，填入正确的路径
vim my_template.yaml
vim my_template.sh
```

### 步骤2：生成批量配置

```bash
biopytools egapx-batch \
    -g genome.fa \
    -y my_template.yaml \
    -s my_template.sh \
    -o egapx_batch_output
```

### 步骤3：执行批量任务

```bash
# 并行执行
bash egapx_batch_output/run_all_parallel.sh 4
```

### 步骤4：收集结果

每个染色体的结果在对应的output目录中：

```bash
# 查看所有结果
ls -lh egapx_batch_output/*/output/

# 合并GFF文件
cat egapx_batch_output/*/output/*.gff > all_genes.gff

# 合并统计信息
cat egapx_batch_output/*/output/*summary.txt > all_summary.txt
```

## 应用场景

### 1. 大型基因组注释

大型基因组（如人类、小麦）需要分染色体处理以降低内存需求：

```bash
biopytools egapx-batch \
    -g large_genome.fa \
    -y template.yaml \
    -s template.sh \
    -o annotation_results
```

### 2. 多样本比较分析

对不同样本的同一染色体分别进行注释：

```bash
# 样本1
biopytools egapx-batch \
    -g sample1.fa \
    -y template.yaml \
    -s template.sh \
    -o sample1_results \
    --locus-prefix Sample1

# 样本2
biopytools egapx-batch \
    -g sample2.fa \
    -y template.yaml \
    -s template.sh \
    -o sample2_results \
    --locus-prefix Sample2
```

### 3. 特定区域分析

只分析特定染色体或scaffold：

```bash
# 只分析主染色体
biopytools egapx-batch \
    -g genome.fa \
    -y template.yaml \
    -s template.sh \
    -o main_chr_results \
    --chr-prefix Chr
```

## 故障排除

### 1. "基因组文件不存在"

```
错误: 基因组文件不存在: genome.fa
```

**解决方案：**
- 检查`-g`参数指定的文件路径
- 确认文件格式为`.fa`, `.fasta`, `.fa.gz`或`.fasta.gz`

### 2. "未拆分出任何序列"

```
错误: 未拆分出任何序列
```

**可能原因：**
- 使用了错误的`--chr-prefix`
- 基因组文件格式不正确

**解决方案：**
```bash
# 不使用chr-prefix参数
biopytools egapx-batch -g genome.fa -y template.yaml -s template.sh -o output

# 或检查染色体ID格式
grep "^>" genome.fa | head -10
```

### 3. "YAML模板不存在"

```
错误: YAML模板不存在: template.yaml
```

**解决方案：**
- 使用绝对路径
- 参考默认模板：`~/software/scripts/62.EGAPx_temple.yaml`

### 4. "脚本模板不存在"

```
错误: 脚本模板不存在: template.sh
```

**解决方案：**
- 使用绝对路径
- 参考默认模板：`~/software/scripts/63.EGAPx_temple.sh`

### 5. "EGAPx路径不存在"

```
错误: EGAPx路径不存在: /path/to/egapx
```

**解决方案：**
- 使用`--egapx`参数指定正确路径
- 或确保默认路径存在：`/share/org/YZWL/yzwl_lixg/software/EGAPX_v.0.4.1-alpha/egapx`

## 最佳实践

### 1. 模板文件准备

- 使用绝对路径避免路径问题
- 确保模板中的路径可以被正确替换
- 在测试数据上验证模板的正确性

### 2. 命名规范

- `locus_tag_prefix`使用有意义的名称：
  - `Target` - 目标基因
  - `Gene` - 通用基因
  - `Sample1` - 特定样本
  - 项目名称

- `report_name`描述项目：
  - `EGAPx` - 通用
  - `MyProject` - 项目名称
  - `SampleAnnotation` - 样本注释

### 3. 染色体过滤

使用`--chr-prefix`只处理需要的序列：

| 示例 | 效果 |
|------|------|
| `--chr-prefix Chr` | 只处理Chr01, Chr02... |
| `--chr-prefix NC` | 只处理NC开头的scaffold |
| `--chr-prefix scaffold_` | 只处理scaffold_开头的序列 |

### 4. 并行执行策略

- 根据可用内存和CPU选择并行度
- 每个EGAPx任务通常需要16-32GB内存
- 推荐并行度：可用内存 / 单任务内存需求

```bash
# 示例：128GB内存，每个任务需要32GB
bash run_all_parallel.sh 4  # 128/32 = 4
```

### 5. 结果整理

EGAPx运行完成后，需要合并结果：

```bash
# 合并GFF3基因注释
cat */output/*.gff3 > all_genes.gff3

# 合并GBK文件
cat */output/*.gbk > all_genes.gbk

# 提取基因统计
for dir in Chr*/; do
    echo "=== $dir ==="
    grep "Total genes" "$dir"output/*summary.txt
done
```

## 技术细节

### 拆分原理

工具使用awk按FASTA序列头拆分基因组：

1. 读取基因组FASTA文件
2. 遇到`>`开头的行时，提取序列ID
3. 清理序列ID中的特殊字符（`:`, `|`, `/`, `\`替换为`_`）
4. 应用前缀过滤（如果指定）
5. 将序列写入独立的文件

### 路径替换逻辑

脚本会智能替换以下路径：

1. **YAML配置中的路径**：
   - `genome`: 替换为对应的染色体序列文件
   - `locus_tag_prefix`: 替换为前缀_序列名

2. **Shell脚本中的路径**：
   - YAML文件路径: 替换为新的YAML配置
   - `-w` 工作目录: 替换为对应染色体的work目录
   - `-o` 输出目录: 替换为对应染色体的output目录
   - `-r` 报告名称: 替换为报告名_序列名

### 软链接管理

为每个染色体目录创建EGAPx程序的软链接：

- 方便独立运行
- 避免路径问题
- 执行后自动清理

## 版本历史

- **v1.0.0** (2025-12-25)
  - 初始版本发布
  - 支持按序列拆分基因组
  - 支持生成独立配置和脚本
  - 支持并行执行脚本生成
  - 支持自定义命名和过滤

## 相关工具

- `biopytools gtx-joint` - GTX Joint Calling命令生成
- `biopytools annovar` - ANNOVAR变异注释
- `biopytools vcf-renamer` - VCF样品名称重命名
