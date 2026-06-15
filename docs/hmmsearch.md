# HMMsearch 分析工具

## 功能概述|Overview

HMMsearch分析工具用于运行hmmsearch并处理结果，生成易读的CSV/Excel表格，并提取匹配的蛋白序列和domain序列。

**两种工作模式**：
- **模式1**: 处理已有的hmmsearch domtblout输出文件
- **模式2**: 运行hmmsearch并自动处理结果

## 主要特性|Key Features

- 运行hmmsearch进行蛋白domain搜索
- 解析hmmsearch domtblout格式输出
- 生成易读的CSV和Excel表格
- 支持E-value和分数过滤（运行时和结果后处理）
- 提取匹配的蛋白序列
- 提取domain序列
- 生成统计摘要

## 快速开始|Quick Start

### 模式1: 处理已有domtblout文件（不指定-p）

```bash
# 基本用法：处理已有domtblout文件
biopytools hmmsearch -i results.domtblout

# 使用E-value过滤结果
biopytools hmmsearch -i results.domtblout -e 1e-10
```

### 模式2: 运行hmmsearch并处理结果（指定-p）

```bash
# 基本用法：运行hmmsearch
biopytools hmmsearch -i NB-ARC.hmm -p proteins.fa

# 使用cut_tc阈值
biopytools hmmsearch -i NB-ARC.hmm -p proteins.fa --cut-tc

# 使用自定义E-value阈值
biopytools hmmsearch -i NB-ARC.hmm -p proteins.fa --evalue-cutoff 1e-5

# 使用多个线程
biopytools hmmsearch -i NB-ARC.hmm -p proteins.fa -t 24
```

## 参数说明|Parameters

### 工作模式选择|Working Mode Selection

| 参数 | 描述 | 示例 |
|------|------|------|
| `-i, --input` | 主输入文件（自动判断模式）<br>不指定-p：模式1（处理domtblout）<br>指定-p：模式2（运行hmmsearch） | `-i results.domtblout`<br>`-i NB-ARC.hmm -p proteins.fa` |
| `-p, --protein-fasta` | 蛋白序列FASTA文件<br>模式2必需，模式1可选（提取序列时需要） | `-p proteins.fa` |

### 输出配置|Output Configuration

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-o, --output-dir` | `./hmmsearch_output` | 输出目录 |
| `--output-prefix` | `hmmsearch_results` | 输出文件前缀 |

### HMMsearch软件配置|HMMsearch Software Configuration

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--hmmsearch-path` | `/share/org/YZWL/yzwl_lixg/miniforge3/envs/resistify_v.1.3.0/bin/hmmsearch` | hmmsearch程序路径 |
| `-t, --threads` | `12` | 线程数 |

### HMMsearch运行参数|HMMsearch Run Parameters (模式2)

| 参数 | 描述 |
|------|------|
| `--evalue-cutoff` | hmmsearch报告E-value阈值 |
| `--score-cutoff` | hmmsearch报告分数阈值 |
| `--cut-tc` | 使用模型的TC trusted cutoff |
| `--cut-ga` | 使用模型的GA gathering cutoff |
| `--cut-nc` | 使用模型的NC noise cutoff |

### 结果过滤参数|Result Filtering Parameters (后处理)

| 参数 | 描述 |
|------|------|
| `-e, --evalue-threshold` | Domain E-value阈值（保留小于等于该值的） |
| `-s, --score-threshold` | Domain分数阈值（保留大于等于该值的） |

### 序列提取选项|Sequence Extraction Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--no-extract-proteins` | False | 不提取蛋白序列 |
| `--no-extract-domains` | False | 不提取domain序列 |

### 输出格式选项|Output Format Options

| 参数 | 描述 |
|------|------|
| `--no-csv` | 不输出CSV文件 |
| `--no-excel` | 不输出Excel文件 |

## 输出文件|Output Files

```
hmmsearch_output/
├── hmmsearch_results.domtblout      # hmmsearch原始输出
├── hmmsearch_results.csv            # 结果表格CSV格式
├── hmmsearch_results.xlsx           # 结果表格Excel格式
├── hmmsearch_results_proteins.fa    # 匹配的蛋白序列
├── hmmsearch_results_domains.fa     # Domain序列
└── hmmsearch_analysis.log           # 运行日志
```

## 表格列说明|Table Columns

输出的CSV/Excel文件包含以下列：

| 列名 | 描述 |
|------|------|
| target_name | 目标基因名称 |
| target_accession | 目标基因登录号 |
| target_length | 目标基因长度 |
| query_name | 查询（HMM）名称 |
| query_accession | 查询登录号 |
| query_length | 查询长度 |
| full_evalue | 全序列E-value |
| full_score | 全序列分数 |
| full_bias | 全序列bias |
| domain_number | Domain编号 |
| domain_total | Domain总数 |
| domain_evalue | Domain E-value |
| domain_score | Domain分数 |
| domain_bias | Domain bias |
| hmm_from | HMM模型起始位置 |
| hmm_to | HMM模型终止位置 |
| align_from | 比对起始位置 |
| align_to | 比对终止位置 |
| env_from | Envelope起始位置 |
| env_to | Envelope终止位置 |
| accuracy | 准确度 |
| description | 描述 |

## 使用示例|Usage Examples

### 示例1：处理已有domtblout文件（模式1）

```bash
# 只解析结果，不提取序列
biopytools hmmsearch -i results.domtblout -o output

# 解析结果并提取序列
biopytools hmmsearch -i results.domtblout -p proteins.fa -o output
```

### 示例2：运行hmmsearch分析NB-ARC domain（模式2）

```bash
# 使用TC cutoff
biopytools hmmsearch -i PF00931.hmm -p Ccu.pep.fa --cut-tc -t 12
```

### 示例3：使用自定义E-value阈值

```bash
# 运行时E-value阈值
biopytools hmmsearch -i PF00931.hmm -p Ccu.pep.fa --evalue-cutoff 1e-5

# 结果后处理E-value阈值
biopytools hmmsearch -i results.domtblout -p proteins.fa -e 1e-10
```

### 示例4：只提取domain序列，不提取完整蛋白

```bash
biopytools hmmsearch -i PF00931.hmm -p Ccu.pep.fa --no-extract-proteins
```

### 示例5：使用GA gathering cutoff

```bash
biopytools hmmsearch -i PF00931.hmm -p Ccu.pep.fa --cut-ga
```

## 注意事项|Important Notes

1. **模式自动判断**: 通过是否指定`-p`参数自动判断工作模式
   - 不指定`-p`：模式1（处理domtblout文件）
   - 指定`-p`：模式2（运行hmmsearch）
2. **蛋白序列**: 模式2必需，模式1可选（仅在提取序列时需要）
3. **Cut选项**: --cut-tc、--cut-ga、--cut-nc互斥，只能选择一个
4. **断点续传**: 如果domtblout文件已存在，会跳过hmmsearch运行
5. **阈值区别**:
   - `--evalue-cutoff`: hmmsearch运行时的报告阈值
   - `-e/--evalue-threshold`: 结果后处理过滤阈值
6. **线程数**: 默认12，遵循开发规范

## 完整工作流程示例|Complete Workflow Example

```bash
# 运行hmmsearch分析NB-ARC domain，使用TC cutoff，12线程
biopytools hmmsearch \
    -i /path/to/PF00931.hmm \
    -p /path/to/proteins.fa \
    --cut-tc \
    -t 12 \
    -o NB_ARC_analysis

# 查看结果
# - 输出表格：NB_ARC_analysis/hmmsearch_results.csv/.xlsx
# - 蛋白序列：NB_ARC_analysis/hmmsearch_results_proteins.fa
# - Domain序列：NB_ARC_analysis/hmmsearch_results_domains.fa
# - 统计日志：NB_ARC_analysis/hmmsearch_analysis.log
```
