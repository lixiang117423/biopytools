# K-mer工具集 (KmerTools)

**基于kmtricks和RocksDB的kmer分析工具集 | K-mer Analysis Toolkit Based on kmtricks and RocksDB**

## 功能概述 | Overview

KmerTools是一个基于kmtricks和RocksDB的高效kmer分析工具集，提供从kmer库构建、提取、查询到最终矩阵生成和样本筛选的完整流程。

## 主要特性 | Key Features

- **完整流程支持**: build → extract → query → matrix → screener
- **高效查询**: 基于RocksDB的高性能kmer查询，支持反向互补kmer自动检测
- **Seqkit风格**: 采用seqkit风格的子命令设计
- **样本筛选**: 基于连续0数量的智能样本筛选功能
- **可配置参数**: kmer大小、线程数、阈值等参数灵活可配

## 快速开始 | Quick Start

### 安装依赖 | Install Dependencies

```bash
# 安装kmtricks
conda install -c bioconda kmtricks

# 安装python依赖
pip install python-rocksdb pyfastx

# 可选: bgzip (用于压缩矩阵)
conda install -c bioconda htslib
```

### 基本用法 | Basic Usage

#### 1. 构建kmer数据库 | Build Kmer Database

```bash
# 构建kmer库：kmtricks + aggregate + rocksdb导入
biopytools kmertools build -i fastq_dir -o output_dir
```

#### 2. 从FASTA提取kmer | Extract Kmers from FASTA

```bash
biopytools kmertools extract \
    -i gene.fa \
    --kmer-output gene.kmer.txt \
    --kmer-pos-output gene.kmer.pos.txt
```

#### 3. 从RocksDB查询kmer | Query Kmers from RocksDB

```bash
biopytools kmertools query \
    -r ./kmer_rocksdb \
    -i gene.kmer.txt \
    -o gene.query.txt
```

#### 4. 生成矩阵 | Generate Matrix

```bash
biopytools kmertools matrix \
    -i gene.query.txt \
    -p gene.kmer.pos.txt \
    -o gene.matrix.txt
```

#### 5. 筛选样本 | Screen Samples

```bash
biopytools kmertools screener \
    -g gene.txt \
    -s sample.txt \
    -d ./matrices \
    -o ./results
```

## 子命令详解 | Subcommands Details

### 1. build - 构建kmer数据库 | Build Kmer Database

**功能**: 运行kmtricks pipeline、aggregate并导入到RocksDB

**用法**:
```bash
biopytools kmertools build [OPTIONS]
```

**选项**:

| 参数 | 简写 | 默认值 | 描述 |
|------|------|--------|------|
| `--input-dir` | `-i` | 必需 | 输入FASTQ目录 |
| `--fof-file` | | 自动生成 | FOF文件路径 |
| `--kmer-size` | `-k` | 51 | Kmer大小 |
| `--hard-min` | | 2 | 最小丰度 |
| `--recurrence-min` | | 2 | 最小出现次数 |
| `--run-dir` | | output/kmtricks_run | kmtricks运行目录 |
| `--rocksdb-path` | `-r` | output/kmer_rocksdb | RocksDB数据库路径 |
| `--header-file` | | 自动生成 | 样本名文件 |
| `--header-db-key` | | kmer_header | RocksDB header key |

**输出文件**:
```
kmertools_output/
├── input.fof                # FOF文件
├── kmtricks_run/            # kmtricks运行目录
├── kmer_matrix.txt.gz      # 压缩的kmer矩阵
├── header.txt              # 样本名文件
└── kmer_rocksdb/           # RocksDB数据库
```

### 2. extract - 提取kmer | Extract Kmers

**功能**: 从FASTA文件提取kmer和位置信息，生成canonical kmer

**用法**:
```bash
biopytools kmertools extract [OPTIONS]
```

**选项**:

| 参数 | 简写 | 默认值 | 描述 |
|------|------|--------|------|
| `--fasta-file` | `-i` | 必需 | 输入FASTA文件 |
| `--kmer-size` | `-k` | 51 | Kmer大小 |
| `--kmer-output` | | 必需 | 输出kmer列表文件 |
| `--kmer-pos-output` | | 必需 | 输出kmer位置文件 |

**输出文件**:
- `*.kmer.txt`: kmer列表（每行一个canonical kmer）
- `*.kmer.pos.txt`: kmer位置信息（kmer, sample, position）

**特点**:
- 自动生成canonical kmer（取kmer和反向互补序列中字典序较小者）
- 支持多序列FASTA文件
- 记录kmer在序列中的起始位置

### 3. query - 查询kmer | Query Kmers

**功能**: 从RocksDB数据库查询kmer的存在情况

**用法**:
```bash
biopytools kmertools query [OPTIONS]
```

**选项**:

| 参数 | 简写 | 默认值 | 描述 |
|------|------|--------|------|
| `--rocksdb-path` | `-r` | 必需 | RocksDB数据库路径 |
| `--kmer-input` | `-i` | 必需 | 输入kmer文件 |
| `--query-output` | `-o` | 必需 | 查询结果输出文件 |
| `--header-db-key` | | kmer_header | RocksDB header key |
| `--bloom-bits` | | 15 | Bloom filter bits per key |

**输出格式**:
```
ID	sample1	sample2	sample3	...
kmer1	1	0	1	...
kmer2	1	1	1	...
kmer3	0	1	0	...
```

**特点**:
- 自动检测反向互补kmer
- 支持批量查询（multi_get）
- 高性能Bloom filter索引

### 4. matrix - 生成矩阵 | Generate Matrix

**功能**: 生成gene×sample存在/缺失矩阵

**用法**:
```bash
biopytools kmertools matrix [OPTIONS]
```

**选项**:

| 参数 | 简写 | 默认值 | 描述 |
|------|------|--------|------|
| `--query-result` | `-i` | 必需 | 查询结果文件 |
| `--pos-file` | `-p` | 必需 | kmer位置文件 |
| `--matrix-output` | `-o` | 必需 | 矩阵输出文件 |

**输出文件**:
- `*.matrix.txt`: gene×sample矩阵
- `*.matrix.heatmap.txt`: 转置后的矩阵（用于热图）

**处理步骤**:
1. 添加位置信息到查询结果
2. 生成gene_pos格式的矩阵
3. 转置矩阵（样本为行，kmer为列）

### 5. screener - 筛选样本 | Screen Samples

**功能**: 基于连续0数量筛选样本

**用法**:
```bash
biopytools kmertools screener [OPTIONS]
```

**选项**:

| 参数 | 简写 | 默认值 | 描述 |
|------|------|--------|------|
| `--gene-list` | `-g` | 必需 | 基因列表文件 |
| `--sample-list` | `-s` | 必需 | 样本列表文件 |
| `--matrix-dir` | `-d` | 必需 | 矩阵文件目录 |
| `--output-dir` | `-o` | 必需 | 输出目录 |
| `--max-zeros` | | 51 | 最大连续0数量 |
| `--matrix-suffix` | | .kmer.matrix.heatmap.cluster.txt | 矩阵文件后缀 |
| `--delimiter` | | \t | 矩阵文件分隔符 |

**筛选逻辑**:
- 全1样本：所有位点均为1的样本直接保留
- 连续0限制：连续0的长度不超过max_zeros的样本保留
- 用于识别可能的假阴性或低覆盖区域

**输出文件**:
- `screening_results.csv`: 筛选结果汇总（gene, sample_names, sample_count）
- `gene_sample_matrix.csv`: 样本-基因矩阵（sample为行，gene为列）

## 工作流程 | Workflow

完整的工作流程如下：

```
步骤1: 构建kmer数据库 (build)
  ├─> 生成FOF文件
  ├─> kmtricks pipeline (kmer计数)
  ├─> kmtricks aggregate (矩阵聚合)
  ├─> bgzip压缩矩阵
  └─> 导入RocksDB

步骤2: 提取kmer (extract)
  └─> 从基因FASTA提取kmer

步骤3: 查询kmer (query)
  └─> 从RocksDB查询kmer存在情况

步骤4: 生成矩阵 (matrix)
  ├─> 添加位置信息
  ├─> 生成gene×sample矩阵
  └─> 转置矩阵

步骤5: 筛选样本 (screener)
  ├─> 筛选符合条件的样本
  └─> 生成样本-基因矩阵
```

## 使用示例 | Usage Examples

### 示例1：完整流程分析 | Example 1: Complete Workflow

```bash
# 1. 构建kmer数据库
biopytools kmertools build \
    -i ./clean_fastq \
    -k 51 \
    -o ./kmertools_results

# 2. 提取基因kmer
biopytools kmertools extract \
    -i Pi63.fa \
    --kmer-output Pi63.kmer.txt \
    --kmer-pos-output Pi63.pos.txt

# 3. 查询kmer
biopytools kmertools query \
    -r ./kmertools_results/kmer_rocksdb \
    -i Pi63.kmer.txt \
    -o Pi63.query.txt

# 4. 生成矩阵
biopytools kmertools matrix \
    -i Pi63.query.txt \
    -p Pi63.pos.txt \
    -o Pi63.matrix.txt

# 5. 筛选样本
biopytools kmertools screener \
    -g gene_list.txt \
    -s sample_list.txt \
    -d ./gene_matrices \
    -o ./screening_results \
    --max-zeros 45
```

### 示例2：批量基因分析 | Example 2: Batch Gene Analysis

```bash
# 为多个基因创建分析脚本
for gene in $(cat gene_list.txt); do
    # 提取kmer
    biopytools kmertools extract \
        -i ${gene}.fa \
        --kmer-output ${gene}.kmer.txt \
        --kmer-pos-output ${gene}.pos.txt

    # 查询kmer
    biopytools kmertools query \
        -r ./kmer_rocksdb \
        -i ${gene}.kmer.txt \
        -o ${gene}.query.txt

    # 生成矩阵
    biopytools kmertools matrix \
        -i ${gene}.query.txt \
        -p ${gene}.pos.txt \
        -o ${gene}.matrix.txt
done
```

## 通用参数 | Common Parameters

所有子命令都支持的通用参数：

| 参数 | 简写 | 默认值 | 描述 |
|------|------|--------|------|
| `--threads` | `-t` | 12 | 线程数 |
| `--output-dir` | `-o` | ./kmertools_output | 输出目录 |

## 注意事项 | Important Notes

1. **依赖软件**: 需要安装kmtricks、python-rocksdb、pyfastx
2. **内存需求**: RocksDB导入需要较大内存，建议16GB以上
3. **线程数**: build步骤建议使用较多线程（如88），默认为12线程
4. **kmer大小**: 常用kmer大小为31或51
5. **样本命名**: FOF文件中样本名不能包含特殊符号

## 故障排除 | Troubleshooting

### 常见问题 | Common Issues

**Q: "kmtricks: command not found"**
```bash
# 安装kmtricks
conda install -c bioconda kmtricks
```

**Q: "No module named 'rocksdb'"**
```bash
# 安装python-rocksdb
pip install python-rocksdb
```

**Q: "No module named 'pyfastx'"**
```bash
# 安装pyfastx
pip install pyfastx
```

**Q: RocksDB导入失败**
- 检查内存是否充足
- 检查矩阵文件是否正确压缩
- 确保bgzip已安装

## 技术细节 | Technical Details

### Canonical Kmer

工具使用canonical kmer表示，即取kmer和其反向互补序列中字典序较小的一个作为标准表示。

### 反向互补查询

查询步骤会自动检测kmer的反向互补序列，确保线状和环状序列都能正确匹配。

### 样本筛选逻辑

- **全1样本**: 所有位点均为1的样本直接保留
- **连续0限制**: 连续0的长度不超过指定阈值的样本保留
- **用途**: 识别可能的假阴性或低覆盖区域

## 命令对比 | Command Comparison

### 原始脚本 vs KmerTools

| 步骤 | 原始脚本 | KmerTools |
|------|---------|-----------|
| 构建库 | kmtricks pipeline + aggregate + import_rocksdb_rice.py | `kmertools build` |
| 提取kmer | get_kmer_from_fasta.py | `kmertools extract` |
| 查询kmer | search_kmer_from_rocksdb_rc.py | `kmertools query` |
| 生成矩阵 | kmer.matrix.add.py + get_kmer_matrix.py | `kmertools matrix` |
| 筛选样本 | kmer_matrix_screener.py | `kmertools screener` |

## 参考资源 | References

- [kmtricks官方文档](https://github.com/tlemane/kmtricks)
- [RocksDB文档](https://github.com/facebook/rocksdb)
- [seqkit](https://github.com/shenwei356/seqkit) - 命令风格参考
- 原始脚本参考：马老师kmer脚本

## 版本历史 | Version History

| 版本 | 日期 | 主要变更 |
|------|------|----------|
| 1.0.0 | 2026-01-29 | 初始版本，采用seqkit风格的子命令设计 |

## 许可证 | License

本项目采用MIT许可证 - 详见 [LICENSE](LICENSE) 文件
