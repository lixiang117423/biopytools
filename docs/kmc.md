# KMC K-mer Analysis Tool | KMC K-mer分析工具

## 功能介绍 | Features

基于KMC (K-mer Counter) 的k-mer分析工具，支持：

- **k-mer统计** (count): 为每个样本统计k-mer丰度
- **丰度矩阵构建** (matrix): 构建跨样本的k-mer丰度矩阵
- **k-mer查询** (query): 查询特定k-mer在各样本中的丰度
- **增量更新** (add): 添加新样本到现有矩阵

K-mer analysis tool based on KMC (K-mer Counter), supporting:

- **k-mer counting** (count): Count k-mer abundance for each sample
- **abundance matrix building** (matrix): Build cross-sample k-mer abundance matrix
- **k-mer query** (query): Query specific k-mer abundance across samples
- **incremental update** (add): Add new samples to existing matrix

## 安装 | Installation

确保KMC已安装在指定路径：

```bash
# KMC默认路径
/share/org/YZWL/yzwl_lixg/miniforge3/envs/kmc_v.3.2.4/bin
```

Ensure KMC is installed at the specified path:

```bash
# Default KMC path
/share/org/YZWL/yzwl_lixg/miniforge3/envs/kmc_v.3.2.4/bin
```

## 使用方法 | Usage

### 完整工作流程 | Complete Workflow

```bash
# 步骤1: 统计k-mer（建立数据库）
biopytools kmc count -i sample1.fq -i sample2.fq -i sample3.fq -k 21 -o kmc_output

# 步骤2: 构建丰度矩阵
biopytools kmc matrix -o kmc_output

# 步骤3: 查询k-mer
biopytools kmc query -q AAAAAAAAAAAAAAAAAAAAAAA -o kmc_output

# 步骤4: 添加新样本（可选）
biopytools kmc add -i new_sample.fq -n NewSample -o kmc_output

# 步骤5: 再次查询
biopytools kmc query -q AAAAAAAAAAAAAAAAAAAAAAA -o kmc_output
```

### 1. k-mer统计 | k-mer Counting

统计单个或多个样本的k-mer：

```bash
# 单个样本 | Single sample
biopytools kmc count -i sample1.fq -k 21 -o kmc_output

# 多个样本 | Multiple samples
biopytools kmc count -i sample1.fq -i sample2.fq -i sample3.fq -k 21 -o kmc_output

# 自定义样本名 | Custom sample names
biopytools kmc count -i s1.fq -i s2.fq -n Sample1 -n Sample2 -k 21 -o kmc_output
```

Count k-mers for single or multiple samples:

```bash
# Single sample
biopytools kmc count -i sample1.fq -k 21 -o kmc_output

# Multiple samples
biopytools kmc count -i sample1.fq -i sample2.fq -i sample3.fq -k 21 -o kmc_output

# Custom sample names
biopytools kmc count -i s1.fq -i s2.fq -n Sample1 -n Sample2 -k 21 -o kmc_output
```

### 2. 丰度矩阵构建 | Abundance Matrix Building

**注意：需要先运行 `count` 命令建立数据库**

基于已存在的 KMC 数据库自动构建跨样本丰度矩阵：

**方式1：不指定输入目录（使用-o作为输入和输出）|Method 1: Don't specify input directory (use -o for both input and output)**

```bash
# 会自动扫描 kmc_output/kmc_databases/ 目录
biopytools kmc matrix -o kmc_output
```

**方式2：指定输入目录（推荐）|Method 2: Specify input directory (recommended)**

```bash
# -i 指定包含kmc_databases的目录（即count步骤的-o参数）
biopytools kmc matrix -i kmc_database -o kmc_database
```

Build cross-sample k-mer abundance matrix from existing KMC databases:

**Note: Run `count` command first to build databases**

**Method 1: Don't specify input directory (use -o for both input and output)**

```bash
# Auto-scan kmc_output/kmc_databases/ directory
biopytools kmc matrix -o kmc_output
```

**Method 2: Specify input directory (recommended)**

```bash
# -i specifies directory containing kmc_databases (i.e., -o from count step)
biopytools kmc matrix -i kmc_database -o kmc_database
```

**工作原理|How it works:**
- 如果指定 `-i/--input-dir`：扫描 `input_dir/kmc_databases/` 目录|If `-i/--input-dir` specified: scan `input_dir/kmc_databases/` directory
- 如果未指定 `-i`：扫描 `output_dir/kmc_databases/` 目录|If `-i` not specified: scan `output_dir/kmc_databases/` directory
- 查找所有 `.kmc_pre` 和 `.kmc_suf` 文件|Find all `.kmc_pre` and `.kmc_suf` files
- 提取样本名称并构建矩阵|Extract sample names and build matrix

**输出文件 | Output files:**

- `abundance_matrix.h5`: 丰度矩阵 | Abundance matrix
- `kmer_dictionary.h5`: k-mer字典 | K-mer dictionary

### 3. k-mer查询 | k-mer Query

查询特定k-mer在各样本中的丰度：

```bash
biopytools kmc query -q AAAAAAAAAAAAAAAAAAAAAAAAA -o kmc_output
```

Query specific k-mer abundance across samples:

```bash
biopytools kmc query -q AAAAAAAAAAAAAAAAAAAAAAAAA -o kmc_output
```

### 4. 增量更新 | Incremental Update

**方式1：添加新样本（推荐）|Method 1: Add new samples (recommended)**

自动为新样本建库并更新矩阵：

```bash
biopytools kmc add -i new_sample1.fq -i new_sample2.fq -n New1 -n New2 -o kmc_output
```

**方式2：仅更新矩阵|Method 2: Update matrix only**

当新样本的数据库已存在时，仅更新矩阵：

```bash
biopytools kmc add -o kmc_output
```

Add new samples to existing matrix:

**Method 1: Add new samples (recommended)**

Automatically build databases and update matrix:

```bash
biopytools kmc add -i new_sample1.fq -i new_sample2.fq -n New1 -n New2 -o kmc_output
```

**Method 2: Update matrix only**

When new sample databases already exist:

```bash
biopytools kmc add -o kmc_output
```

## 参数说明 | Parameters

### 核心参数 | Core Parameters

| 参数 | Parameter | 默认值 | Default | 说明 | Description |
|------|-----------|--------|---------|------|-------------|
| `-m, --mode` | | count | count | 操作模式 | Operation mode |
| `-i, --input` | | 必需 | required | 输入文件(可多次使用) | Input files |
| `-n, --sample-names` | | - | - | 样本名称(可多次使用) | Sample names |
| `-k, --kmer-size` | | 21 | 21 | k-mer大小 | k-mer size |
| `--min-count` | | 2 | 2 | 最小计数阈值 | Minimum count |
| `--max-count` | | - | - | 最大计数阈值 | Maximum count |

### 路径参数 | Path Parameters

| 参数 | Parameter | 默认值 | Default | 说明 | Description |
|------|-----------|--------|---------|------|-------------|
| `-o, --output-dir` | | ./kmc_output | - | 输出目录 | Output directory |
| `--tmp-dir` | | ./kmc_tmp | - | 临时文件目录 | Temp directory |
| `--kmc-path` | | /share/.../bin | - | KMC软件路径 | KMC path |

### 处理参数 | Processing Parameters

| 参数 | Parameter | 默认值 | Default | 说明 | Description |
|------|-----------|--------|---------|------|-------------|
| `-t, --threads` | | 12 | 12 | 线程数 | Thread count |
| `--memory-limit` | | - | - | 内存限制(如: 12G) | Memory limit |

### 矩阵参数 | Matrix Parameters

| 参数 | Parameter | 默认值 | Default | 说明 | Description |
|------|-----------|--------|---------|------|-------------|
| `--matrix-format` | | hdf5 | hdf5 | 存储格式 | Storage format |
| `--dense-matrix` | | False | - | 使用密集矩阵 | Use dense matrix |

## Python API使用 | Python API Usage

```python
from biopytools.kmc import KMCManager, KMCConfig

# 方法1: 使用KMCManager | Method 1: Use KMCManager
manager = KMCManager(
    mode='count',
    input_files=['sample1.fq', 'sample2.fq'],
    kmer_size=21,
    output_dir='./kmc_output'
)
manager.run_analysis()

# 方法2: 直接使用各模块 | Method 2: Use modules directly
from biopytools.kmc import KMCCounter, KMCMatrixBuilder, KMCQuery

config = KMCConfig(
    input_files=['sample1.fq', 'sample2.fq'],
    kmer_size=21,
    output_dir='./kmc_output'
)

# 统计k-mer | Count k-mers
counter = KMCCounter(config)
counter.run()

# 构建矩阵 | Build matrix
matrix_builder = KMCMatrixBuilder(config)
matrix_builder.build_matrix(['sample1', 'sample2'])

# 查询k-mer | Query k-mer
query = KMCQuery(config)
abundances = query.query_kmer('AAAAAAAAAAAAAAAAAAAAAAA')
```

## 输出文件说明 | Output Files

### 目录结构 | Directory Structure

```
kmc_output/
├── kmc_databases/          # KMC数据库 | KMC databases
│   ├── sample1.kmc_pre
│   ├── sample1.kmc_suf
│   ├── sample2.kmc_pre
│   └── sample2.kmc_suf
├── abundance_matrix.h5     # 丰度矩阵 | Abundance matrix
├── kmer_dictionary.h5      # k-mer字典 | K-mer dictionary
├── global_kmers.txt        # 全局k-mer列表 | Global k-mer list
└── kmc_analysis.log        # 日志文件 | Log file
```

### HDF5矩阵结构 | HDF5 Matrix Structure

**abundance_matrix.h5:**

```python
# 密集矩阵 | Dense matrix
abundance          # shape: (n_kmers, n_samples)

# 稀疏矩阵 | Sparse matrix (默认 | default)
kmer_id            # k-mer ID数组
sample_id          # 样本ID数组
abundance          # 丰度数组
```

**kmer_dictionary.h5:**

```python
kmer               # k-mer序列数组
```

## 性能优化建议 | Performance Optimization

### 1. k-mer长度选择 | k-mer Size Selection

- **21-31**: 短 reads，基因组组装 | Short reads, genome assembly
- **51-71**: 长 reads，变异检测 | Long reads, variant detection
- **推荐值 | Recommended**: 21 (通用 | general purpose)

### 2. 内存和线程 | Memory and Threads

```bash
# 高内存服务器 | High-memory server
biopytools kmc -m count -i *.fq -k 21 -t 24 --memory-limit 64G

# 低内存环境 | Low-memory environment
biopytools kmc -m count -i *.fq -k 21 -t 4 --memory-limit 8G
```

### 3. 稀疏vs密集矩阵 | Sparse vs Dense Matrix

- **稀疏 (默认)**: 适合k-mer稀疏的数据 | For sparse k-mer data
- **密集**: 适合k-mer密集的小数据集 | For dense small datasets

## 常见问题 | FAQ

### Q1: 如何选择k-mer长度？

**A:** 取决于数据类型和目的：
- 测序错误检测：使用较短的k-mer (15-21)
- 基因组组装：使用中等长度k-mer (21-31)
- 物种鉴定：使用较长k-mer (31-51)

**A:** Depends on data type and purpose:
- Sequencing error detection: use shorter k-mers (15-21)
- Genome assembly: use medium k-mers (21-31)
- Species identification: use longer k-mers (31-51)

### Q2: 如何处理大量样本？

**A:** 对于 >1000 样本：
1. 分批处理 (batch_size=100)
2. 使用稀疏矩阵 (默认)
3. 确保足够的临时存储空间

**A:** For >1000 samples:
1. Process in batches (batch_size=100)
2. Use sparse matrix (default)
3. Ensure sufficient temp storage

### Q3: 如何导出矩阵为TSV？

**A:** 使用Python API:

```python
from biopytools.kmc import KMCQuery

query = KMCQuery(config)
query.export_matrix_to_tsv('output.tsv')
```

## 参考资源 | Resources

- [KMC官方文档](https://github.com/refresh-bio/KMC)
- [KMC API文档](https://github.com/refresh-bio/KMC/wiki)
- [biopytools文档](../README.md)

## 版本历史 | Version History

| 版本 | Version | 日期 | Date | 说明 | Notes |
|------|----------|------|------|------|-------|
| 1.0.0 | 1.0.0 | 2025-01-22 | 初始版本 | Initial release |
