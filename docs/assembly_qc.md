# 基因组组装质量评估工具 (assembly-qc)

## 功能描述

基因组组装质量评估工具提供基因组组装质量的综合评估，包括：

- **核心评估**（基于基因组本身，不依赖外部数据）
  - BUSCO完整性评估
  - LAI指数评估（LTR组装完整性）

- **可选评估**（基于外部数据）
  - QV质量值计算（基于k-mer分析）
  - Mapping评估（基于NGS/TGS数据比对）

- **输出报告**
  - HTML综合质量报告
  - 发表用表格（TSV/Excel格式）

## 使用方法

### 基本用法

```bash
# 核心评估（只基于基因组，不需要额外数据）
biopytools assembly-qc \
    --genome genome.fa \
    --lineage embryophyta_odb10 \
    -o qc_results
```

### 添加QV评估

```bash
biopytools assembly-qc \
    --genome genome.fa \
    --lineage embryophyta_odb10 \
    --qv-reads ./illumina_reads \
    -o qc_results
```

### 添加Mapping评估

```bash
biopytools assembly-qc \
    --genome genome.fa \
    --lineage embryophyta_odb10 \
    --mapping-reads ./ngs_reads \
    -o qc_results
```

### 完整评估（复用reads数据）

```bash
biopytools assembly-qc \
    --genome genome.fa \
    --lineage embryophyta_odb10 \
    --reads-for-all ./illumina_reads \
    --enable-qv \
    --enable-mapping \
    -o qc_results
```

## 参数说明

### 必需参数

| 参数 | 简写 | 说明 |
|------|------|------|
| `--genome` | `-g` | 基因组FASTA文件 |
| `--lineage` | `-l` | BUSCO数据集名称（如embryophyta_odb10）或完整路径 |
| `--output-dir` | `-o` | 输出目录（默认：./assembly_qc_output） |

### 可选参数

#### 样品信息

| 参数 | 简写 | 默认值 | 说明 |
|------|------|--------|------|
| `--sample-name` | `-s` | genome_sample | 样品名称 |

#### 核心评估控制

| 参数 | 说明 |
|------|------|
| `--skip-busco` | 跳过BUSCO评估 |
| `--skip-lai` | 跳过LAI评估 |

#### QV评估参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--enable-qv` | False | 启用QV评估 |
| `--qv-reads` | None | 用于QV的reads目录 |
| `--qv-kmer-size` | None | k-mer大小（None表示自动选择） |
| `--qv-threads` | 24 | QV计算线程数 |

#### Mapping评估参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--enable-mapping` | False | 启用Mapping评估 |
| `--mapping-reads` | None | 用于Mapping的reads目录 |
| `--mapping-type` | ngs | Mapping类型（ngs/pacbio/ont） |
| `--mapping-pattern` | _1.clean.fq.gz | FASTQ文件匹配模式 |
| `--mapping-threads` | 64 | Mapping线程数 |

#### 数据复用

| 参数 | 说明 |
|------|------|
| `--reads-for-all` | 同时用于QV和Mapping的reads目录 |

#### 报告参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--no-html` | False | 不生成HTML报告 |
| `--no-table` | False | 不生成表格 |
| `--table-format` | both | 表格格式（tsv/xlsx/both） |

#### 线程参数

| 参数 | 简写 | 默认值 | 说明 |
|------|------|--------|------|
| `--threads` | `-t` | 64 | 全局线程数（所有评估步骤都使用全部线程）|

**线程使用规则**：
- 由于工作流是**串行执行**的，每个评估步骤都使用全部线程以最大化性能
- BUSCO评估：64线程
- LAI评估：64线程
- QV评估：64线程
- Mapping评估：64线程

## 输出文件

### 目录结构

```
qc_results/
├── 01.busco_evaluation/        # BUSCO评估输出
├── 02.lai_evaluation/          # LAI评估输出
├── 03.qv_evaluation/           # QV评估输出
├── 04.mapping_evaluation/      # Mapping评估输出
├── assembly_qc_report.html     # 综合质量报告
├── assembly_qc_table.tsv       # 发表用表格（TSV）
└── assembly_qc_table.xlsx      # 发表用表格（Excel）
```

### 表格格式

输出表格包含以下列：

| 列名 | 说明 |
|------|------|
| Sample | 样品名称 |
| Assembly_Size_(Mb) | 组装大小 |
| Contig_N50_(Mb) | Contig N50 |
| Scaffold_N50_(Mb) | Scaffold N50 |
| GC_% | GC含量 |
| BUSCO_Complete_% | BUSCO完整度 |
| BUSCO_Single_% | 单拷贝完整比例 |
| BUSCO_Duplicated_% | 多拷贝完整比例 |
| BUSCO_Fragmented_% | 碎片化比例 |
| BUSCO_Missing_% | 缺失比例 |
| LAI | LAI指数 |
| QV | QV质量值 |
| Mapped_Reads_% | 比对率 |
| Mean_Coverage_(X) | 平均覆盖度 |
| Coverage_≥10X_% | ≥10X覆盖比例 |

## 使用示例

### 示例1：植物基因组评估

```bash
biopytools assembly-qc \
    --genome plant_genome.fa \
    --lineage embryophyta_odb10 \
    --sample-name Arabidopsis_thaliana \
    -o plant_qc_results
```

### 示例2：完整质量评估（包含QV和Mapping）

```bash
biopytools assembly-qc \
    --genome genome.fa \
    --lineage embryophyta_odb10 \
    --sample-name Sample01 \
    --reads-for-all ./illumina_reads \
    --enable-qv \
    --enable-mapping \
    --threads 128 \
    -o qc_results
```

### 示例3：仅Mapping评估（多样品）

```bash
biopytools assembly-qc \
    --genome genome.fa \
    --lineage embryophyta_odb10 \
    --mapping-reads ./ngs_data \
    --mapping-type ngs \
    --mapping-pattern "_1.clean.fq.gz" \
    -o mapping_qc_results
```

## 注意事项

1. **BUSCO数据集路径**：默认使用 `/share/org/YZWL/yzwl_lixg/database/busco`，如需修改请设置环境变量 `BUSCO_DATASET_PATH`

2. **Mapping评估支持多样品**：会在指定目录下自动发现所有符合命名规则的FASTQ文件对（`_1.clean.fq.gz` / `_2.clean.fq.gz`）

3. **线程数设置**：
   - 使用统一的 `--threads` 参数
   - 由于工作流是**串行执行**的，每个评估步骤都使用全部线程以最大化性能

4. **断点续传**：默认启用，已完成的步骤会自动跳过

## 软件依赖

### 核心依赖

- BUSCO (busco)
- LAI工具包
  - LTR_harvest (gt)
  - LTR_FINDER_parallel
  - LTR_retriever
  - LAI

### 可选依赖

- Merqury (QV评估)
- BWA (Mapping评估)
- Samtools (Mapping评估)
- Bedtools (Mapping评估)

## 常见问题

### Q1: BUSCO评估失败

**原因**：BUSCO数据集路径不正确或谱系名称错误

**解决**：
```bash
# 检查BUSCO数据集
ls /share/org/YZWL/yzwl_lixg/database/busco

# 使用正确的谱系名称
biopytools busco --list-datasets  # 查看可用谱系
```

### Q2: Mapping评估未发现样品

**原因**：FASTQ文件命名不符合模式

**解决**：确保文件命名匹配默认模式 `_1.clean.fq.gz`，或使用 `--mapping-pattern` 指定正确模式

```bash
# 查看实际文件名
ls ./ngs_data/

# 使用正确的模式
biopytools assembly-qc ... --mapping-pattern "_R1.fastq.gz"
```

### Q3: Excel表格生成失败

**原因**：缺少pandas或openpyxl包

**解决**：
```bash
pip install pandas openpyxl
```

或只生成TSV格式：
```bash
biopytools assembly-qc ... --table-format tsv
```

## 版本信息

- 当前版本：1.0.0
- 最后更新：2026-03-11
