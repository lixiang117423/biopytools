# RNA-Bloom 转录组从头组装模块

**专业的无参考转录组组装工具 | Professional De Novo Transcriptome Assembly Tool**

## 功能概述 | Overview

RNA-Bloom转录组从头组装模块是一个基于RNA-Bloom软件的强大转录组组装工具，专门用于无参考基因组情况下的转录序列从头组装。该模块支持多种测序数据类型，包括短reads、长reads以及单细胞RNA-seq数据，具有内存效率高、组装速度快、适用范围广等特点，适用于非模式生物、转录组研究和新基因发现等多种应用场景。

## 主要特性 | Key Features

- **多数据类型支持**：支持bulk RNA-seq、单细胞RNA-seq、长reads测序数据
- **高效内存管理**：基于Bloom filter的k-mer计数，内存占用小
- **灵活输入模式**：支持paired-end、single-end、混合数据组装
- **链特异性组装**：支持strand-specific数据组装
- **单细胞混合组装**：优化的单细胞pooled assembly模式
- **长reads支持**：支持ONT和PacBio长reads数据
- **参考引导组装**：可选择使用参考转录本引导组装
- **自动质量控制**：内置序列验证和依赖检查
- **详细日志记录**：完整的组装过程日志和错误追踪

## 快速开始 | Quick Start

### 基本用法 | Basic Usage

```bash
# 短reads双端组装（最常用）
biopytools rnabloom \
    --left reads_1.fq \
    --right reads_2.fq \
    -o ./assembly_results

# 短reads单端组装
biopytools rnabloom \
    --sef reads.fq \
    -o ./assembly_results

# 长reads组装
biopytools rnabloom \
    --long long_reads.fq \
    -o ./assembly_results
```

### 高级用法 | Advanced Usage

```bash
# 链特异性数据组装
biopytools rnabloom \
    --left reads_2.fq \
    --right reads_1.fq \
    --stranded \
    --revcomp-right \
    -t 24 \
    -o ./stranded_assembly

# 长reads + 短reads混合组装
biopytools rnabloom \
    --long long_reads.fq \
    --sef short_reads.fq \
    -t 24 \
    -o ./hybrid_assembly

# 单细胞混合组装
biopytools rnabloom \
    --cell-list cells.txt \
    -t 24 \
    -o ./single_cell_assembly

# PacBio长reads组装
biopytools rnabloom \
    --long pacbio_reads.fq \
    --pacbio \
    -t 24 \
    -o ./pacbio_assembly

# 参考引导组装
biopytools rnabloom \
    --left reads_1.fq \
    --right reads_2.fq \
    --ref reference_transcripts.fa \
    -o ./guided_assembly
```

## 参数说明 | Parameters

### 必需参数（至少指定一种）| Required Parameters (Specify at Least One)

| 参数 | 描述 | 示例 |
|------|------|------|
| `--left, -1` | 左端reads文件（paired-end）| `--left reads_1.fq` |
| `--right, -2` | 右端reads文件（paired-end）| `--right reads_2.fq` |
| `--sef` | 单端正向reads文件 | `--sef single_end.fq` |
| `--ser` | 单端反向reads文件 | `--ser single_reverse.fq` |
| `--long` | 长reads文件（ONT/PacBio）| `--long ont_reads.fq` |
| `--cell-list` | 单细胞列表文件 | `--cell-list cells.txt` |

**注意**：
- paired-end模式：必须同时指定`--left`和`--right`
- 单细胞模式：使用`--cell-list`，不能与其他输入参数混用

### 输出参数 | Output Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-o, --output-dir` | **必需** | 输出目录路径 |

### 处理参数 | Processing Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-t, --threads` | `12` | 线程数 |
| `--rnabloom-path` | `rnabloom` | RNA-Bloom工具路径 |

### Bloom Filter配置 | Bloom Filter Configuration

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--mem` | `自动` | Bloom filter总大小（GB）|
| `--fpr` | `自动` | 假阳性率（0-1之间）|
| `--nk, --num-kmers` | `自动` | 唯一k-mer数量 |

**说明**：如果不指定这些参数，RNA-Bloom会自动使用ntCard估算最佳值。

### 数据类型配置 | Data Type Configuration

| 参数 | 描述 | 适用场景 |
|------|------|----------|
| `--stranded` | 链特异性数据 | strand-specific RNA-seq |
| `--revcomp-left` | 反向互补左端reads | F2R1 orientation |
| `--revcomp-right` | 反向互补右端reads | F1R2 orientation |
| `--pacbio` | PacBio数据 | PacBio cDNA/Direct RNA |

### 参考引导组装 | Reference-Guided Assembly

| 参数 | 描述 |
|------|------|
| `--ref, --reference` | 参考转录本FASTA文件 |

### 输出选项 | Output Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--min-length` | `200` | 最小转录本长度（bp）|
| `--uracil` | `False` | 输出尿嘧啶(U)而非胸腺嘧啶(T)|
| `--no-nr` | `False` | 不导出去冗余转录本 |

### 处理控制 | Processing Control

| 参数 | 描述 | 阶段说明 |
|------|------|----------|
| `--stage` | 停止阶段（1-3）| 1: 构建图<br>2: 组装片段/纠正reads<br>3: 组装转录本 |

## 输入文件格式 | Input File Formats

### 短reads FASTQ格式 | Short-Read FASTQ Format

标准FASTQ格式的测序数据：
```
@READ_ID
GATCGATCGATCGATCGATCG
+
IIIIIIIIIIIIIIIIIIIIIII
```

**要求**：
- 可以是gzip压缩格式（.fq.gz, .fastq.gz）
- 支持多文件输入（用空格分隔）

### 长reads FASTQ格式 | Long-Read FASTQ Format

ONT或PacBio长reads数据：
```
@read_id
GATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
+
####################################
```

**注意**：
- 长reads不能使用纯数字ID（如1, 2, 3）
- 如需重命名，使用`seqtk rename`工具

### 单细胞列表文件 | Single-Cell List File

单细胞pooled assembly的输入文件格式：

```
#name left right
cell1 /path/to/cell1/left.fq /path/to/cell1/right.fq
cell2 /path/to/cell2/left.fq /path/to/cell2/right.fq
cell3 /path/to/cell3/left.fq /path/to/cell3/right.fq
```

**说明**：
- 第一行为头部，以`#`开头
- 列之间用空格或Tab分隔
- 同一个cell可以有多行（用于多个read文件）
- 也支持单端reads（添加`sef`和`ser`列）

### 参考转录本文件 | Reference Transcript File

标准FASTA格式的参考转录本序列：
```
>transcript1
ATGGCGATCGATCGATCGATCGATCGATCGATCGATCG
>transcript2
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
```

## 输出结果 | Output Results

### 输出文件 | Output Files

| 文件名 | 描述 |
|--------|------|
| `rnabloom.transcripts.fa` | 主要转录本（≥min_length）|
| `rnabloom.transcripts.short.fa` | 短转录本（<min_length）|
| `rnabloom.transcripts.nr.fa` | 去冗余转录本 |
| `rnabloom_assembly.log` | 完整日志文件 |

### 结果解读 | Result Interpretation

**主要转录本文件**：
- 包含所有长度≥阈值的组装转录本
- 推荐用于下游分析

**去冗余转录本**：
- 经过冗余去除的转录本集合
- 适合用于转录本注释和功能分析

## 使用示例 | Usage Examples

### 示例1：植物转录组组装 | Example 1: Plant Transcriptome Assembly

```bash
# 非模式植物转录组组装
biopytools rnabloom \
    --left plant_R1.fastq.gz \
    --right plant_R2.fastq.gz \
    --revcomp-right \
    -t 24 \
    -o ./plant_transcriptome
```

### 示例2：单细胞转录组组装 | Example 2: Single-Cell Transcriptome Assembly

```bash
# 准备单细胞列表文件 cells.txt
cat > cells.txt << EOF
#name left right
cell1 /data/cell1_R1.fq.gz /data/cell1_R2.fq.gz
cell2 /data/cell2_R1.fq.gz /data/cell2_R2.fq.gz
cell3 /data/cell3_R1.fq.gz /data/cell3_R2.fq.gz
EOF

# 运行单细胞混合组装
biopytools rnabloom \
    --cell-list cells.txt \
    -t 24 \
    -o ./single_cell_assembly
```

### 示例3：链特异性转录组组装 | Example 3: Stranded Transcriptome Assembly

```bash
# F2R1 orientation的链特异性数据
biopytools rnabloom \
    --left reads_2.fq \
    --right reads_1.fq \
    --stranded \
    --revcomp-right \
    -t 24 \
    -o ./stranded_transcriptome
```

### 示例4：ONT长reads组装 | Example 4: ONT Long-Read Assembly

```bash
# Oxford Nanopore长reads组装
biopytools rnabloom \
    --long ont_reads.fastq.gz \
    -t 24 \
    -o ./ont_assembly
```

### 示例5：PacBio长reads组装 | Example 5: PacBio Long-Read Assembly

```bash
# PacBio cDNA reads组装
biopytools rnabloom \
    --long pacbio_reads.fastq.gz \
    --pacbio \
    -t 24 \
    -o ./pacbio_assembly
```

### 示例6：混合组装（长reads + 短reads）| Example 6: Hybrid Assembly

```bash
# 使用短readspolish长reads组装
biopytools rnabloom \
    --long long_reads.fq \
    --sef short_reads.fq \
    -t 24 \
    -o ./hybrid_assembly
```

### 示例7：参考引导组装 | Example 7: Reference-Guided Assembly

```bash
# 使用近缘物种转录本引导
biopytools rnabloom \
    --left reads_R1.fq.gz \
    --right reads_R2.fq.gz \
    --ref related_species_transcripts.fa \
    -t 24 \
    -o ./guided_assembly
```

### 示例8：自定义Bloom Filter大小 | Example 8: Custom Bloom Filter Size

```bash
# 设置Bloom filter大小和假阳性率
biopytools rnabloom \
    --left reads_R1.fq.gz \
    --right reads_R2.fq.gz \
    --mem 16 \
    --fpr 0.01 \
    -t 24 \
    -o ./custom_memory_assembly
```

## 系统要求 | System Requirements

### 依赖软件 | Dependencies

**必需软件**：
- **Java** (版本 11 或 17)
  - 下载：https://www.oracle.com/java/technologies/downloads/
- **minimap2** (≥2.22)
  - 安装：`conda install -c bioconda minimap2`

**可选软件**：
- **ntCard** (≥1.2.1)
  - 安装：`conda install -c bioconda ntcard`
  - 用于自动估算Bloom filter大小
- **Racon** (仅长readspolish需要)
  - 安装：`conda install -c bioconda racon`

### 安装RNA-Bloom | Installing RNA-Bloom

**方法1：使用conda（推荐）**
```bash
# 使用conda安装
conda install -c bioconda rnabloom

# 或使用mamba
mamba install -c bioconda rnabloom
```

**方法2：从GitHub下载**
```bash
# 下载最新发布版本
wget https://github.com/bcgsc/RNA-Bloom/releases/download/v2.0.1/rnabloom_v2.0.1.tar.gz

# 解压
tar -xzf rnabloom_v2.0.1.tar.gz

# 使用
java -jar RNA-Bloom.jar ...
```

### 硬件建议 | Hardware Recommendations

| 数据类型 | CPU | 内存 | 存储 |
|----------|-----|------|------|
| 短reads | 8-16核 | 8-32 GB | 输入文件大小的5倍 |
| 长reads | 8-16核 | 16-64 GB | 输入文件大小的10倍 |
| 单细胞 | 16-32核 | 32-128 GB | 输入文件大小的10倍 |

## 注意事项 | Important Notes

### 数据质量控制 | Data Quality Control

1. **Adapter修剪**：
   - 长reads数据建议先用Porechop修剪adapters
   - 短reads建议用fastp进行质量控制

2. **Read ID重命名**：
   - 长reads不能使用纯数字ID
   - 使用`seqtk rename`重命名

3. **文件格式**：
   - 支持FASTQ和FASTA格式
   - 支持gzip压缩格式

### 参数选择建议 | Parameter Selection Guidelines

**Bloom Filter大小**：
- 小数据集（<10 GB reads）：使用默认值
- 大数据集（>10 GB reads）：考虑设置`--mem`
- 低假阳性率需求：设置`--fpr 0.01`

**链特异性数据**：
- F2R1 orientation：使用`--stranded --revcomp-right`
- F1R2 orientation：使用`--stranded --revcomp-left`
- 非链特异性：不需要额外参数

**长reads类型**：
- ONT数据：默认设置
- PacBio数据：添加`--pacbio`参数
- Direct RNA：添加`--stranded`参数

### 常见问题 | Common Issues

**Q: 内存不足错误**
```bash
# 减小Bloom filter大小
biopytools rnabloom ... --mem 8

# 或设置假阳性率
biopytools rnabloom ... --fpr 0.05
```

**Q: 未组装出转录本**
```bash
# 检查输入文件质量
# 检查reads是否足够
# 降低最小长度阈值
biopytools rnabloom ... --min-length 150
```

**Q: 速度慢**
```bash
# 增加线程数
biopytools rnabloom ... -t 32
```

## 故障排除 | Troubleshooting

### 依赖问题 | Dependency Issues

**Java未找到**
```bash
# 检查Java版本
java -version

# 安装Java 11或17
conda install openjdk=11
```

**minimap2未找到**
```bash
# 安装minimap2
conda install -c bioconda minimap2

# 验证安装
minimap2 --version
```

### 输入问题 | Input Issues

**文件格式错误**
```bash
# 验证FASTQ格式
zcat test.fq.gz | head -4

# 检查reads ID
zcat long_reads.fq.gz | head -1
```

**单细胞列表格式错误**
```bash
# 正确格式示例
#name left right
cell1 /path/to/R1.fq /path/to/R2.fq
```

### 输出问题 | Output Issues

**转录本数量少**
```bash
# 检查输入数据质量
# 降低最小长度阈值
biopytools rnabloom ... --min-length 100

# 使用参考引导
biopytools rnabloom ... --ref reference.fa
```

## 与其他模块配合 | Integration with Other Modules

### 转录组分析流程 | Transcriptome Analysis Pipeline

```bash
# 1. 质量控制
biopytools fastp -1 raw_R1.fq -2 raw_R2.fq -o clean_R1.fq -O clean_R2.fq

# 2. 转录组组装
biopytools rnabloom --left clean_R1.fq --right clean_R2.fq -o ./assembly

# 3. 转录本注释
biopytools interproscan -i assembly/rnabloom.transcripts.fa -o ./annotation

# 4. 功能分析
biopytools blast -q assembly/rnabloom.transcripts.fa -d nr_database -o blast_results
```

### 与dual_rnaseq的对比 | Comparison with dual_rnaseq

| 模块 | 适用场景 | 参考基因组需求 |
|------|----------|----------------|
| **rnabloom** | 无参考转录组组装 | 不需要 |
| **dual_rnaseq** | 双物种RNA-seq定量 | 需要 |

### 下游分析建议 | Downstream Analysis Recommendations

组装后的转录本可以用于：
- 转录本功能注释（InterProScan）
- 同源搜索（BLAST）
- ORF预测和蛋白序列提取
- 表达量定量（Salmon/Kallisto）
- 差异表达分析

## 引用信息 | Citation

如果在学术研究中使用RNA-Bloom，请引用：

**Long-read RNA-seq assembly**:
```
Ka Ming Nip, Saber Hafezqorani, Kristina K. Gagalova, Readman Chiu,
Chen Yang, René L. Warren, and Inanc Birol. Reference-free assembly
of long-read transcriptome sequencing data with RNA-Bloom2.
Nature Communications. 2023 May 22;14(1):2940. doi: 10.1038/s41467-023-38553-y
```

**Short-read RNA-seq assembly**:
```
Ka Ming Nip, Readman Chiu, Chen Yang, Justin Chu, Hamid Mohamadi,
René L. Warren, and Inanc Birol. RNA-Bloom enables reference-free
and reference-guided sequence assembly for single-cell transcriptomes.
Genome Research. 2020 Aug;30(8):1191-1200.
doi: 10.1101/gr.260174.119
```

## 许可证 | License

本项目采用MIT许可证 - 详见 [LICENSE](LICENSE) 文件

**注意**：RNA-Bloom软件本身遵循特定的许可证条款。

## 相关资源 | Related Resources

- [RNA-Bloom GitHub](https://github.com/bcgsc/RNA-Bloom)
- [RNA-Bloom Documentation](https://bcgsc.github.io/RNA-Bloom/)
- [Bioconda Recipe](https://anaconda.org/bioconda/rnabloom)
- [minimap2 GitHub](https://github.com/lh3/minimap2)
- [ntCard GitHub](https://github.com/bcgsc/ntCard)

---

**版本信息**: RNA-Bloom模块版本 1.0.0 | Module Version 1.0.0
