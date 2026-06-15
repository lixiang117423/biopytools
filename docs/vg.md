# VG变异图分析工具 | VG Variation Graph Analysis Tool

版本 | Version: 1.0.0
作者 | Author: Xiang LI
日期 | Date: 2026-04-02

## 概述 | Overview

VG工具是基于**VG (Variation Graph)** 工具包的Python封装，使用conda环境调用VG，提供完整的变异图分析流程。

The VG tool is a Python wrapper for **VG (Variation Graph)** toolkit, using conda environment to invoke VG, providing a complete variation graph analysis pipeline.

## 功能特点 | Features

- 🧬 **变异图构建**: 从VCF和参考基因组构建变异图|Build variation graph from VCF and reference
- 📊 **多格式索引**: 支持XG、GCSA、GBWT、GIRAFFE等多种索引格式|Multiple index formats: XG, GCSA, GBWT, GIRAFFE
- 🎯 **快速比对**: 使用Giraffe进行快速的序列比对|Fast read alignment with Giraffe
- 💾 **VCF导出**: 从变异图导出VCF文件|Export VCF from variation graph
- ⚡ **高性能**: 支持多线程并行计算|Multi-threaded parallel computing
- 🔧 **Conda调用**: 使用conda run调用VG|Use conda run to invoke VG

## 分析流程 | Analysis Pipeline

### 1. Construct - 变异图构建 | Graph Construction

从参考基因组和VCF文件构建变异图。

Build variation graph from reference genome and VCF file.

**输入文件 | Input Files:**
- 参考基因组FASTA文件 | Reference FASTA file
- VCF文件（建议压缩并索引）| VCF file (compressed and indexed recommended)

**输出文件 | Output Files:**
- `{output}.vg` - VG格式变异图 | VG format variation graph

### 2. Index - 索引创建 | Index Creation

为变异图创建索引，用于后续比对。

Create indexes for variation graph for subsequent alignment.

**索引类型 | Index Types:**
- `XG` - 扩展的图索引 | Extended graph index
- `GCSA` - 通用压缩后缀数组 | Generalized compressed suffix array
- `GBWT` - 图Burrows-Wheeler变换 | Graph Burrows-Wheeler transform
- `GIRAFFE` - Giraffe比对索引（包含.min和.dist）| Giraffe alignment indexes (.min and .dist)

### 3. Giraffe - 序列比对 | Read Alignment

使用Giraffe进行快速的序列比对。

Use Giraffe for fast read alignment.

**输入文件 | Input Files:**
- 索引的图文件（前缀）| Indexed graph file (prefix)
- FASTQ格式的reads | FASTQ format reads

**输出文件 | Output Files:**
- `{output}.gam` - GAM格式比对结果 | GAM format alignments
- `{output}.gaf` - GAF格式比对结果 | GAF format alignments

### 4. Deconstruct - VCF导出 | VCF Export

从变异图导出VCF文件。

Export VCF file from variation graph.

**输入文件 | Input Files:**
- VG格式图文件 | VG format graph file

**输出文件 | Output Files:**
- `{output}.vcf` - VCF格式变异文件 | VCF format variant file

## 安装和使用 | Installation and Usage

### 前置要求 | Prerequisites

#### VG安装 | VG Installation

```bash
# 创建conda环境并安装VG|Create conda environment and install VG
conda create -n vg_v.1.7.0 -c conda-forge -c bioconda vg
conda activate vg_v.1.7.0

# 验证安装|Verify installation
vg --version
```

#### Python环境要求 | Python Environment Requirements

```bash
# biopytools会自动处理Python依赖
# 确保biopytools已安装
pip install biopytools
```

### 基本用法 | Basic Usage

#### Construct - 构建变异图 | Build Variation Graph

```bash
# 基本用法|Basic usage
biopytools vg construct -r ref.fa -v variants.vcf.gz -o graph.vg

# 指定染色体区域|Specify chromosome region
biopytools vg construct -r ref.fa -v variants.vcf.gz -o graph.vg -R chr1

# 保存alt等位基因路径|Save alt allele paths
biopytools vg construct -r ref.fa -v variants.vcf.gz -o graph.vg --alt-paths

# 增加线程数|Increase threads
biopytools vg construct -r ref.fa -v variants.vcf.gz -o graph.vg -t 24
```

#### Index - 创建索引 | Create Indexes

```bash
# 创建XG索引|Create XG index
biopytools vg index -i graph.vg -o graph --xg

# 创建GIRAFFE索引（用于比对）|Create GIRAFFE indexes (for alignment)
biopytools vg index -i graph.vg -o graph --giraffe

# 创建所有索引|Create all indexes
biopytools vg index -i graph.vg -o graph --xg --gcsa --gbwt --giraffe

# 自定义k-mer大小|Custom k-mer size
biopytools vg index -i graph.vg -o graph --gcsa -k 31
```

#### Giraffe - 序列比对 | Read Alignment

```bash
# 单端测序比对|Single-end alignment
biopytools vg giraffe -g graph -f reads.fq -o alignments.gam

# 双端测序比对|Paired-end alignment
biopytools vg giraffe -g graph -f reads1.fq -f2 reads2.fq -o alignments.gam

# 输出GAF格式|Output GAF format
biopytools vg giraffe -g graph -f reads.fq -o alignments.gaf --format GAF

# 指定片段长度|Specify fragment length
biopytools vg giraffe -g graph -f reads.fq -o alignments.gam -l 500 -s 50

# 增加线程数|Increase threads
biopytools vg giraffe -g graph -f reads.fq -o alignments.gam -t 24
```

#### Deconstruct - 导出VCF | Export VCF

```bash
# 基本用法|Basic usage
biopytools vg deconstruct -i graph.vg -r ref_path -o output.vcf

# 指定样本|Specify samples
biopytools vg deconstruct -i graph.vg -r ref_path -o output.vcf -s sample1 -s sample2

# 增加线程数|Increase threads
biopytools vg deconstruct -i graph.vg -r ref_path -o output.vcf -t 24
```

### 参数说明 | Parameter Description

#### 全局参数 | Global Parameters

| 参数 | 默认值 | 说明 | 示例 |
|------|--------|------|------|
| `--vg-env` | `vg_v.1.7.0` | VG conda环境名称|VG conda environment name | `--vg-env my_vg` |
| `--log-level` | `INFO` | 日志级别|Log level | `--log-level DEBUG` |

#### Construct参数 | Construct Parameters

| 参数 | 默认值 | 说明 | 示例 |
|------|--------|------|------|
| `-r, --reference` | 必需 | 参考基因组FASTA文件|Reference FASTA file | `-r ref.fa` |
| `-v, --vcf` | 必需 | VCF文件|VCF file | `-v variants.vcf.gz` |
| `-o, --output` | 必需 | 输出VG文件|Output VG file | `-o graph.vg` |
| `-R, --region` | None | 指定染色体区域|Specify chromosome region | `-R chr1` |
| `-t, --threads` | `12` | 线程数|Number of threads | `-t 24` |
| `--alt-paths` | `False` | 保存alt等位基因路径|Save alt allele paths | `--alt-paths` |
| `--no-progress` | `False` | 不显示进度|Do not show progress | `--no-progress` |

#### Index参数 | Index Parameters

| 参数 | 默认值 | 说明 | 示例 |
|------|--------|------|------|
| `-i, --input` | 必需 | 输入图文件|Input graph file | `-i graph.vg` |
| `-o, --output` | 必需 | 输出前缀|Output prefix | `-o graph` |
| `--xg` | `False` | 创建XG索引|Create XG index | `--xg` |
| `--gcsa` | `False` | 创建GCSA索引|Create GCSA index | `--gcsa` |
| `--gbwt` | `False` | 创建GBWT索引|Create GBWT index | `--gbwt` |
| `--giraffe` | `False` | 创建GIRAFFE索引|Create GIRAFFE indexes | `--giraffe` |
| `-k, --kmer-size` | `16` | GCSA k-mer大小|GCSA k-mer size | `-k 31` |
| `-t, --threads` | `12` | 线程数|Number of threads | `-t 24` |

#### Giraffe参数 | Giraffe Parameters

| 参数 | 默认值 | 说明 | 示例 |
|------|--------|------|------|
| `-g, --graph` | 必需 | 图文件前缀（索引）|Graph file prefix (indexed) | `-g graph` |
| `-f, --reads` | 必需 | 输入reads文件|Input reads file | `-f reads.fq` |
| `-o, --output` | 必需 | 输出GAM文件|Output GAM file | `-o alignments.gam` |
| `-f2, --reads2` | None | 第二个reads文件|Second reads file | `-f2 reads2.fq` |
| `-t, --threads` | `12` | 线程数|Number of threads | `-t 24` |
| `-l, --fragment-length` | `0` | 片段长度（0=自动）|Fragment length (0=auto) | `-l 500` |
| `-s, --fragment-std-dev` | `0` | 片段长度标准差（0=自动）|Fragment length std dev (0=auto) | `-s 50` |
| `--min-identity` | `0.0` | 最小相似度|Min identity | `--min-identity 0.95` |
| `--format` | `GAM` | 输出格式|Output format | `--format GAF` |

#### Deconstruct参数 | Deconstruct Parameters

| 参数 | 默认值 | 说明 | 示例 |
|------|--------|------|------|
| `-i, --input` | 必需 | 输入图文件|Input graph file | `-i graph.vg` |
| `-o, --output` | 必需 | 输出VCF文件|Output VCF file | `-o output.vcf` |
| `-r, --reference-path` | 必需 | 参考路径名称|Reference path name | `-r ref_path` |
| `-s, --samples` | None | 样本列表（可多次使用）|Sample list (can be used multiple times) | `-s sample1 -s sample2` |
| `-t, --threads` | `12` | 线程数|Number of threads | `-t 24` |

## 输出文件 | Output Files

### Construct输出 | Construct Output

- `{output}.vg` - VG格式变异图|VG format variation graph

### Index输出 | Index Output

- `{output}.xg` - XG索引|XG index
- `{output}.gcsa` - GCSA索引|GCSA index
- `{output}.gbwt` - GBWT索引|GBWT index
- `{output}.min` - Minimizer索引|Minimizer index (Giraffe)
- `{output}.dist` - Distance索引|Distance index (Giraffe)

### Giraffe输出 | Giraffe Output

- `{output}.gam` - GAM格式比对结果|GAM format alignments
- `{output}.gaf` - GAF格式比对结果|GAF format alignments

### Deconstruct输出 | Deconstruct Output

- `{output}.vcf` - VCF格式变异文件|VCF format variant file

## 常见问题 | FAQ

### 1. VCF文件需要索引吗？

**问题|Question**: VCF文件需要索引吗？

**答|Answer**: 建议使用bgzip压缩并创建tabix索引。

**解决方案|Solution**:
```bash
# 压缩VCF|Compress VCF
bgzip -c variants.vcf > variants.vcf.gz

# 创建索引|Create index
tabix -p vcf variants.vcf.gz
```

### 2. Giraffe需要哪些索引文件？

**问题|Question**: Giraffe需要哪些索引文件？

**答|Answer**: Giraffe需要XG、minimizer和distance索引。

**解决方案|Solution**:
```bash
# 使用--giraffe选项自动创建所有必需索引
biopytools vg index -i graph.vg -o graph --giraffe
```

### 3. 如何提高比对速度？

**问题|Question**: 如何提高比对速度？

**答|Answer**:
```bash
# 增加线程数|Increase threads
biopytools vg giraffe -g graph -f reads.fq -o alignments.gam -t 48

# 使用GIRAFFE索引|Use GIRAFFE indexes（已优化|already optimized）
```

### 4. 内存不足怎么办？

**问题|Question**: 内存不足怎么办？

**答|Answer**:
```bash
# 构建时指定染色体区域，减少内存使用|Specify chromosome region during construction
biopytools vg construct -r ref.fa -v variants.vcf.gz -o graph.vg -R chr1

# 减少线程数|Reduce threads
biopytools vg giraffe -g graph -f reads.fq -o alignments.gam -t 8
```

### 5. 如何验证变异图是否正确？

**问题|Question**: 如何验证变异图是否正确？

**答|Answer**:
```bash
# 查看图统计信息|View graph statistics
biopytools vg stats -i graph.vg

# 转换为GFA格式查看|Convert to GFA for viewing
vg view -a graph.vg > graph.gfa
```

## 性能优化建议 | Performance Optimization Recommendations

### 1. Construct优化 | Construct Optimization

```bash
# 分染色体构建|Construct by chromosome
for chr in {1..22} X Y; do
    biopytools vg construct -r ref.fa -v chr${chr}.vcf.gz -o chr${chr}.vg -R ${chr}
done

# 然后合并图|Then merge graphs
vg combine -o merged.vg chr*.vg
```

### 2. Index优化 | Index Optimization

```bash
# 只创建需要的索引|Only create needed indexes
# 对于Giraffe比对，只需要GIRAFFE索引
biopytools vg index -i graph.vg -o graph --giraffe
```

### 3. Giraffe优化 | Giraffe Optimization

```bash
# 批量比对|Batch alignment
for reads in *.fq; do
    biopytools vg giraffe -g graph -f $reads -o ${reads%.fq}.gam -t 24 &
done
wait
```

## 依赖项 | Dependencies

### 系统依赖 | System Dependencies

- conda >= 4.0
- VG >= 1.70.0

### Python依赖 | Python Dependencies

- Python >= 3.8
- pathlib
- dataclasses
- logging
- subprocess

## 参考资源 | References

### VG相关

- **VG GitHub**: https://github.com/vgteam/vg
- **VG文档**: https://vgteam.github.io/vg/
- **VG Wiki**: https://github.com/vgteam/vg/wiki
- **Giraffe论文**: https://www.science.org/doi/10.1126/science.abg8871

### 变异图方法

- **变异图综述**: https://www.nature.com/articles/s41587-020-0472-3
- **GFA格式**: https://github.com/GFA-spec/GFA-spec/
- **泛基因组分析**: https://www.nature.com/articles/s41579-020-00452-w

## 更新日志 | Changelog

### v1.0.0 (2026-04-02)

- 初始版本发布
- 支持construct、index、giraffe、deconstruct四个子命令
- 使用conda run调用VG
- 完整的日志分离（stdout/stderr）
- 自动环境检测和验证
- 结果文件自动检查

## 许可证 | License

本工具遵循biopytools项目的许可证。

This tool follows the biopytools project license.

## 联系方式 | Contact

作者|Author: Xiang LI <lixiang117423@gmail.com>

项目地址|Project: https://github.com/your-org/biopytools
