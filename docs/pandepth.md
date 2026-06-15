# PanDepth 覆盖度计算工具

**超快速基因组覆盖度计算工具 | Ultra-fast Genome Coverage Calculation Tool**

## 功能概述 | Overview

PanDepth 覆盖度计算工具是基于 PanDepth 软件封装的 Python 工具，提供超快速高效的基因组覆盖度计算功能。支持染色体级别、基因级别、特定区域和滑动窗口的覆盖度分析，适用于各种基因组测序数据分析场景。

## 主要特性 | Key Features

- **超高性能**: 比传统工具快数倍，支持多线程加速
- **多种分析模式**: 染色体/基因/区域/窗口覆盖度分析
- **批量处理**: 支持单个BAM文件或BAM文件目录批量处理
- **灵活过滤**: 可配置的比对质量、深度和FLAG过滤
- **多格式支持**: 支持 SAM/BAM/CRAM/PAF 格式
- **无需排序**: 不需要对BAM文件进行排序和索引
- **详细日志**: 完整的处理过程日志和错误追踪

## 快速开始 | Quick Start

### 基本用法 | Basic Usage

```bash
# 单个BAM文件染色体覆盖度分析
biopytools pandepth -i sample.bam -o coverage_results

# 批量处理目录中的所有BAM文件
biopytools pandepth -i bam_files/ -o coverage_results -t 24
```

### 高级用法 | Advanced Usage

```bash
# 基因覆盖度分析（使用GFF/GTF文件）
biopytools pandepth \
    -i sample.bam \
    -g genes.gff \
    -o gene_coverage \
    -t 24

# 特定区域覆盖度分析（使用BED文件）
biopytools pandepth \
    -i sample.bam \
    -b regions.bed \
    -o region_coverage

# 滑动窗口覆盖度分析
biopytools pandepth \
    -i sample.bam \
    -w 1000 \
    -o window_coverage
```

## 参数说明 | Parameters

### 必需参数 | Required Parameters

| 参数 | 描述 | 示例 |
|------|------|------|
| `-i, --input` | 输入BAM文件或BAM文件目录 | `-i sample.bam` 或 `-i bam_files/` |
| `-o, --output` | 输出目录 | `-o coverage_results` |

### 目标区域选项 | Target Region Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-g, --gff` | None | GFF/GTF文件用于基因覆盖度分析 |
| `-f, --feature` | `CDS` | GFF/GTF特征类型 (CDS或exon) |
| `-b, --bed` | None | BED文件用于特定区域覆盖度 |
| `-w, --window` | None | 滑动窗口大小(bp) |

**注意**: 只能指定一个目标区域选项 (-g, -b, -w)

### 过滤选项 | Filter Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-q, --min-mapq` | `0` | 最小比对质量 |
| `-d, --min-depth` | `1` | 最小深度用于统计 |
| `-x, --exclude-flag` | `1796` | 排除reads的FLAG标志 |

### 其他选项 | Other Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-t, --threads` | `12` | 线程数 |
| `-r, --reference` | None | 参考基因组文件(用于CRAM解码或GC计算) |
| `-c, --enable-gc` | False | 启用GC含量计算 |
| `-a, --all-sites` | False | 输出所有位点深度 |
| `--pandepth-path` | `/share/org/YZWL/yzwl_lixg/software/PanDepth-2.26-Linux-x86_64/pandepth` | PanDepth程序路径 |

## 使用示例 | Usage Examples

### 示例1: 染色体级别覆盖度分析

```bash
# 单个样本
biopytools pandepth -i sample.bam -o chr_coverage

# 批量处理
biopytools pandepth -i bam_directory/ -o chr_coverage -t 24
```

**输出文件**: `{sample_name}.chr.stat.gz`

### 示例2: 基因覆盖度分析

```bash
# 使用CDS特征
biopytools pandepth \
    -i sample.bam \
    -g annotation.gff \
    -f CDS \
    -o gene_coverage

# 使用exon特征
biopytools pandepth \
    -i sample.bam \
    -g annotation.gtf \
    -f exon \
    -o gene_coverage
```

**输出文件**: `{sample_name}.gene.stat.gz`

### 示例3: 特定区域覆盖度分析

```bash
# 使用BED文件定义区域
biopytools pandepth \
    -i sample.bam \
    -b target_regions.bed \
    -o region_coverage
```

**输出文件**: `{sample_name}.bed.stat.gz`

### 示例4: 滑动窗口覆盖度分析

```bash
# 1kb窗口
biopytools pandepth \
    -i sample.bam \
    -w 1000 \
    -o window_coverage

# 10kb窗口
biopytools pandepth \
    -i sample.bam \
    -w 10000 \
    -o window_coverage
```

**输出文件**: `{sample_name}.win.stat.gz`

### 示例5: 高质量reads过滤

```bash
# 使用高质量过滤
biopytools pandepth \
    -i sample.bam \
    -o high_quality_coverage \
    -q 20 \
    -d 5
```

### 示例6: 启用GC含量计算

```bash
biopytools pandepth \
    -i sample.bam \
    -o gc_coverage \
    -r genome.fa \
    -c
```

## 输出文件格式 | Output File Format

### 染色体统计文件 (.chr.stat.gz)

```
#Chr	Length	CoveredSite	TotalDepth	Coverage(%)	MeanDepth
Chr01	332615375	331690145	16945298412	99.72	50.95
Chr02	177319215	176853774	8368043949	99.74	47.19
```

### 基因统计文件 (.gene.stat.gz)

```
#Chr	Start	End	GeneID	Length	CoveredSite	TotalDepth	Coverage(%)	MeanDepth
Chr01	16621	30840	Gene001	819	819	38297	100.00	46.76
Chr01	66931	67440	Gene002	510	510	23804	100.00	46.67
```

### 区域统计文件 (.bed.stat.gz)

```
#Chr	Start	End	RegionID	Length	CoveredSite	TotalDepth	Coverage(%)	MeanDepth
Chr01	16621	16662	Chr01_16621_16662	42	42	1935	100.00	46.07
Chr01	16838	16947	Chr01_16838_16947	110	110	5522	100.00	50.20
```

### 窗口统计文件 (.win.stat.gz)

```
#Chr	Start	End	Length	CoveredSite	TotalDepth	Coverage(%)	MeanDepth
Chr01	1	100	100	0	0	0.00	0.00
Chr01	101	200	100	0	0	0.00	0.00
Chr01	201	300	100	26	26	26.00	0.26
```

## 性能优势 | Performance Advantages

相比传统工具（samtools, bedtools等）：

- **速度**: 比samtools depth快3-5倍
- **内存**: 内存占用更低
- **便捷**: 无需BAM文件排序和索引
- **准确**: 与samtools depth结果完全一致

## 系统要求 | System Requirements

### 依赖软件 | Dependencies

- **PanDepth** (v2.26或更新)
- **Python** (版本 3.7+)
- **Python包**: `click`, `pathlib`, `subprocess`

### 硬件建议 | Hardware Recommendations

- **CPU**: 多核处理器（推荐4核以上）
- **RAM**: 最少4GB（大基因组推荐8GB以上）
- **存储**: 根据BAM文件大小预留足够空间

## 注意事项 | Important Notes

1. **目标区域选项**: 只能指定一个(-g, -b, -w)，不能同时使用
2. **文件格式**: 支持BAM/SAM/CRAM/PAF格式
3. **排序要求**: 不需要对BAM文件进行排序
4. **索引要求**: 不需要索引文件，但如果有索引可以利用多线程加速
5. **内存使用**: 大基因组文件可能需要较多内存

## 故障排除 | Troubleshooting

### 常见问题 | Common Issues

**Q: "PanDepth not found" 错误**
```bash
# 检查PanDepth程序是否存在
ls -l /share/org/YZWL/yzwl_lixg/software/PanDepth-2.26-Linux-x86_64/pandepth

# 使用--pandepth-path指定正确路径
biopytools pandepth ... --pandepth-path /path/to/pandepth
```

**Q: 找不到BAM文件**
```bash
# 检查输入路径
ls -l your_bam_directory/*.bam

# 确保使用正确的路径
biopytools pandepth -i /full/path/to/bam_file.bam -o output
```

**Q: 内存不足**
```bash
# 减少线程数
biopytools pandepth -i sample.bam -o output -t 4

# 或分批处理BAM文件
```

## 相关资源 | Related Resources

- [PanDepth GitHub](https://github.com/HuiyangYu/PanDepth)
- [PanDepth 论文](https://doi.org/10.1093/bib/bbae197)
- [PMID: 38701418](https://pubmed.ncbi.nlm.nih.gov/38701418/)

## 引用信息 | Citation

如果在学术研究中使用此工具，请引用PanDepth相关文献：

```
Yu, H., et al. (2024).
PanDepth: an ultrafast and efficient tool for coverage calculation.
Briefings in Bioinformatics, 25(3), bbae197.
PMID: 38701418
```

## 许可证 | License

本项目采用MIT许可证

**注意**: PanDepth软件本身采用MIT-0许可证
