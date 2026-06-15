# HiFi+Hi-C 基因组组装与挂载流程模块

**专业的植物基因组完整组装工具 | Professional Plant Genome Complete Assembly Tool**

## 功能概述 | Overview

HiFi+Hi-C工作流模块是一个整合了4个子模块的完整基因组组装流程，专门用于植物基因组的从头组装。该模块采用最新的HiFi长读长测序技术和Hi-C染色体构象捕获技术，结合参考基因组引导命名，实现了从原始测序数据到染色体级别高质量基因组的自动化组装。

该流程整合了：
1. **hifi_hic**: HiFi组装（可选NGS polish）
2. **haphic**: Hi-C染色体挂载与聚类
3. **rename_chromosomes**: 参考基因组引导的染色体标准化命名
4. **hic_heatmap**: 全基因组染色体间Hi-C热图生成

## 主要特性 | Key Features

- **完整流程整合**: 4个步骤无缝集成，自动化数据传递
- **NGS Polish支持**: 可选NGS数据polish，提高组装准确性
- **参考基因组命名**: 基于minimap2比对的智能染色体命名
- **断点续传**: 支持从任意断点恢复执行
- **灵活控制**: 可跳过任意步骤，支持部分执行
- **详细日志**: 完整的执行日志和错误追踪
- **双语支持**: 中英文双语注释和日志

## 快速开始 | Quick Start

### 基本用法 | Basic Usage

```bash
# 标准流程（无NGS polish）
biopytools hifi-hic-workflow \
    --hifi hifi_reads.fq.gz \
    --hic-r1 hic_R1.fq.gz \
    --hic-r2 hic_R2.fq.gz \
    --ref reference.fa \
    -o ./workflow_output

# 使用NGS polish
biopytools hifi-hic-workflow \
    --hifi hifi_reads.fq.gz \
    --hic-r1 hic_R1.fq.gz \
    --hic-r2 hic_R2.fq.gz \
    --ref reference.fa \
    --use-ngs-polish \
    --ngs-data ngs_data_dir \
    -o ./workflow_output
```

### 高级用法 | Advanced Usage

```bash
# 自定义参数
biopytools hifi-hic-workflow \
    --hifi hifi_reads.fq.gz \
    --hic-r1 hic_R1.fq.gz \
    --hic-r2 hic_R2.fq.gz \
    --ref reference.fa \
    -o ./workflow_output \
    -p sample1 \
    -t 64 \
    --nchrs 20 \
    --naming-min-identity 90.0 \
    --heatmap-resolution 500000

# 跳过某些步骤
biopytools hifi-hic-workflow \
    --hifi hifi_reads.fq.gz \
    --hic-r1 hic_R1.fq.gz \
    --hic-r2 hic_R2.fq.gz \
    --ref reference.fa \
    -o ./workflow_output \
    --skip-heatmap

# 强制重新运行
biopytools hifi-hic-workflow \
    --hifi hifi_reads.fq.gz \
    --hic-r1 hic_R1.fq.gz \
    --hic-r2 hic_R2.fq.gz \
    --ref reference.fa \
    -o ./workflow_output \
    --force
```

## 参数说明 | Parameters

### 必需参数 | Required Parameters

| 参数 | 描述 | 示例 |
|------|------|------|
| `--hifi` | HiFi reads文件 | `--hifi hifi_reads.fq.gz` |
| `--hic-r1` | Hi-C R1文件 | `--hic-r1 hic_R1.fq.gz` |
| `--hic-r2` | Hi-C R2文件 | `--hic-r2 hic_R2.fq.gz` |
| `--ref`, `--reference` | 参考基因组（仅用于命名） | `--ref reference.fa` |
| `-o`, `--output` | 输出目录 | `-o ./output` |

### 全局参数 | Global Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-p`, `--prefix` | `genome_sample` | 样本前缀 |
| `-t`, `--threads` | `64` | 线程数 |
| `-v`, `--verbose` | `False` | 显示详细日志 |

### 流程控制 | Workflow Control

| 参数 | 描述 |
|------|------|
| `--skip-hifi-hic` | 跳过HiFi组装 |
| `--skip-haphic` | 跳过Hi-C挂载 |
| `--skip-rename` | 跳过重命名 |
| `--skip-heatmap` | 跳过热图 |
| `--no-resume` | 禁用断点续传 |
| `--force` | 强制重新运行所有步骤 |

### Step 1: HiFi组装参数 | Step 1: HiFi Assembly Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--genome-size` | `1.45g` | 预估基因组大小 |
| `--n-hap` | `2` | 倍性 |
| `--purge-level` | `None` | Purge level (0/1/2/3) |
| `--hom-cov` | `None` | Homozygous coverage |

### NGS Polish参数 | NGS Polish Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--use-ngs-polish` | `False` | 启用NGS polish |
| `--ngs-data` | `None` | NGS数据目录 |
| `--ngs-high-cov` | `95.0` | 高质量覆盖度阈值(%) |
| `--ngs-pattern` | `_1.clean.fq.gz` | NGS文件匹配模式 |

### Step 2: HapHiC参数 | Step 2: HapHiC Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--nchrs` | `自动统计` | 染色体数量 |
| `--haphic-bin` | `haphic` | HapHiC可执行文件 |
| `--bwa-bin` | `bwa` | BWA可执行文件 |
| `--samtools-bin` | `samtools` | Samtools可执行文件 |

### Step 3: 染色体重命名参数 | Step 3: Chromosome Rename Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--rename-keep-all` | `True` | 保留所有序列 |
| `--naming-min-identity` | `80.0` | 最小序列一致性(%) |
| `--naming-min-coverage` | `80.0` | 最小覆盖度(%) |
| `--naming-minimap2-preset` | `asm5` | minimap2预设 |

### Step 4: Hi-C热图参数 | Step 4: Hi-C Heatmap Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--hicpro-enzyme` | `MboI` | 限制性内切酶 |
| `--heatmap-resolution` | `100000` | 热图分辨率(bp) |
| `--heatmap-colormap` | `YlOrRd` | 颜色方案 |
| `--heatmap-format` | `pdf` | 输出格式 |

## 输入文件格式 | Input File Formats

### HiFi Reads
标准FASTQ格式的HiFi长读长测序数据：
```
@read_id
GATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
+
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
```

### Hi-C Reads
双端测序FASTQ文件，支持gzip压缩。

### 参考基因组
标准FASTA格式的参考基因组序列（用于命名指导）：
```
>chr1
ATGGCGATCGATCGATCGATCGATCGATCGATCGATCG
>chr2
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
```

## 输出结果 | Output Results

### 目录结构 | Directory Structure

```
{output_dir}/
├── 00.pipeline_info/          # 流程信息和报告
│   ├── hifi_hic_info.txt
│   ├── haphic_info.txt
│   ├── rename_info.txt
│   └── workflow_report.txt
├── 01.hifi_assembly/         # HiFi组装结果
│   └── {prefix}/
│       ├── 01.raw_output/    # Hifiasm原始输出
│       ├── 02.fasta/         # FASTA转换结果
│       ├── 03.ngs_polish/    # NGS polish结果（可选）
│       └── 04.statistics/    # 组装统计
├── 02.hic_scaffolding/       # HapHiC挂载结果
│   ├── 00.mapping/           # BWA比对结果
│   ├── 01.clustering/        # 聚类分析
│   ├── 02.ordering/          # 排序结果
│   ├── 03.correction/        # 组装校正
│   └── 04.build/             # 最终scaffolds
│       └── {prefix}.fa
├── 03.chromosome_rename/     # 染色体重命名结果
│   ├── alignment/            # minimap2比对
│   ├── naming_report.txt     # 命名报告
│   └── {prefix}.renamed.fa   # 重命名基因组
├── 04.hic_heatmap/           # Hi-C热图结果
│   ├── hicpro_output/        # HiCPro分析结果
│   ├── plot/                 # 热图文件
│   │   └── {prefix}_hic_heatmap.{format}
│   └── logs/
└── logs/                     # 流程日志
    └── workflow_*.log
```

### 主要输出文件 | Main Output Files

| 文件 | 描述 |
|------|------|
| `01.hifi_assembly/{prefix}/02.fasta/{prefix}.primary.fa` | HiFi组装primary contigs |
| `01.hifi_assembly/{prefix}/03.ngs_polish/04.reassembly/02.fasta/{prefix}.primary.fa` | NGS polish后基因组（如果启用） |
| `02.hic_scaffolding/04.build/{prefix}.fa` | HapHiC挂载后的染色体级别scaffolds |
| `03.chromosome_rename/{prefix}.renamed.fa` | 标准命名的最终基因组 |
| `04.hic_heatmap/plot/{prefix}_hic_heatmap.pdf` | 全基因组染色体间Hi-C热图 |

## 使用示例 | Usage Examples

### 示例1: 标准植物基因组组装

```bash
biopytools hifi-hic-workflow \
    --hifi pacbio_hifi.fq.gz \
    --hic-r1 arri_hic_R1.fq.gz \
    --hic-r2 arri_hic_R2.fq.gz \
    --ref related_species.fa \
    -o ./plant_genome \
    -p plant1 \
    -t 64
```

### 示例2: 使用NGS Polish提高组装质量

```bash
biopytools hifi-hic-workflow \
    --hifi pacbio_hifi.fq.gz \
    --hic-r1 arri_hic_R1.fq.gz \
    --hic-r2 arri_hic_R2.fq.gz \
    --ref related_species.fa \
    --use-ngs-polish \
    --ngs-data /data/illumina_reads \
    --ngs-high-cov 90.0 \
    -o ./plant_genome_polished
```

### 示例3: 仅执行部分步骤

```bash
# 如果HiFi组装已完成，从HapHiC开始
biopytools hifi-hic-workflow \
    --hifi pacbio_hifi.fq.gz \
    --hic-r1 arri_hic_R1.fq.gz \
    --hic-r2 arri_hic_R2.fq.gz \
    --ref related_species.fa \
    --skip-hifi-hic \
    -o ./workflow_output
```

### 示例4: 自定义染色体命名参数

```bash
biopytools hifi-hic-workflow \
    --hifi pacbio_hifi.fq.gz \
    --hic-r1 arri_hic_R1.fq.gz \
    --hic-r2 arri_hic_R2.fq.gz \
    --ref related_species.fa \
    --naming-min-identity 90.0 \
    --naming-min-coverage 85.0 \
    --naming-minimap2-preset asm10 \
    -o ./workflow_output
```

### 示例5: 生成高分辨率Hi-C热图

```bash
biopytools hifi-hic-workflow \
    --hifi pacbio_hifi.fq.gz \
    --hic-r1 arri_hic_R1.fq.gz \
    --hic-r2 arri_hic_R2.fq.gz \
    --ref related_species.fa \
    --heatmap-resolution 50000 \
    --heatmap-colormap RdYlBu \
    --heatmap-format png \
    -o ./workflow_output
```

## 流程说明 | Workflow Details

### Step 1: HiFi组装 (hifi_hic)

**功能**: 使用hifiasm进行HiFi长读长基因组组装

**输入**:
- HiFi reads (必需)
- Hi-C reads (可选，用于更好的单倍型分离)
- NGS reads (可选，用于NGS polish)

**流程**:
1. Hifiasm组装生成GFA文件
2. GFA转FASTA格式
3. 生成contig-reads映射文件
4. (可选) BWA比对NGS数据
5. (可选) 覆盖度过滤筛选高质量contigs
6. (可选) 提取高质量HiFi reads
7. (可选) 使用筛选后的reads重新组装

**输出**: primary.fa (或polished.fa)

### Step 2: Hi-C挂载 (haphic)

**功能**: 使用Hi-C数据进行染色体级别scaffolding

**输入**:
- HiFi组装的基因组
- Hi-C reads (FASTQ格式)
- 染色体数量

**流程**:
1. BWA比对Hi-C reads到基因组
2. Hi-C数据过滤
3. 染色体聚类
4. contig排序和定向
5. 组装校正
6. 生成最终scaffolds

**输出**: {prefix}.fa (染色体级别scaffolds)

### Step 3: 染色体重命名 (rename_chromosomes)

**功能**: 参考基因组引导的染色体标准化命名

**输入**:
- HapHiC挂载的scaffolds
- 参考基因组
- 命名阈值

**流程**:
1. 使用minimap2将scaffolds比对到参考基因组
2. 分析比对结果，提取最佳匹配
3. 根据匹配结果重命名序列
4. 生成命名报告

**命名规则**:
- 一致性≥80% 且 覆盖度≥80% → 使用参考染色体名称
- 否则 → 保留原名

**输出**: {prefix}.renamed.fa

### Step 4: Hi-C热图 (hic_heatmap)

**功能**: 生成全基因组染色体间Hi-C热图

**输入**:
- 重命名后的基因组
- Hi-C reads

**流程**:
1. HiCPro比对和分析
2. 生成接触矩阵
3. PlotHiC可视化

**输出**: {prefix}_hic_heatmap.{format}

## 系统要求 | System Requirements

### 依赖软件 | Dependencies

**必需软件**:
- **hifiasm** (≥0.19)
  - 安装: `conda install -c bioconda hifiasm`
- **bwa** (≥0.7.17)
  - 安装: `conda install -c bioconda bwa`
- **samtools** (≥1.15)
  - 安装: `conda install -c bioconda samtools`
- **minimap2** (≥2.24)
  - 安装: `conda install -c bioconda minimap2`

**Python依赖**:
- Python ≥ 3.8
- pandas
- numpy

**可选软件**:
- **ntCard** (用于NGS polish)
  - 安装: `conda install -c bioconda ntcard`
- **HiC-Pro** (用于Hi-C热图)
  - 安装: `conda install -c bioconda hicpro`
- **PlotHiC** (用于热图可视化)

### 硬件建议 | Hardware Recommendations

| 数据类型 | CPU | 内存 | 存储 |
|----------|-----|------|------|
| HiFi (100x) | 32-64核 | 128-256 GB | HiFi数据大小的10倍 |
| Hi-C (150x) | 32-64核 | 64-128 GB | Hi-C数据大小的10倍 |
| 完整流程 | 64-128核 | 256-512 GB | 输入数据总和的15倍 |

## 注意事项 | Important Notes

### 数据质量控制

1. **HiFi数据**:
   - 推荐长度N50 ≥ 20 kb
   - 推荐平均准确性 ≥ QV30
   - 推荐测序深度 ≥ 80x

2. **Hi-C数据**:
   - 推荐测序深度 ≥ 150x
   - 推荐insert size ≥ 10 kb
   - 需要成对的R1/R2文件

3. **NGS数据** (可选):
   - 推荐测序深度 ≥ 100x
   - 推荐读长 ≥ 150 bp
   - 推荐insert size = 350 bp

### 参数选择建议

**基因组大小 (--genome-size)**:
- 小于1G: 使用 "m" 单位 (例如: 800m)
- 1-10G: 使用 "g" 单位 (例如: 1.45g)
- 大于10G: 使用 "g" 单位 (例如: 25g)

**倍性 (--n-hap)**:
- 二倍体: 2
- 多倍体: 相应倍性 (例如: 6)

**Purge level (--purge-level)**:
- 0: 不purge (保留所有冗余)
- 1: 轻度purge (保留主要冗余)
- 2/3: 激进purge (默认unzip模式)

**命名阈值**:
- 严格命名: identity ≥ 90%, coverage ≥ 85%
- 标准命名: identity ≥ 80%, coverage ≥ 80% (默认)
- 宽松命名: identity ≥ 70%, coverage ≥ 70%

### 常见问题

**Q: 内存不足错误**
```bash
# 减少线程数
biopytools hifi-hic-workflow ... -t 32

# 或禁用某些内存密集型步骤
biopytools hifi-hic-workflow ... --skip-heatmap
```

**Q: Hi-C挂载效果差**
```bash
# 检查Hi-C数据质量
# 增加Hi-C测序深度
# 调整聚类参数
biopytools hifi-hic-workflow ... --nchrs 30
```

**Q: 染色体命名不理想**
```bash
# 调整命名阈值
biopytools hifi-hic-workflow ... \
    --naming-min-identity 90.0 \
    --naming-min-coverage 85.0
```

## 故障排除 | Troubleshooting

### Step 1: HiFi组装问题

**Hifiasm运行失败**
- 检查HiFi数据格式和完整性
- 增加内存分配
- 检查线程数是否过多

**NGS polish失败**
- 检查NGS数据目录结构
- 检查文件命名模式
- 降低high_cov阈值

### Step 2: HapHiC问题

**BWA比对失败**
- 检查Hi-C FASTQ文件格式
- 检查磁盘空间
- 减少线程数

**聚类效果差**
- 检查Hi-C数据质量
- 调整聚类参数 (--min-inflation, --max-inflation)
- 检查染色体数估计

### Step 3: 染色体重命名问题

**minimap2比对失败**
- 检查参考基因组格式
- 检查输入基因组格式
- 调整minimap2预设

**命名结果不理想**
- 调整命名阈值
- 检查参考基因组质量
- 查看命名报告

### Step 4: Hi-C热图问题

**HiCPro运行失败**
- 检查Hi-C数据
- 检查磁盘空间
- 增加内存限制

**PlotHiC无输出**
- 检查HiCPro输出
- 检查分辨率设置
- 查看日志文件

## 性能优化建议 | Performance Optimization

### 并行化

- **HiFi组装**: 使用多线程加速 (--threads)
- **BWA比对**: 每个线程约需2GB内存
- **HapHiC**: 调整processes参数优化速度
- **HiCPro**: 根据内存调整max_memory

### 内存优化

- **Hifiasm**: 减少线程数降低内存需求
- **BWA**: 限制并发比对任务
- **HiCPro**: 调整SORT_RAM参数

### 磁盘I/O

- 使用SSD存储临时文件
- 减少日志详细程度
- 定期清理中间文件

## 与其他模块配合 | Integration with Other Modules

### 完整基因组分析流程

```bash
# 1. HiFi+Hi-C组装 (本模块)
biopytools hifi-hic-workflow ...

# 2. 基因组注释
biopytools braker -g renamed.fa -o annotation

# 3. 重复序列注释
biopytools repeatmask -g renamed.fa -o repeats

# 4. 基因组评估
biopytools busco -i renamed.fa -l embryophyta

# 5. 共线性分析
biopytools ragtag ...
```

### 与单独模块的对比

| 模块 | 适用场景 | 优势 |
|------|----------|------|
| **hifi_hic_workflow** | 完整流程 | 一步到位，自动化程度高 |
| **hifi_hic** | 仅需HiFi组装 | 灵活控制组装参数 |
| **haphic** | 已有基因组，需挂载 | 专注于scaffolding |
| **rename_chromosomes** | 仅需重命名 | 快速命名 |
| **hic_heatmap** | 仅需热图 | 灵活可视化 |

## 开发信息 | Development Information

### 代码结构

```
biopytools/hifi_hic_workflow/
├── __init__.py           # 模块导出
├── config.py             # 配置管理
├── data_transfer.py      # 数据传递
├── reference_namer.py    # 参考基因组命名
├── workflow.py           # 流程编排
└── main.py               # 入口和CLI
```

### 开发规范遵循

本模块严格遵循以下开发规范：
- **Section 12.2**: 输出目录结构 (00_pipeline_info, 01_, 02_, etc.)
- **Section 13.3**: Conda环境软件调用 (build_conda_command, timeout=60)
- **双语注释**: 中文|English
- **类型提示**: 完整的类型标注
- **错误处理**: 完善的异常处理和日志记录

## 引用信息 | Citation

如果在学术研究中使用本工作流，请引用各子模块：

**HiFi Assembly (hifiasm)**:
```
Cheng H, Concepcion GT, Feng X, et al.
Haplotype-resolved de novo assembly using hifiasm.
Nature Methods. 2021. doi: 10.1038/s41592-021-01243-7
```

**Hi-C Scaffolding (HapHiC)**:
```
Zhou C, Wu Y, Sevim I, et al.
HapHiC: scaffolding diploid genomes with Hi-C data.
Nature Plants. 2019. doi: 10.1038/s41477-019-0500-6
```

**Hi-C Heatmap (HiC-Pro + PlotHiC)**:
```
Servant N, Varoquaux N, Lajoie BR, et al.
HiC-Pro: an optimized and flexible pipeline for Hi-C data processing.
Genome Biology. 2015. doi: 10.1186/s13059-015-0630-x
```

## 许可证 | License

本项目采用MIT许可证 - 详见 [LICENSE](LICENSE) 文件

**注意**: 各子模块遵循各自的许可证条款。

## 相关资源 | Related Resources

- [hifiasm GitHub](https://github.com/chhylp123/hifiasm)
- [HapHiC GitHub](https://github.com/dbrandan/HapHiC)
- [HiC-Pro GitHub](https://github.com/nservant/HiC-Pro)
- [minimap2 GitHub](https://github.com/lh3/minimap2)

---

**版本信息**: hifi_hic_workflow模块版本 1.0.0 | Module Version 1.0.0
