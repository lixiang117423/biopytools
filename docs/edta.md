# EDTA 转座子注释分析模块

**全基因组转座子鉴定、分类和注释工具 | Whole-genome Transposon Identification, Classification, and Annotation Tool**

## 功能概述 | Overview

EDTA转座子注释分析模块基于EDTA软件，提供全基因组转座子的从头鉴定、分类和注释功能。支持单基因组和泛基因组两种分析模式，能够识别LTR、LINE、SINE、TIR、Helitron等各类转座子元件，并生成高质量的转座子注释库。

## 主要特性 | Key Features

- **完整的注释流程**: 转座子识别→分类→过滤→全基因组注释四步骤自动化
- **两种分析模式**: 单基因组注释和泛基因组注释
- **多种转座子类型**: 支持LTR、LINE、SINE、TIR、Helitron等所有主要类型
- **灵活步骤控制**: 支持单独运行任意步骤或完整流程
- **依赖自动检查**: 自动检测EDTA及所有依赖软件
- **标准化日志**: 完整的处理过程日志和进度追踪
- **参数验证**: 智能的输入参数验证和错误提示

## 快速开始 | Quick Start

### 单基因组注释 | Single-genome Annotation

```bash
# 基本用法
biopytools edta -i plant.fa --anno 1 -t 24

# 使用CDS文件提高注释质量
biopytools edta -i plant.fa -c cds.fa --anno 1 -t 24

# 使用已筛选的TE库
biopytools edta -i plant.fa --curatedlib curated_TE.fa --anno 1 -t 24
```

### 泛基因组注释 | Pan-genome Annotation

```bash
# 基本用法
biopytools panedta -i genomes.txt -c cds.fa -t 24

# 使用筛选库
biopytools panedta -i genomes.txt -c cds.fa -l curated_TE.fa -t 24
```

## 参数说明 | Parameters

### 单基因组注释参数 | Single-genome Annotation Parameters

#### 必需参数 | Required Parameters

| 参数 | 描述 | 示例 |
|------|------|------|
| `-i, --genome` | 基因组FASTA文件 | `-i plant.fa` |

#### 基本参数 | Basic Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--species` | `others` | 物种类型(Rice/Maize/others) |
| `--step` | `all` | 运行步骤(all/filter/final/anno) |
| `--overwrite` | `0` | 覆盖已有结果(0/1) |
| `-t, --threads` | `12` | 线程数 |
| `-o, --output-dir` | `./edta_output` | 输出目录 |

#### 输入文件参数 | Input File Parameters

| 参数 | 描述 |
|------|------|
| `--cds` | CDS序列文件，用于去除基因序列 |
| `--curatedlib` | 已筛选的TE库，100%信任 |
| `--rmlib` | RepeatModeler库，增强LINE敏感性 |
| `--exclude` | 排除区域BED文件 |
| `--rmout` | RepeatMasker输出文件 |

#### 分析参数 | Analysis Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--sensitive` | `0` | 使用RepeatModeler识别剩余TE |
| `--anno` | `0` | 执行全基因组TE注释 |
| `--maxdiv` | `40` | 最大分歧度(0-100) |
| `--evaluate` | `0` | 评估注释一致性 |
| `--force` | `0` | 无TE时使用水稻TE继续 |
| `--u` | `1.3e-8` | 中性突变率 |
| `--debug` | `0` | 保留中间文件 |

### 泛基因组注释参数 | Pan-genome Annotation Parameters

#### 必需参数 | Required Parameters

| 参数 | 描述 | 示例 |
|------|------|------|
| `-i, --genome-list` | 基因组列表文件 | `-i genomes.txt` |

#### 可选参数 | Optional Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-c, --cds` | None | CDS序列文件 |
| `-l, --curatedlib` | None | 已筛选的TE库 |
| `-f, --fl-copy` | `3` | 全长拷贝数阈值 |
| `-a, --anno` | `1` | 执行全基因组注释 |
| `--overwrite` | `0` | 覆盖已有结果 |
| `-t, --threads` | `12` | 线程数 |
| `-o, --output-dir` | `./panedta_output` | 输出目录 |

## 输入文件格式 | Input File Formats

### 基因组FASTA文件 | Genome FASTA File

标准FASTA格式的基因组序列文件，要求：
- 序列名称建议≤13字符
- 仅包含字母、数字和下划线
- 序列名称应简洁唯一

### CDS文件 | CDS File

FASTA格式的编码序列文件：
- 不包含内含子、UTR和转座子
- 用于从TE库中去除基因序列

### 基因组列表文件 | Genome List File

纯文本文件，每行一个基因组路径：
```
/path/to/genome1.fa
/path/to/genome2.fa
/path/to/genome3.fa
```

或者每行包含基因组及其对应的CDS：
```
/path/to/genome1.fa /path/to/genome1_cds.fa
/path/to/genome2.fa /path/to/genome2_cds.fa
```

### BED文件 | BED File

标准BED格式的排除区域文件：
```
chr1 1000 5000 gene1
chr2 2000 6000 gene2
```

## 输出结果 | Output Results

### 单基因组注释输出 | Single-genome Annotation Output

主要输出文件：

| 文件 | 描述 |
|------|------|
| `{genome}.mod.EDTA.TElib.fa` | 非冗余TE库 |
| `{genome}.mod.EDTA.TElib.fa.mod.tsv` | TE库统计表 |
| `{genome}.mod.EDTA.TEanno.gff3` | 全基因组TE注释(使用--anno 1) |
| `{genome}.mod.EDTA.TEanno.sum` | 注释统计摘要(使用--anno 1) |
| `{genome}.mod.MAKER.masked` | 仅屏蔽长TE的基因组(使用--anno 1) |

### 泛基因组注释输出 | Pan-genome Annotation Output

输出目录包含：
- 各基因组的独立注释结果
- 泛基因组TE库
- 使用泛基因组库重新注释的结果

## 使用示例 | Usage Examples

### 示例1: 植物基因组完整注释

```bash
biopytools edta \
    -i Arabidopsis.fa \
    -c arabidopsis_cds.fa \
    --species others \
    --sensitive 1 \
    --anno 1 \
    --evaluate 1 \
    -t 24 \
    -o arabidopsis_edta
```

### 示例2: 使用已筛选TE库

```bash
biopytools edta \
    -i rice.fa \
    -c rice_cds.fa \
    --curatedlib rice_curated_TE.fa \
    --species Rice \
    --anno 1 \
    -t 24
```

### 示例3: 仅构建TE库

```bash
biopytools edta \
    -i plant.fa \
    -c plant_cds.fa \
    --anno 0 \
    -t 24
```

### 示例4: 从特定步骤继续

```bash
# 先运行到filter步骤
biopytools edta -i plant.fa --step filter -t 24

# 从final步骤继续
biopytools edta -i plant.fa --step final --overwrite 0 -t 24
```

### 示例5: 泛基因组注释

```bash
# 创建基因组列表文件
cat > genomes.txt << EOF
/path/to/genome1.fa
/path/to/genome2.fa
/path/to/genome3.fa
EOF

# 运行泛基因组注释
biopytools panedta \
    -i genomes.txt \
    -c shared_cds.fa \
    -l curated_TE.fa \
    -f 3 \
    -t 24
```

## 系统要求 | System Requirements

### 依赖软件 | Dependencies

EDTA模块需要以下软件（通过conda安装EDTA环境时自动包含）：

- **EDTA** (v2.2.2或更新)
- **Perl** (v5.10+)
- **RepeatMasker**
- **RepeatModeler**
- **LTR_retriever**
- **LTR_FINDER_parallel**
- **LTR_HARVEST_parallel**
- **TIR-Learner**
- **HelitronScanner**
- **CD-HIT**
- **GenomeTools**
- **BLAST+**
- **TEsorter**
- **BEDTools**
- **SAMtools**

### 安装EDTA环境 | Install EDTA Environment

```bash
# 使用conda安装EDTA
conda create -n EDTA_v.2.2.2
conda activate EDTA_v.2.2.2
mamba install -c conda-forge -c bioconda edta
```

### 硬件建议 | Hardware Recommendations

- **CPU**: 多核处理器（推荐8核以上）
- **RAM**: 最少8GB（大基因组推荐32GB以上）
- **存储**: 预留基因组文件大小10倍的磁盘空间
- **运行时间**: 数小时到数天（取决于基因组大小）

## 注意事项 | Important Notes

1. **序列名称**: 基因组序列名称应保持简洁（≤13字符），避免特殊字符
2. **CDS文件**: 提供CDS文件可显著提高TE库质量，去除基因污染
3. **筛选库**: 仅提供经过手工筛选的高质量TE库
4. **步骤控制**: 可使用--step参数分步执行，便于调试和恢复
5. **长时间运行**: 完整流程可能需要数小时到数天，建议在screen或tmux中运行
6. **磁盘空间**: 确保有足够的磁盘空间存储中间文件和结果

## 故障排除 | Troubleshooting

### 常见问题 | Common Issues

**Q: "EDTA.pl not found" 错误**
```bash
# 检查EDTA环境
conda activate EDTA_v.2.2.2
which EDTA.pl

# 或指定EDTA路径
biopytools edta -i plant.fa --edta-path /path/to/EDTA
```

**Q: 内存不足错误**
```bash
# 减少线程数
biopytools edta -i plant.fa -t 4
```

**Q: 未找到TE候选**
```bash
# 使用--force 1参数
biopytools edta -i plant.fa --force 1
```

**Q: 想跳过某些步骤**
```bash
# 从特定步骤开始
biopytools edta -i plant.fa --step final
```

## 相关资源 | Related Resources

- [EDTA官方文档](https://github.com/oushujun/EDTA)
- [EDTA Wiki](https://github.com/oushujun/EDTA/wiki)
- [EDTA安装指南](https://github.com/oushujun/EDTA#installation)

## 引用信息 | Citation

如果在研究中使用EDTA，请引用：

```
Ou S., Su W., Liao Y., Chougule K., Agda J. R. A., Hellinga A. J.,
Lugo C. S. B., Elliott T. A., Ware D., Peterson T., Jiang N., Hirsch C. N.
and Hufford M. B. (2019). Benchmarking Transposable Element Annotation
Methods for Creation of a Streamlined, Comprehensive Pipeline.
Genome Biol. 20(1): 275.
```

如果使用panEDTA功能，请额外引用：

```
Ou S., Scheben A., Collins T., Qiu Y., Seetharam A., Menard C.,
Manchanda N., Gent J., Schatz M., Anderson S., Hufford M., Hirsch C. (2024).
Differences in activity and stability drive transposable element variation
in tropical and temperate maize. Genome Research.
```
