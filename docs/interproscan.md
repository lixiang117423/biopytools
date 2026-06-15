# InterProScan 蛋白质功能注释模块

**InterProScan Protein Function Annotation Module**

## 功能概述 | Overview

InterProScan蛋白质功能注释模块是对蛋白质序列进行功能结构域注释和GO术语预测的专业工具。基于InterProScan软件构建，提供蛋白质序列功能域预测、GO术语注释、Pathway分析等完整功能。适用于各种蛋白质组学和基因组学注释研究。

## 主要特性 | Key Features

- **全面的功能注释**: 支持18个主流蛋白质数据库(Pfam, SMART, Gene3D, PANTHER等)的注释
- **GO术语预测**: 自动提取Gene Ontology术语，提供分子功能、生物过程和细胞组分类注释
- **Pathway分析**: 支持Pathway数据库注释，分析蛋白质在代谢通路中的作用
- **灵活的输出格式**: 支持TSV、GFF3、XML、HTML、JSON等多种输出格式
- **离线运行支持**: 支持禁用在线预计算服务，完全本地化运行(无需网络)
- **批量序列处理**: 高效处理大规模蛋白质序列集
- **可配置数据库**: 可选择性运行特定数据库应用，节省计算时间
- **详细日志记录**: 完整的运行过程日志和错误追踪

## 快速开始 | Quick Start

### 基本用法 | Basic Usage

```bash
# 基本蛋白质注释
biopytools interproscan -i proteins.fa -o results

# 指定线程数
biopytools interproscan -i proteins.fa -o results -t 32

# 指定特定数据库
biopytools interproscan -i proteins.fa -o results -appl Pfam,SMART,Gene3D
```

### 高级用法 | Advanced Usage

```bash
# 完整功能注释(GO + Pathway)
biopytools interproscan -i proteins.fa -o results --pathways -t 32

# 仅使用特定数据库快速注释
biopytools interproscan -i proteins.fa -o results -appl Pfam -f TSV -t 16

# 启用在线预计算服务(需要网络)
biopytools interproscan -i proteins.fa -o results --enable-precalc
```

## 参数说明 | Parameters

### 必需参数 | Required Parameters

| 参数 | 描述 | 示例 |
|------|------|------|
| `-i, --input` | 输入蛋白质FASTA文件路径 | `-i proteins.fa` |
| `-o, --output-prefix` | 输出文件前缀(不含扩展名) | `-o results` |

### 软件配置 | Software Configuration

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-a, --interproscan-path` | `/share/org/YZWL/yzwl_lixg/software/InterProScan/v.5.75-106.0/interproscan-5.75-106.0/interproscan.sh` | InterProScan软件安装路径 |
| `-f, --format` | `TSV` | 输出格式 (TSV/GFF3/XML/HTML/JSON/TXT) |
| `-t, --threads` | `24` | 线程数 |

### 处理控制选项 | Processing Control Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--disable-precalc` | `True` | 禁用在线预计算查找服务(解决网络问题) |
| `--enable-precalc` | `False` | 启用在线预计算查找服务(需要网络连接) |
| `--goterms` | `True` | 获取GO术语注释 |
| `--no-goterms` | `False` | 不获取GO术语注释 |
| `--pathways` | `False` | 获取Pathway信息 |
| `--applications` | `全部` | 指定运行的应用(逗号分隔) |
| `--temp-dir` | `自动` | 临时目录路径 |

## 输入文件格式 | Input File Format

### 蛋白质FASTA文件 | Protein FASTA File

标准FASTA格式的蛋白质序列文件:

```fasta
>protein1
MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAG
QEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDL
PSRTVDTKQAQDLARSYGIPFIETSAKTRQRVEDAFYTLVREIRQYRLKKISKEEKTPG
>protein2
MSKEKLFKRAGLRSPPRPPRPSPSSPLKSCSPTRAPRRPRPSQRPPLPGLGPVRRRVRR
```

**文件要求**:
- 标准FASTA格式
- 氨基酸使用单字母代码
- 序列标识符以`>`开头

## 可用数据库应用 | Available Applications

| 应用 | 数据库 | 描述 |
|------|--------|------|
| Pfam | Pfam | 蛋白质家族数据库 |
| SMART | SMART | 功能域注释 |
| Gene3D | Gene3D | 结构域数据库 |
| PANTHER | PANTHER | 蛋白质分类和进化 |
| ProSiteProfiles | ProSite | 蛋白质特征和谱 |
| ProSitePatterns | ProSite | 蛋白质模式 |
| PRINTS | PRINTS | 蛋白质指纹 |
| SUPERFAMILY | SUPERFAMILY | 超家族注释 |
| CDD | CDD | 保守结构域 |
| Coils | Coils | 卷曲螺旋预测 |
| Hamap | Hamap | 高质量自动化注释 |
| TIGRFAMs | TIGRFAMs | 蛋白质家族 |
| PIRSF | PIRSF | 超家族 |
| FunFam | FunFam | 功能家族 |
| SFLD | SFLD | 底物特异性家族 |
| AntiFam | AntiFam | 假基因污染检测 |
| MobiDBLite | MobiDBLite | 固有无序区域 |
| NCBIfam | NCBIfam | NCBI保守结构域 |

## 使用示例 | Usage Examples

### 示例1：基本注释流程 | Example 1: Basic Annotation Pipeline

```bash
# 对植物蛋白质组进行基本注释
biopytools interproscan \
    -i Arabidopsis_proteins.fa \
    -o arabidopsis_annotation \
    -t 32
```

### 示例2：使用特定数据库快速注释 | Example 2: Fast Annotation with Specific Databases

```bash
# 仅使用Pfam数据库进行快速注释
biopytools interproscan \
    -i proteins.fa \
    -o pfam_only_results \
    -appl Pfam \
    -t 16
```

### 示例3：完整功能注释(GO + Pathway) | Example 3: Complete Annotation (GO + Pathway)

```bash
# 完整注释包括GO术语和Pathway信息
biopytools interproscan \
    -i bacterial_proteins.fa \
    -o complete_annotation \
    --pathways \
    -t 48
```

### 示例4：离线运行(禁用在线服务) | Example 4: Offline Mode (Disable Online Service)

```bash
# 完全本地化运行，不需要网络连接
biopytools interproscan \
    -i proteins.fa \
    -o offline_results \
    --disable-precalc \
    -t 24
```

### 示例5：小规模序列快速测试 | Example 5: Small Scale Quick Test

```bash
# 使用少量序列测试参数设置
biopytools interproscan \
    -i test_proteins.fa \
    -o test_results \
    -appl Pfam,SMART \
    -t 8
```

## 输出结果 | Output Results

### TSV格式输出 | TSV Format Output

```
ProteinID	MD5	Length	Database	Accession	Start	Stop	EC	GO	InterProAccession	InterProDescription	Pathways
protein1	e7d735c7e8d4e1f9	358	Pfam	PF00071	4	166	-	GO:0004674;GO:0005524	IPR001245	Protein kinase	-
protein1	e7d735c7e8d4e1f9	358	SMART	SM00222	4	166	-	GO:0004674	IPR001245	Protein kinase	-
protein2	9a8b7c6d5e4f3a2b1	252	Pfam	PF00069	10	150	-	GO:0005506;GO:0016491	IPR018091	Esterase	-
```

### GFF3格式输出 | GFF3 Format Output

```
##gff-version 3
##sequence-region protein1 1 358
protein1	InterProScan	PF00071	4	166	.	.	.	ID=protein1_PF00071;Target=PF00071 4 166
protein1	InterProScan	SM00222	4	166	.	.	.	ID=protein1_SM00222;Target=SM00222 4 166
protein2	InterProScan	PF00069	10	150	.	.	.	ID=protein2_PF00069;Target=PF00069 10 150
```

### 输出字段说明 | Output Field Descriptions

| 字段 | 描述 |
|------|------|
| ProteinID | 蛋白质序列标识符 |
| MD5 | 序列的MD5校验和 |
| Length | 蛋白质序列长度 |
| Database | 注释数据库 |
| Accession | 数据库登录号 |
| Start/Stop | 结构域起始/结束位置 |
| EC | EC酶编号 |
| GO | Gene Ontology术语 |
| InterProAccession | InterPro登录号 |
| InterProDescription | InterPro描述 |
| Pathways | 代谢通路信息 |

## 系统要求 | System Requirements

### 依赖软件 | Dependencies

- **InterProScan** (版本 5.75-106.0 或更新)
  - 下载地址: https://github.com/ebi-pf-team/interproscan
- **Java** (版本 11 或更新)
  - InterProScan运行环境
- **Python** (版本 3.7+)
- **Python包**:
  - `click` - 命令行界面
  - `pathlib` - 路径处理

### 硬件建议 | Hardware Recommendations

- **CPU**: 多核处理器(推荐8核以上)
- **RAM**: 最少8GB(大规模蛋白质组推荐32GB以上)
- **存储**: 预留输入文件大小10倍的磁盘空间
- **网络**: 如需使用在线预计算服务，需要稳定网络

## 注意事项 | Important Notes

1. **内存使用**: 大规模蛋白质序列注释会消耗大量内存，建议分批处理
2. **网络连接**: 默认禁用在线预计算服务，避免网络问题导致运行失败
3. **输出格式**: TSV格式适合下游分析，GFF3格式适合基因组浏览器可视化
4. **数据库选择**: 选择特定数据库可显著减少运行时间
5. **线程数**: 建议根据可用CPU核心数和内存大小调整线程数

## 故障排除 | Troubleshooting

### 常见问题 | Common Issues

**Q: "Java heap space" 错误**
```bash
# 增加Java内存限制
export JAVA_OPTS="-Xmx100G"
biopytools interproscan ...
```

**Q: 网络连接超时错误**
```bash
# 禁用在线预计算服务(默认已禁用)
biopytools interproscan ... --disable-precalc
```

**Q: 内存不足错误**
```bash
# 减少线程数
biopytools interproscan ... -t 8

# 或分批处理序列
# 将大文件拆分为小文件后分别运行
```

**Q: 运行时间过长**
```bash
# 仅使用主要数据库
biopytools interproscan ... -appl Pfam,SMART,Gene3D

# 减少序列数量
```

## 结果解读指南 | Result Interpretation Guide

### GO术语分类 | GO Term Categories

| GO类别 | 描述 | 示例 |
|--------|------|------|
| **分子功能(Molecular Function)** | 蛋白质的分子活动 | ATP结合、酶活性 |
| **生物过程(Biological Process)** | 参与的生物过程 | 细胞分裂、信号转导 |
| **细胞组分(Cellular Component)** | 蛋白质定位的细胞结构 | 细胞核、线粒体 |

### 结构域重要性评估 | Domain Importance Assessment

1. **高置信度**: 多个数据库共同注释的结构域
2. **中置信度**: 单个数据库注释，但有GO术语支持
3. **低置信度**: 仅单个数据库注释

## 相关资源 | Related Resources

- [InterProScan官方文档](https://github.com/ebi-pf-team/interproscan/wiki)
- [InterPro数据库](https://www.ebi.ac.uk/interpro/)
- [GO术语数据库](http://geneontology.org/)
- [蛋白质功能注释最佳实践](https://www.nature.com/articles/nrg3933)

## 引用信息 | Citation

如果在学术研究中使用此工具，请引用InterProScan相关文献:

```
Jones, P. et al. (2014)
InterProScan 5: genome-scale protein function classification.
Bioinformatics, 30(9), 1236-1240.
```
