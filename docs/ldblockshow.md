# 连锁不平衡热图分析工具

**基于VCF文件生成连锁不平衡(LD)热图的专业工具 | Professional Tool for Generating Linkage Disequilibrium (LD) Heatmap from VCF Files**

## 功能概述 | Overview

连锁不平衡热图分析工具是基于LDBlockShow软件的Python封装，用于从VCF文件快速生成连锁不平衡热图。该工具支持多种LD度量方法（D'和R²）、LD block检测、GWAS结果可视化和基因组注释展示，适用于连锁不平衡分析、GWAS结果可视化和基因组区域LD结构分析等多种研究场景。

## 主要特性 | Key Features

- **高效计算**: 基于C++实现，处理速度快，内存占用低
- **多种LD度量**: 支持D'、R²或两者同时计算
- **灵活Block检测**: Gabriel方法、Solid spine、自定义cutoff等多种检测策略
- **集成可视化**: 支持GWAS P值和GFF3基因注释同时展示
- **多种输出格式**: SVG、PNG、PDF格式输出
- **质量控制**: 支持MAF、Missing、HWE、Het过滤
- **亚群分析**: 支持亚群样本的LD分析

## 快速开始 | Quick Start

### 基本用法 | Basic Usage

```bash
# 基本LD热图分析
biopytools ldblockshow \
    -i variants.vcf.gz \
    -o results/ \
    -r chr1:1000000-2000000

# 使用R²统计量
biopytools ldblockshow \
    -i variants.vcf.gz \
    -o results/ \
    -r chr1:1000000-2000000 \
    --sele-var 2
```

### 高级用法 | Advanced Usage

```bash
# 结合GWAS和基因注释的可视化
biopytools ldblockshow \
    -i variants.vcf.gz \
    -o results/ \
    -r chr1:1000000-2000000 \
    --in-gwas gwas.txt \
    --in-gff annotation.gff \
    --out-png \
    --sele-var 4

# 亚群分析
biopytools ldblockshow \
    -i variants.vcf.gz \
    -o results/ \
    -r chr1:1000000-2000000 \
    --sub-pop subgroup_samples.txt \
    --maf 0.01 \
    --miss 0.2
```

## 参数说明 | Parameters

### 必需参数 | Required Parameters

| 参数 | 描述 | 示例 |
|------|------|------|
| `-i, --vcf-file` | VCF变异文件路径 | `-i variants.vcf.gz` |
| `-o, --output-dir` | 输出目录(自动创建)；每 region 产物落在 目录/<label>.* | `-o results/` |
| `-r, --region` | 单个分析区域，格式chr:start-end（与 `-b` 二选一） | `-r chr1:1000000-2000000` |
| `-b, --bed` | 基因组BED文件(每行 chrom start end [name])，等价多个 -r | `-b peaks.bed` |

### LD度量选择 | LD Statistic Selection

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--sele-var` | `1` | LD度量统计量：1=D', 2=R², 3=R²和D', 4=D'和R² |

### 过滤参数 | Filter Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--maf` | `0.05` | 最小次要等位基因频率 |
| `--miss` | `0.25` | 最大缺失率 |
| `--hwe` | `0.0` | Hardy-Weinberg平衡P值阈值 |
| `--het` | `1.0` | 最大杂合率 |

### Block检测参数 | Block Detection Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--block-type` | `1` | Block检测方法：1=Gabriel(PLINK), 2=SolidSpine, 3=BlockCut, 4=FixBlock, 5=NoBlock |
| `--block-cut` | `0.85:0.90` | BlockType3的cutoff（格式：cutoff:ratio） |
| `--fix-block` | - | 固定block文件路径（用于BlockType=4） |

### 可视化参数 | Visualization Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--in-gwas` | - | GWAS P值文件（格式：chr site Pvalue） |
| `--in-gff` | - | GFF3注释文件路径 |
| `--mer-min-snp-num` | `50` | 合并网格的最小SNP数 |

### 输出格式 | Output Format

| 参数 | 描述 |
|------|------|
| `--out-png` | 输出PNG格式图像 |
| `--out-pdf` | 输出PDF格式图像 |

### 其他参数 | Other Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--sub-pop` | - | 亚群样本文件路径 |
| `--tag-snp-cut` | `0.80` | TagSNP的LD cutoff |
| `--ldblockshow-path` | `/share/org/YZWL/yzwl_lixg/software/LDBlockShow/bin/LDBlockShow` | LDBlockShow可执行文件路径 |

## 输入文件格式 | Input File Formats

### VCF文件 | VCF File

标准VCF格式的变异文件，支持gzip压缩：

```vcf
##fileformat=VCFv4.2
##contig=<ID=chr1,length=249250621>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample1
chr1	10000	rs1234	A	G	45.2	PASS	DP=30	GT:DP	1/1:30
chr1	10100	rs5678	T	C	62.8	PASS	DP=25	GT:DP	0/1:25
```

### GWAS P值文件 | GWAS P-value File

三列文本文件，包含染色体、位点和P值：

```
chr1	1001000	3.5e-08
chr1	1005000	2.1e-05
chr1	1010000	1.2e-03
```

### GFF3注释文件 | GFF3 Annotation File

标准GFF3格式：

```gff3
##gff-version 3
chr1	RefSeq	gene	1000000	2000000	.	+	.	ID=gene1;Name=GENE1
chr1	RefSeq	mRNA	1000000	1500000	.	+	.	ID=transcript1;Parent=gene1
chr1	RefSeq	exon	1000000	1000500	.	+	.	ID=exon1;Parent=transcript1
chr1	RefSeq	CDS	1000100	1000400	.	+	0	ID=cds1;Parent=transcript1
```

### 亚群样本文件 | Subgroup Sample File

每行一个样本名：

```
sample1
sample2
sample3
```

## 输出结果 | Output Results

### 输出文件说明 | Output Files Description

| 文件 | 描述 |
|------|------|
| `{prefix}.svg` | LD热图SVG格式图像 |
| `{prefix}.png` | LD热图PNG格式图像（使用--out-png） |
| `{prefix}.pdf` | LD热图PDF格式图像（使用--out-pdf） |
| `{prefix}.site.gz` | 过滤后保留的SNP位点列表 |
| `{prefix}.blocks.gz` | 检测到的LD blocks |
| `{prefix}.TriangleV.gz` | 区域配对LD值（R²/D'） |

### 输出目录结构 | Output Directory Structure

```
results/
├── chr1_1000000_2000000.svg              # LD热图SVG图像
├── chr1_1000000_2000000.png              # LD热图PNG图像（可选）
├── chr1_1000000_2000000.pdf              # LD热图PDF图像（可选）
├── chr1_1000000_2000000.site.gz          # SNP位点文件
├── chr1_1000000_2000000.blocks.gz        # LD block文件
├── chr1_1000000_2000000.TriangleV.gz     # 配对LD值文件
└── ldblockshow_processing_*.log            # 分析日志
```

## 使用示例 | Usage Examples

### 示例1：基本LD热图分析

```bash
biopytools ldblockshow \
    -i /data/population.vcf.gz \
    -o /results/ \
    -r chr1:1000000-2000000 \
    --out-png
```

### 示例2：结合GWAS结果的可视化

```bash
biopytools ldblockshow \
    -i /data/population.vcf.gz \
    -o /results/ \
    -r chr5:25000000-26000000 \
    --in-gwas /data/gwas_results.txt \
    --out-png \
    --sele-var 4
```

### 示例3：添加基因注释信息

```bash
biopytools ldblockshow \
    -i /data/population.vcf.gz \
    -o /results/ \
    -r chr3:50000000-51000000 \
    --in-gff /data/annotation.gff \
    --in-gwas /data/gwas_results.txt \
    --out-png \
    --sele-var 4
```

### 示例4：亚群分析

```bash
# 创建亚群样本文件
cat > subgroup.txt << EOF
sample1
sample2
sample3
EOF

biopytools ldblockshow \
    -i /data/population.vcf.gz \
    -o /results/ \
    -r chr2:10000000-11000000 \
    --sub-pop subgroup.txt \
    --maf 0.01 \
    --out-png
```

### 示例5：自定义Block检测

```bash
biopytools ldblockshow \
    -i /data/population.vcf.gz \
    -o /results/ \
    -r chr4:30000000-31000000 \
    --block-type 3 \
    --block-cut 0.90:0.95 \
    --out-png
```

## 分析流程 | Analysis Pipeline

1. **VCF文件读取**: 读取并解析VCF格式变异文件
2. **质量控制**: 应用MAF、Missing、HWE、Het过滤
3. **LD计算**: 计算SNP配对之间的LD值（D'和/或R²）
4. **Block检测**: 根据选定方法检测LD blocks
5. **图像生成**: 生成LD热图SVG/PNG/PDF图像
6. **结果输出**: 输出SNP位点、blocks和配对LD值文件

## 系统要求 | System Requirements

### 依赖软件 | Dependencies

- **LDBlockShow** (版本 1.41+)
  - 下载地址: https://github.com/hewm2008/LDBlockShow
  - 已预装在: `/share/org/YZWL/yzwl_lixg/software/LDBlockShow/bin/`
- **Python** (版本 3.7+)
- **Python包**:
  - `click` - 命令行界面

### 软件安装 | Software Installation

LDBlockShow已预装在系统中。如需自行安装：

```bash
# 克隆仓库
git clone https://github.com/hewm2008/LDBlockShow.git
cd LDBlockShow
chmod 755 configure
./configure
make
mv LDBlockShow bin/
```

### 硬件建议 | Hardware Recommendations

- **CPU**: 多核处理器（推荐4核以上）
- **RAM**: 最少4GB（大样本推荐16GB以上）
- **存储**: 预留VCF文件大小5倍的磁盘空间

## 注意事项 | Important Notes

1. **区域格式**: Region参数格式必须为chr:start-end，使用冒号和短横线分隔
2. **VCF格式**: 建议使用gzip压缩的VCF文件以节省存储空间
3. **SNP数量**: 区域内SNP数量过多（>1000）会导致图像非常大，建议使用--mer-min-snp-num参数
4. **LD度量选择**:
   - D'适合小样本量和罕见变异
   - R²更常用，适合大多数情况
   - 选择3或4可同时计算两种度量
5. **Block类型**: Gabriel方法（BlockType=1）是最常用的方法

## 故障排除 | Troubleshooting

### 常见问题 | Common Issues

**Q: "LDBlockShow: command not found" 错误**
```bash
# 检查LDBlockShow路径
ls -la /share/org/YZWL/yzwl_lixg/software/LDBlockShow/bin/LDBlockShow

# 使用--ldblockshow-path指定正确路径
biopytools ldblockshow ... --ldblockshow-path /correct/path/to/LDBlockShow
```

**Q: VCF文件格式错误**
```bash
# 检查VCF文件格式
zcat variants.vcf.gz | head -20

# 确保有正确的header行（#CHROM开头的行）
```

**Q: 区域内SNP数量太少**
```bash
# 降低MAF阈值以包含更多SNP
biopytools ldblockshow ... --maf 0.01
```

**Q: 输出SVG文件无法打开**
```bash
# 使用--out-pdf生成PDF文件
biopytools ldblockshow ... --out-pdf

# 或使用--out-png生成PNG文件
biopytools ldblockshow ... --out-png
```

## 结果解读指南 | Result Interpretation Guide

### LD热图解读

- **红色区域**: 高LD (R²/D'接近1)
- **白色区域**: 无LD (R²/D'接近0)
- **渐变色**: 中等LD程度
- **黑色边框**: 检测到的LD blocks边界

### LD Block类型说明

1. **Gabriel方法 (BlockType=1)**: 基于置信区间的方法，最常用
2. **Solid spine (BlockType=2)**: 基于连续强LD的方法
3. **BlockCut (BlockType=3)**: 自定义cutoff的方法
4. **FixBlock (BlockType=4)**: 使用预定义的block区域
5. **NoBlock (BlockType=5)**: 不检测blocks

### GWAS可视化说明

- **散点图**: GWAS P值（-log10转换）
- **水平线**: 显著性阈值线（使用ShowLDSVG工具设置）
- **峰值点**: 最显著的SNP位点

## 引用信息 | Citation

如果在学术研究中使用此工具，请引用原始LDBlockShow文章：

```
Wang, H., et al. (2020). LDBlockShow: a fast and effective tool for
linkage disequilibrium (LD) heatmap and block visualization.
Briefings in Bioinformatics, 21(4), 1321-1327.

PMID: 33126273 | DOI: 10.1093/bib/bbaa227
```

## 相关资源 | Related Resources

- [LDBlockShow GitHub](https://github.com/hewm2008/LDBlockShow)
- [LDBlockShow Manual (English)](https://github.com/hewm2008/LDBlockShow/blob/main/LDBlockShow_Manual_English.pdf)
- [LDBlockShow Manual (Chinese)](https://github.com/hewm2008/LDBlockShow/blob/main/LDBlockShow_Manual_Chinese.pdf)

## 许可证 | License

本项目采用MIT许可证 - 详见 [LICENSE](LICENSE) 文件

**注意**: LDBlockShow软件本身遵循其原始许可证。

## 更新日志 | Changelog

### Version 1.0.0 (2026-02-09)

- 初始版本发布|Initial release
- 完整的LDBlockShow Python封装|Complete LDBlockShow Python wrapper
- CLI接口支持|CLI interface support
- 符合biopytools开发规范|Compliant with biopytools development standards
