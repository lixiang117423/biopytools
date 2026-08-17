# 🧬 VCF转PHYLIP格式转换模块

**VCF变异数据到系统发育分析格式的专业转换工具 | Professional VCF to Phylogenetic Format Conversion Tool**

## 📖 功能概述 | Overview

VCF转PHYLIP格式转换模块是一个专业的生物信息学工具，用于将VCF格式的SNP数据转换为多种系统发育分析所需的格式。支持PHYLIP、FASTA、NEXUS和二进制NEXUS格式，具备高效的样本和位点过滤功能，是群体基因组学和系统发育研究的重要预处理工具。

## ✨ 主要特性 | Key Features

- **🔄 多格式输出**: 支持PHYLIP、FASTA、NEXUS、二进制NEXUS四种主流格式
- **⚡ 高效处理**: 多线程并行处理，支持大规模VCF文件快速转换
- **🎯 智能过滤**: 灵活的样本覆盖度过滤和质量控制
- **🧬 IUPAC处理**: 支持杂合子基因型的IUPAC模糊代码处理和随机解析
- **🌳 外群支持**: 可指定外群样本并自动排列在输出矩阵首位
- **📊 位点追踪**: 可输出通过筛选的位点坐标列表，便于质量控制
- **🗂️ 文件管理**: 智能的输出目录管理和文件命名策略
- **📝 详细日志**: 完整的处理过程记录和错误追踪

## 🚀 快速开始 | Quick Start

### 基本用法 | Basic Usage

```bash
# 基本VCF转换（默认PHYLIP格式）
biopytools vcf2phylip -i population_snps.vcf -o phylo_input

# 转换为多种格式
biopytools vcf2phylip -i variants.vcf.gz -o multi_format \
    --fasta --nexus --threads 32

# 指定外群和最小样本数
biopytools vcf2phylip -i cohort.vcf -o analysis \
    -m 8 -g outgroup_sample --fasta
```

### 高级用法 | Advanced Usage

```bash
# 完整格式转换分析
biopytools vcf2phylip \
    -i population_variants.vcf.gz \
    -o comprehensive_analysis \
    --output-prefix pop_study \
    --fasta \
    --nexus \
    --nexus-binary \
    -m 10 \
    -g ancestral_sample \
    --resolve-IUPAC \
    --write-used-sites \
    --threads 64

# 大规模数据处理
biopytools vcf2phylip \
    -i large_cohort.vcf.gz \
    -o batch_conversion \
    --phylip-disable \
    --fasta \
    -m 20 \
    --threads 88
```

## 📋 命令行参数 | Command Line Parameters

### 必需参数 | Required Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-i, --input` | `None` | 📁 输入VCF文件路径（支持gzip压缩） |

### 输出配置参数 | Output Configuration

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-o, --output` | `./converted_output` | 📂 输出目录路径 |
| `--output-prefix` | `None` | 📝 输出文件名前缀 |

### 过滤参数 | Filtering Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-m, --min-samples-locus` | `4` | 🎯 每个位点最少样本数 |
| `-g, --outgroup` | `""` | 🌳 外群样本名称 |

### 输出格式参数 | Output Format Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-p, --phylip-disable` | `False` | ⏭️ 禁用PHYLIP输出 |
| `-f, --fasta` | `False` | ✅ 启用FASTA输出 |
| `-n, --nexus` | `False` | ✅ 启用NEXUS输出 |
| `-b, --nexus-binary` | `False` | ✅ 启用二进制NEXUS输出 |

### 处理选项 | Processing Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-r, --resolve-IUPAC` | `False` | 🎲 随机解析杂合子基因型 |
| `-w, --write-used-sites` | `False` | 📝 保存通过筛选的位点坐标 |
| `-t, --threads` | `88` | 🧵 并行线程数 |

### 步骤说明 | Step Descriptions

| 步骤 | 名称 | 描述 |
|------|------|------|
| **1** | 📖 VCF解析 | 验证VCF格式和提取样本信息 |
| **2** | 🔍 质量过滤 | 基于样本覆盖度过滤SNP位点 |
| **3** | 🧬 基因型处理 | 处理杂合子和缺失数据 |
| **4** | 📄 格式转换 | 生成多种系统发育分析格式 |
| **5** | 🌳 外群处理 | 外群样本识别和重新排序 |
| **6** | 📊 结果输出 | 多格式文件生成和质量报告 |

## 📁 输入文件格式 | Input File Formats

### VCF变异文件 | VCF Variant File

标准VCF格式的变异调用结果：

```vcf
##fileformat=VCFv4.2
##contig=<ID=chr1,length=249250621>
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  Sample1  Sample2  Sample3
chr1    10177   rs367896724     A       AC      100     PASS    AC=1;AF=0.500;AN=2;DP=13 GT      0/1     1/1     0/0
chr1    10352   rs555500075     T       TA      100     PASS    AC=1;AF=0.500;AN=2;DP=10 GT      0/1     0/0     0/1
```

**文件要求**:
- 符合VCF 4.0+标准格式
- 支持gzip压缩格式（.vcf.gz）
- 包含完整的基因型信息（GT字段）
- 样本名称必须唯一且不含特殊字符

**数据质量要求**:
- SNP位点应经过基本质量控制
- 建议预先进行MAF和缺失率过滤
- 避免过多的低质量基因型调用
- 样本覆盖度相对均匀

## 📂 输出文件说明 | Output Files Description

### PHYLIP格式 | PHYLIP Format

标准PHYLIP序列格式：

```phylip
   10  1250
Sample1     ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG...
Sample2     ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG...
Outgroup    ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG...
```

**特点**:
- 适用于RAxML、IQ-TREE、MEGA等软件
- 严格格式限制（序列名≤10字符）或relaxed格式
- 包含样本数和位点数信息头
- 外群样本自动排在首位

### FASTA格式 | FASTA Format

标准FASTA序列格式：

```fasta
>Sample1
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG...
>Sample2
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG...
>Outgroup
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG...
```

**特点**:
- 通用性强，兼容大多数分析软件
- 支持长序列名称
- 易于后续处理和编辑
- 样本名作为序列标识符

### NEXUS格式 | NEXUS Format

NEXUS数据块格式：

```nexus
#NEXUS
BEGIN DATA;
    DIMENSIONS NTAX=10 NCHAR=1250;
    FORMAT DATATYPE=DNA MISSING=? GAP=-;
    MATRIX
Sample1    ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG...
Sample2    ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG...
Outgroup   ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG...
    ;
END;
```

**特点**:
- 适用于MrBayes、BEAST等贝叶斯分析软件
- 包含完整的数据块定义
- 支持字符状态和缺失数据处理
- 可包含外群和分组信息

### 二进制NEXUS格式 | Binary NEXUS Format

二进制编码的NEXUS格式：

```nexus
#NEXUS
BEGIN DATA;
    DIMENSIONS NTAX=10 NCHAR=1250;
    FORMAT DATATYPE=STANDARD MISSING=? GAP=- SYMBOLS="01";
    MATRIX
Sample1    010011010010011010100110101001011001011...
Sample2    010011010010011010100110101001011001011...
Outgroup   010011010010011010100110101001011001011...
    ;
END;
```

**特点**:
- 仅适用于二倍体数据
- SNP位点编码为0/1二进制
- 显著减少文件大小
- 适用于特定的贝叶斯分析

### 位点坐标文件 | Sites Coordinate File

通过筛选的位点坐标列表：

```
chr1    10177
chr1    10352
chr2    54678
chr2    89012
...
```

## 🔧 使用示例 | Usage Examples

### 群体遗传学分析

```bash
# 人类群体SNP数据转换
biopytools vcf2phylip \
    -i human_population.vcf.gz \
    -o human_phylo \
    --fasta \
    --nexus \
    -m 50 \
    -g chimp_outgroup \
    --threads 64

# 输出文件：
# - human_phylo.phy (PHYLIP格式)
# - human_phylo.fas (FASTA格式)
# - human_phylo.nex (NEXUS格式)
# - human_phylo_sites.txt (位点坐标)
```

### 植物系统发育研究

```bash
# 植物种质资源SNP转换
biopytools vcf2phylip \
    -i plant_gwas.vcf \
    -o plant_analysis \
    --output-prefix rice_diversity \
    --phylip-disable \
    --fasta \
    --nexus-binary \
    -m 8 \
    -g wild_species \
    --resolve-IUPAC \
    --write-used-sites \
    --threads 32
```

### 微生物群体分析

```bash
# 细菌菌株SNP核心基因组转换
biopytools vcf2phylip \
    -i bacterial_core_snps.vcf.gz \
    -o bacterial_phylo \
    --fasta \
    --nexus \
    -m 20 \
    -g reference_strain \
    --threads 88
```

## ⚡ 性能优化 | Performance Optimization

### 线程配置建议

| 文件大小 | 样本数 | 推荐线程数 | 内存需求 |
|----------|--------|------------|----------|
| <10MB | <50 | 4-8 | <2GB |
| 10MB-1GB | 50-200 | 16-32 | 2-8GB |
| >1GB | >200 | 32-88 | 8-32GB |

### 内存优化策略

```bash
# 大文件分批处理
biopytools vcf2phylip \
    -i huge_dataset.vcf.gz \
    -o batch_results \
    --threads 32 \
    -m 15

# 使用SSD存储提升I/O性能
export TMPDIR=/ssd/temp
biopytools vcf2phylip -i data.vcf -o results --threads 64
```

### 输出格式选择优化

```bash
# 最小文件体积（仅二进制NEXUS）
biopytools vcf2phylip -i data.vcf -o minimal --nexus-binary

# 最大兼容性（多格式输出）
biopytools vcf2phylip -i data.vcf -o compatible --fasta --nexus

# 快速处理（仅FASTA）
biopytools vcf2phylip -i data.vcf -o fast --phylip-disable --fasta
```

## 🐛 故障排除 | Troubleshooting

### 常见问题 | Common Issues

**内存不足错误**
```bash
# 解决方案1：减少线程数
biopytools vcf2phylip -i large_file.vcf -o output --threads 16

# 解决方案2：提高最小样本数阈值
biopytools vcf2phylip -i large_file.vcf -o output -m 20

# 解决方案3：仅输出必要格式
biopytools vcf2phylip -i large_file.vcf -o output --fasta
```

**VCF格式错误**
```bash
# 预检查VCF文件
bcftools view -H input.vcf | head -5

# 验证样本数量
bcftools query -l input.vcf | wc -l

# 检查基因型字段
bcftools view -H input.vcf | cut -f10 | head -5
```

**外群样本问题**
```bash
# 检查外群样本是否存在
bcftools query -l input.vcf | grep outgroup_name

# 验证样本名拼写
bcftools query -l input.vcf
```

### 性能调试 | Performance Debugging

```bash
# 监控内存使用
htop

# 检查磁盘空间
df -h /tmp

# 监控I/O性能
iotop

# 测试线程性能
for t in 4 8 16 32 64; do
    time biopytools vcf2phylip -i test.vcf -o test_${t} --threads ${t}
done
```

## 📊 质量控制建议 | Quality Control Recommendations

### 输入数据预处理

```bash
# 使用bcftools进行基本过滤
bcftools view -i 'QUAL>30 && MAF>0.05 && INFO/DP>10' input.vcf -o filtered.vcf

# 移除高缺失率样本
vcftools --vcf input.vcf --missing-indv
# 查看missing报告，移除缺失率>20%的样本

# 移除低质量位点
vcftools --vcf input.vcf --missing-site
# 查看位点missing报告，过滤高缺失位点
```

### 转换参数优化

```bash
# 保守参数设置（高质量数据）
biopytools vcf2phylip \
    -i high_quality.vcf \
    -o conservative \
    -m 15 \
    --write-used-sites

# 宽松参数设置（数据稀疏时）
biopytools vcf2phylip \
    -i sparse_data.vcf \
    -o lenient \
    -m 3 \
    --resolve-IUPAC
```

### 结果验证

```bash
# 检查输出文件格式
head -5 output.phy
head -5 output.fas
head -10 output.nex

# 验证序列长度一致性
grep -v "^>" output.fas | awk '{print length($0)}' | sort | uniq

# 统计转换信息
grep -c "^>" output.fas  # 样本数
head -2 output.fas | tail -1 | wc -c  # 序列长度
```

## 🔗 相关文档 | Related Documentation

- [VCF规范文档](https://samtools.github.io/hts-specs/VCFv4.3.pdf)
- [PHYLIP格式说明](http://evolution.genetics.washington.edu/phylip/phylipdoc.html)
- [NEXUS格式规范](https://wiki.chaos.org.uk/NEXUS_format)
- [IUPAC核苷酸代码表](https://www.bioinformatics.org/sms/iupac.html)
- [biopytools系统发育分析](vcf_phylo.md)
- [IQ-TREE系统发育软件](https://iqtree.org/)
- [RAxML软件手册](https://cme.h-its.org/exelixis/web/software/raxml/)

## 📄 许可证 | License

本模块遵循MIT许可证。详细信息请参见LICENSE文件。

## 🤝 贡献指南 | Contributing

欢迎提交Issue和Pull Request来改进本模块。

## 📞 技术支持 | Support

如有技术问题，请联系：
- 邮箱: yzwl_lixg@outlook.com
- 项目地址: https://github.com/your-org/biopytools

---

**最后更新**: 2024年12月17日
**版本**: 2.9.1
**作者**: biopytools开发团队

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input, -i` | 必填 | Path | 输入VCF文件路径｜Input VCF file path |
| `--output, -o` | `./converted_output` | Path | 输出目录｜Output directory |
| `--output-prefix` | — | str | 输出文件名前缀｜Output filename prefix |
| `--min-samples-locus, -m` | `4` | int | 位点最少样本数｜Minimum samples per locus |
| `--outgroup, -g` | `` | str | 外群样本名称｜Outgroup sample name |
| `--phylip-disable, -p` | — |  | 禁用PHYLIP输出｜Disable PHYLIP output |
| `--fasta, -f` | — |  | 启用FASTA输出｜Enable FASTA output |
| `--nexus, -n` | — |  | 启用NEXUS输出｜Enable NEXUS output |
| `--nexus-binary, -b` | — |  | 启用二进制NEXUS输出｜Enable binary NEXUS output |
| `--resolve-IUPAC, -r` | — |  | 随机解析杂合子基因型｜Resolve heterozygous genotypes |
| `--write-used-sites, -w` | — |  | 保存筛选通过的位点坐标｜Save used sites coordinates |
| `--threads, -t` | `12` | int | 线程数｜Number of threads |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 输入VCF文件路径｜Input VCF file path |
| `-o, --output` | `./converted_output` |  | 输出目录｜Output directory |
| `--output-prefix` | — |  | 输出文件名前缀｜Output filename prefix |
| `-m, --min-samples-locus` | `4` | int | 位点最少样本数｜Minimum samples per locus |
| `-g, --outgroup` | `` |  | 外群样本名称｜Outgroup sample name |
| `-p, --phylip-disable` | — | store_true | 禁用PHYLIP输出｜Disable PHYLIP output |
| `-f, --fasta` | — | store_true | 启用FASTA输出｜Enable FASTA output |
| `-n, --nexus` | — | store_true | 启用NEXUS输出｜Enable NEXUS output |
| `-b, --nexus-binary` | — | store_true | 启用二进制NEXUS输出｜Enable binary NEXUS output |
| `-r, --resolve-IUPAC` | — | store_true | 随机解析杂合子基因型｜Resolve heterozygous genotypes |
| `-w, --write-used-sites` | — | store_true | 保存筛选通过的位点坐标｜Save used sites coordinates |
| `-t, --threads` | `88` | int | 线程数｜Number of threads |
| `-v, --version` | — | version |  |

<!-- END PARAMS:auto -->
