# 🧬 PopLD 连锁不平衡衰减分析模块

**基于PopLDdecay的高效LD衰减分析工具 | Efficient LD Decay Analysis Tool Based on PopLDdecay**

## 📖 功能概述 | Overview

PopLD模块是基于PopLDdecay软件的连锁不平衡（LD）衰减分析工具，支持单群体和多群体LD衰减分析。模块能自动识别样本分组文件格式，按分组进行批量分析，并生成合并结果文件，适用于各种群体遗传学研究和进化分析。

## ✨ 主要特性 | Key Features

- **🎯 自动分组识别**: 智能检测样本文件格式（单列/两列），自动切换分析模式
- **📋 批量分组分析**: 一次性分析多个亚群，自动为每个分组生成独立结果
- **🔀 结果自动合并**: 将所有分组结果合并为一个文件，添加分组标签列
- **⚡ 高性能计算**: 支持两种算法模式，灵活平衡速度和内存
- **📊 灵活质控**: MAF、杂合率、缺失率等多维质量控制
- **🧬 多种输出格式**: 支持8种输出类型，满足不同分析需求
- **🔧 EHH分析**: 支持扩展单倍型纯合度（EHH）分析

## 🚀 快速开始 | Quick Start

### 基本用法 | Basic Usage

```bash
# 单群体分析（不指定样本文件）
biopytools popld -i snp.vcf.gz -o LD_result

# 多群体分组分析（两列样本文件）
biopytools popld -i snp.vcf.gz -s groups.txt -o LD_result
```

### 高级用法 | Advanced Usage

```bash
# 自定义质控参数
biopytools popld -i snp.vcf.gz -o LD_result \
    --maf 0.01 --het 0.9 --miss 0.1

# 使用算法2（可能需要更多内存）
biopytools popld -i snp.vcf.gz -o LD_result --method 2

# 更大距离范围分析
biopytools popld -i snp.vcf.gz -o LD_result -d 500
```

## 📋 参数说明 | Parameters

### 必需参数 | Required Parameters

| 参数 | 描述 | 示例 |
|------|------|------|
| `-i, --input-vcf` | 输入VCF文件路径（支持压缩） | `-i snp.vcf.gz` |
| `-o, --output-stat` | 输出统计文件前缀 | `-o LD_result` |

### 软件配置 | Software Configuration

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-p, --poplddecay-path` | `/share/org/YZWL/yzwl_lixg/.../PopLDdecay` | PopLDdecay软件路径 |

### 分组样本配置 | Sample Group Configuration

| 参数 | 描述 | 示例 |
|------|------|------|
| `-s, --sub-pop` | 样本列表文件（支持两种格式） | `-s groups.txt` |

**样本文件格式说明**：

**格式1：单列格式**（单群体分析）
```text
Sample1
Sample2
Sample3
```

**格式2：两列格式**（多群体分组分析）
```text
Sample1	GroupA
Sample2	GroupA
Sample3	GroupB
Sample4	GroupB
```

### 处理参数 | Processing Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-d, --max-dist` | `300` | SNP间最大距离(kb) |
| `-m, --maf` | `0.005` | 最小次等位基因频率 |
| `--het` | `0.88` | 最大杂合率 |
| `--miss` | `0.25` | 最大缺失率 |
| `--method` | `1` | 算法方法 (1=快速, 2=可能需要更多内存) |
| `--out-type` | `1` | 输出类型 (1-8) |
| `--ehh` | `None` | EHH起始位点 (格式: chr:position) |

### 其他选项 | Other Options

| 参数 | 描述 |
|------|------|
| `--out-filter-snp` | 输出最终用于计算的SNP |

## 📁 输入文件格式 | Input File Formats

### VCF文件要求 | VCF File Requirements

支持标准VCF格式文件（压缩或未压缩）：

```vcf
##fileformat=VCFv4.2
##contig=<ID=1,length=249250621>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample1	Sample2	Sample3
1	14370	rs6054257	G	A	29	PASS	.	GT	0/0	1/0	1/1
1	17330	.	T	A	3	q10	.	GT	0/0	0/1	0/0
```

**文件要求**:
- 包含完整的基因型信息（GT字段）
- 至少包含2个样本
- 支持gzip压缩格式（`.vcf.gz`）

## 📂 输出文件结构 | Output File Structure

### 单群体分析输出 | Single Population Output

```
LD_result.stat.gz    # LD衰减统计结果
popld.log            # 分析日志文件
```

### 多群体分组分析输出 | Multi-Group Population Output

```
LD_result.GroupA.stat.gz      # GroupA分组的统计结果
LD_result.GroupB.stat.gz      # GroupB分组的统计结果
LD_result.GroupC.stat.gz      # GroupC分组的统计结果
...
LD_result.groups.list         # 多群体绘图列表文件
LD_result.merged.stat.gz      # 合并结果文件（含分组标签）
popld.log                     # 分析日志文件
```

## 📊 输出文件格式 | Output File Formats

### 单群体统计文件 (*.stat.gz)

```text
#Dist	Mean_r^2	Mean_D'	Sum_r^2	Sum_D'	NumberPairs
1	0.7550	NA	13974.3413	NA	18510
2	0.7249	NA	10588.0531	NA	14607
3	0.7342	NA	11508.9986	NA	15676
```

### 合并结果文件 (*.merged.stat.gz)

```text
#Dist	Mean_r^2	Mean_D'	Sum_r^2	Sum_D'	NumberPairs	Group
1	0.7550	NA	13974.3413	NA	18510	GroupA
1	0.7234	NA	12345.6789	NA	17080	GroupB
2	0.7249	NA	10588.0531	NA	14607	GroupA
```

**列说明**：
- `Dist`: SNP对之间的距离（kb）
- `Mean_r^2`: 该距离下所有SNP对的r²平均值
- `Mean_D'`: 该距离下所有SNP对的D'平均值
- `Sum_r^2`: 该距离下所有SNP对的r²总和
- `Sum_D'`: 该距离下所有SNP对的D'总和
- `NumberPairs`: 该距离下的SNP对数量
- `Group`: 分组名称（仅合并文件）

## 💡 使用示例 | Usage Examples

### 示例1：基础单群体分析 | Example 1: Basic Single Population Analysis

```bash
# 分析整个VCF文件的LD衰减
biopytools popld \
    -i population_snps.vcf.gz \
    -o population_LD
```

### 示例2：亚群分组分析 | Example 2: Subgroup Analysis

```bash
# 按群体分组进行LD分析
biopytools popld \
    -i population_snps.vcf.gz \
    -s sample_groups.txt \
    -o subgroup_LD
```

**sample_groups.txt 格式**：
```text
Wild_01	Wild
Wild_02	Wild
Cultivar_01	Cultivated
Cultivar_02	Cultivated
```

### 示例3：严格质控分析 | Example 3: Strict QC Analysis

```bash
# 使用严格的质控参数
biopytools popld \
    -i population_snps.vcf.gz \
    -s groups.txt \
    -o strict_LD \
    --maf 0.05 \
    --het 0.9 \
    --miss 0.1
```

### 示例4：大距离范围分析 | Example 4: Large Distance Range Analysis

```bash
# 分析更大距离范围的LD衰减
biopytools popld \
    -i population_snps.vcf.gz \
    -s groups.txt \
    -o long_range_LD \
    -d 1000
```

### 示例5：EHH分析 | Example 5: EHH Analysis

```bash
# 对特定位点进行EHH分析
biopytools popld \
    -i population_snps.vcf.gz \
    -o ehh_result \
    --ehh chr1:5000000
```

## 🎨 绘图说明 | Plotting Instructions

### 单群体绘图 | Single Population Plotting

分析完成后，使用PopLDdecay自带的Perl脚本绘图：

```bash
# 设置PopLDdecay路径
export POPLDPATH=/path/to/PopLDdecay

# 绘制单群体LD衰减图
perl $POPLDPATH/bin/Plot_OnePop.pl \
    -inFile population_LD.stat.gz \
    -output population_LD
```

**输出文件**：
- `population_LD.png` - PNG格式图形
- `population_LD.pdf` - PDF格式图形

### 多群体分组绘图 | Multi-Group Plotting

```bash
# 方法1：使用合并结果绘图（需要自定义脚本）
# 合并结果文件已包含Group列，可使用R/Python自定义绘图

# 方法2：使用PopLDdecay的多群体绘图
perl $POPLDPATH/bin/Plot_MutiPop.pl \
    -inList subgroup_LD.groups.list \
    -output subgroup_LD
```

### R语言绘图示例 | R Plotting Example

```r
# 读取合并结果
library(ggplot2)
data <- read.table("LD_result.merged.stat.gz", header=T, sep="\t")

# 绘制LD衰减曲线
ggplot(data, aes(x=Dist, y=Mean_r2, color=Group)) +
    geom_line(size=1) +
    geom_point(size=1) +
    labs(x="Distance (kb)", y="Mean r²", title="LD Decay Analysis") +
    theme_bw()
```

## 🔧 系统要求 | System Requirements

### 依赖软件 | Dependencies

- **PopLDdecay** (版本 3.43+)
  - 下载地址: https://github.com/hewm2008/PopLDdecay
- **Perl** (用于绘图脚本)
- **Python** (版本 3.7+)
- **Python包**:
  - `click` - 命令行界面

### 安装PopLDdecay | Installing PopLDdecay

```bash
# 方法1：使用Git克隆
git clone https://github.com/hewm2008/PopLDdecay.git
cd PopLDdecay
chmod 755 configure
./configure
make
mv PopLDdecay bin/

# 方法2：下载预编译版本
wget https://github.com/hewm2008/PopLDdecay/archive/v3.43.tar.gz
tar -xzf v3.43.tar.gz
cd PopLDdecay-3.43/src
make
```

## 📊 OutType参数详解 | OutType Parameter Details

PopLDdecay支持8种输出类型：

| OutType | 描述 | 输出内容 |
|---------|------|----------|
| 1 | 最快 | 仅 r² 距离统计 |
| 2 | 标准 | r² 和 D' 距离统计 |
| 3 | 包含成对 | r² + D' + 成对LD结果 |
| 4 | 多种统计 | r², D', 数量统计 |
| 5 | r²详细 | r² + 数量统计 |
| 6 | 成对详细 | r² + D' + 成对LD |
| 7 | 成对扩展 | r² + D' + LOD + 成对LD |
| 8 | 最完整 | r² + D' + LOD + CI + 成对LD |

**建议**：
- 日常分析使用：`--out-type 1`（最快）
- 需要D'值：`--out-type 2`
- 详细分析：`--out-type 4`

## ⚠️ 注意事项 | Important Notes

### 1. 样本文件格式

- 单列格式：每个样本一行，进行单群体分析
- 两列格式：样本ID + 分组名（Tab或空格分隔），自动进行多群体分组分析

### 2. 输出类型选择

- OutType越大，输出文件越大，计算时间越长
- 对于大样本数据集，建议使用OutType 1或2

### 3. 内存使用

- Method 1：默认算法，内存占用小
- Method 2：可能需要更多内存，适用于某些特殊场景

### 4. 质控参数

- MAF过滤掉稀有变异，通常使用0.01-0.05
- Het参数控制杂合子比例，自交物种建议调低（0.5-0.7）
- Miss参数控制缺失率，建议0.1-0.2

## 🐛 故障排除 | Troubleshooting

### 常见问题 | Common Issues

**Q: "PopLDdecay: command not found" 错误**
```bash
# 检查PopLDdecay安装
which PopLDdecay

# 指定正确路径
biopytools popld -i input.vcf -o result \
    -p /path/to/PopLDdecay
```

**Q: 样本文件格式问题**
```bash
# 确保样本文件使用UTF-8编码
file sample_groups.txt

# 检查是否包含特殊字符
cat -A sample_groups.txt | head -5
```

**Q: VCF文件样本名不匹配**
```bash
# 检查VCF文件中的样本名
zcat input.vcf.gz | grep "^#CHROM"

# 确保样本ID完全一致（区分大小写）
```

**Q: 内存不足**
```bash
# 减少MaxDist参数
biopytools popld -i input.vcf -o result -d 100

# 或增加MAF阈值减少SNP数量
biopytools popld -i input.vcf -o result --maf 0.05
```

**Q: 分组分析结果为空**
```bash
# 检查样本列表文件格式
head -5 sample_groups.txt

# 确保使用两列格式（Tab分隔）
awk '{print NF}' sample_groups.txt | sort -u
# 应该只有数字2
```

## 📚 参考文献 | References

1. **PopLDdecay软件**: Zhang et al. (2018) PopLDdecay: a fast and effective tool for linkage disequilibrium decay analysis based on variant call format files. *Bioinformatics* 35: 3658-3664.

2. **软件主页**: https://github.com/hewm2008/PopLDdecay

3. **LD衰减分析原理**: Slatkin M. (2008) Linkage disequilibrium - understanding the evolutionary past and selecting the most valuable populations for conservation. *Heridity* 100: 421-423.

## 📄 许可证 | License

本项目采用MIT许可证 - 详见 [LICENSE](../LICENSE) 文件

PopLDdecay软件遵循其原始许可证条款。

---

## 🔬 引用信息 | Citation

如果在学术研究中使用此工具，请引用：

```
Zhang, C., Dong, S., Xu, J., et al. (2019)
PopLDdecay: a fast and effective tool for linkage disequilibrium decay analysis based on variant call format files.
Bioinformatics, 35(5): 875-877.
DOI: 10.1093/bioinformatics/bty875
```
