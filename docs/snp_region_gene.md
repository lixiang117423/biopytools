# SNP区域基因提取工具

## 功能概述

根据SNP位置和指定区域，提取相关基因的CDS和蛋白序列。

## 主要特性

- 智能识别SNP文件格式（支持表头、Chr:Pos格式）
- 灵活的区域定义（上游、下游距离可配置）
- 自动识别SNP特征（启动子、外显子、内含子）
- 考虑正负链对启动子定义的影响
- 多转录本全部输出（保留可变剪切信息）
- 输出基因-序列对应关系列表

## 快速开始

### 基本用法

```bash
# 提取SNP上下游100kb内的基因
biopytools snp-region-gene \
    -i snp_region.txt \
    -g annotation.gff3 \
    -G genome.fa \
    -l 100000 \
    -r 100000 \
    -o output
```

### 高级用法

```bash
# 自定义启动子区域和工具路径
biopytools snp-region-gene \
    -i snp_region.txt \
    -g annotation.gff3 \
    -G genome.fa \
    -l 50000 \
    -r 50000 \
    -p 3000 \
    --gffread-path /path/to/gffread \
    --seqkit-path /path/to/seqkit \
    -o results
```

### 只提取SNP位置对应的基因

```bash
# left=0, right=0时只提取SNP落在基因范围内的序列
biopytools snp-region-gene \
    -i snp_region.txt \
    -g annotation.gff3 \
    -G genome.fa \
    -l 0 \
    -r 0 \
    -o snp_genes
```

## 参数说明

### 必需参数

| 参数 | 描述 | 示例 |
|------|------|------|
| `-i, --snp` | SNP位置文件（格式：Chr01:24770） | `snp_region.txt` |
| `-g, --gff` | GFF3注释文件 | `annotation.gff3` |
| `-G, --genome` | 基因组FASTA文件 | `genome.fa` |

### 可选参数

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-l, --left` | `0` | SNP上游距离（bp） |
| `-r, --right` | `0` | SNP下游距离（bp） |
| `-p, --promoter` | `2000` | 启动子区域距离（bp） |
| `-o, --output` | `./snp_region_output` | 输出文件前缀 |
| `--gffread-path` | `gffread` | gffread程序路径 |
| `--seqkit-path` | `seqkit` | seqkit程序路径 |
| `--keep-temp` | `False` | 保留临时文件 |

## 输入文件格式

### SNP位置文件

支持带表头的格式：

```
SNP
Chr01:24770
Chr01:24863
Chr01:36202
```

也支持无表头格式：

```
Chr01:24770
Chr01:24863
Chr01:36202
```

**说明**：
- 自动跳过以#开头的注释行
- 自动识别包含"SNP"、"POSITION"等关键词的表头行
- 格式：`染色体:坐标`，用冒号分隔
- 自动去重

### GFF3注释文件

标准GFF3格式，包含gene、mRNA、exon、CDS特征。

## 输出文件

### 1. {prefix}_cds.fasta

区域内基因的CDS序列。

示例：
```
>Glsoj01G0000100.1
ATGGCT...
```

### 2. {prefix}_protein.fasta

区域内基因的蛋白序列。

示例：
```
>Glsoj01G0000100.1
MAV...
```

### 3. {prefix}_gene_list.txt

SNP-基因对应关系及特征信息。

列说明：
- `SNP_Position`: SNP标识（Chr01:24770）
- `Chromosome`: 染色体
- `SNP_Pos`: SNP坐标
- `Gene_ID`: 基因ID
- `mRNA_ID`: 转录本ID
- `Strand`: 正负链（+/-）
- `Feature`: SNP特征（promoter/exon/intron/intergenic）
- `Distance`: 距离值
  - 正链：SNP坐标 - gene_start
  - 负链：gene_end - SNP坐标
  - 负值表示在启动子区

示例：
```
SNP_Position	Chromosome	SNP_Pos	Gene_ID	mRNA_ID	Strand	Feature	Distance
Chr01:24770	Chr01	24770	Glsoj01G0000100	Glsoj01G0000100.1	+	exon	1234
Chr01:24770	Chr01	24770	Glsoj01G0000200	Glsoj01G0000200.1	-	intron	-567
Chr01:24863	Chr01	24863	Glsoj01G0000300	Glsoj01G0000300.1	+	promoter	-1500
```

## SNP特征判断

### 优先级

启动子 > 外显子 > 内含子 > 基因间区

### 启动子定义

- **正链（+）**: gene_start - promoter_length 到 gene_start
- **负链（-）**: gene_end 到 gene_end + promoter_length

### 基因匹配规则

1. **left=0, right=0**：
   - SNP落在基因的 [start, end] 范围内（包含启动子）就算该基因

2. **left>0 或 right>0**：
   - 基因与 [SNP-left, SNP+right] 区间有重叠就算

## 处理流程

```
步骤1: 解析SNP文件（智能识别表头）
步骤2: 解析GFF3，建立基因索引
步骤3: 使用gffread提取所有CDS和蛋白序列
步骤4: 对每个SNP查找相关基因
步骤5: 使用seqkit提取指定mRNA的序列
步骤6: 输出结果文件
步骤7: 清理临时文件
```

## 系统要求

### 依赖软件

- **gffread**（Cufflinks套件）
  - 用于从GFF3和基因组提取CDS和蛋白序列
  - 下载：http://ccb.jhu.edu/software/stringtie/gff

- **seqkit**
  - 用于根据ID快速提取序列
  - 下载：https://github.com/shenwei356/seqkit

### Python依赖

- Python 3.7+
- click（命令行接口）

## 使用示例

### 示例1：提取SNP上下游基因

```bash
biopytools snp-region-gene \
    -i snp_list.txt \
    -g genome.gff3 \
    -G genome.fa \
    -l 100000 \
    -r 100000 \
    -o region_genes
```

输出：
- `region_genes_cds.fasta`
- `region_genes_protein.fasta`
- `region_genes_gene_list.txt`
- `region_genes.log`

### 示例2：只提取SNP位置基因

```bash
biopytools snp-region-gene \
    -i snp_list.txt \
    -g genome.gff3 \
    -G genome.fa \
    -o snp_exact_genes
```

### 示例3：自定义启动子区域

```bash
# 定义启动子为上游3kb
biopytools snp-region-gene \
    -i snp_list.txt \
    -g genome.gff3 \
    -G genome.fa \
    -l 50000 \
    -r 50000 \
    -p 3000 \
    -o custom_promoter
```

## 注意事项

1. **文件格式**：
   - 确保GFF3文件格式正确
   - 基因组FASTA文件与GFF3文件染色体名称一致

2. **工具路径**：
   - 确保gffread和seqkit在系统PATH中
   - 或使用--gffread-path和--seqkit-path指定路径

3. **内存使用**：
   - 大基因组可能需要较多内存
   - 临时文件大小约等于全基因组CDS和蛋白序列

4. **特征判断**：
   - 启动子定义考虑了正负链
   - 启动子重叠优先于外显子/内含子

## 故障排除

### 问题1：gffread命令未找到

```bash
# 解决方法1：安装gffread
# Ubuntu/Debian
sudo apt-get install gffread

# 解决方法2：指定路径
biopytools snp-region-gene ... --gffread-path /path/to/gffread
```

### 问题2：seqkit命令未找到

```bash
# 安装seqkit
# macOS
brew install seqkit

# Linux
wget https://github.com/shenwei356/seqkit/releases/download/v2.8.0/seqkit_linux_amd64.tar.gz
tar -xzf seqkit_linux_amd64.tar.gz
sudo mv seqkit /usr/local/bin/

# 或指定路径
biopytools snp-region-gene ... --seqkit-path /path/to/seqkit
```

### 问题3：没有找到任何基因

可能原因：
- SNP文件格式不正确
- left和right距离太小
- 染色体名称不匹配

检查方法：
```bash
# 检查染色体名称
grep "^>" genome.fa | head
grep -v "^#" genome.gff3 | head -5
cut -d: -f1 snp_list.txt | sort -u
```

## 版本历史

| 版本 | 日期 | 说明 |
|------|------|------|
| 1.0.0 | 2026-01-09 | 初始版本 |
