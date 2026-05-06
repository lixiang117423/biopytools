# GFF文件ID规范化工具

## 功能概述

`biopytools gff-renamer` 是一个用于重命名GFF/GFF3文件中基因、转录本及子特征ID的工具。它将来自不同注释流程（NCBI Gnomon、EGAPx、BRAKER、EviAnn等）的非标准ID统一为规范格式，并输出新旧ID的映射关系。

## 核心功能

- 🔄 **多源格式兼容** - 支持NCBI Gnomon、EGAPx、BRAKER、EviAnn等主流注释流程的ID格式
- 🧹 **AGAT清洗集成** - 默认使用AGAT清洗GFF文件，自动补充UTR等缺失特征
- 📊 **多命名格式** - 提供standard、simple、compact三种命名格式
- 🔀 **并行处理** - 支持多线程加速大文件处理
- 📝 **映射文件输出** - 可选生成mRNA新旧ID映射文件
- 🧬 **染色体映射** - 支持自定义染色体名称映射

## 安装依赖

本工具默认启用AGAT清洗，需要预先安装AGAT：

```bash
# 使用conda安装
conda install -c bioconda agat
```

## 使用方法

### 基本用法

```bash
biopytools gff-renamer \
    -i complete.genomic.gff \
    -o output.gff \
    -p CDRT \
    -s Ii
```

### 包含UTR特征

```bash
biopytools gff-renamer \
    -i complete.genomic.gff \
    -o output.gff \
    -p CDRT \
    -s Ii \
    --include-utr
```

### 输出mRNA映射文件

```bash
biopytools gff-renamer \
    -i complete.genomic.gff \
    -o output.gff \
    -p CDRT \
    -s Ii \
    --output-mrna-mapping
```

### 使用自定义染色体映射

```bash
biopytools gff-renamer \
    -i complete.genomic.gff \
    -o output.gff \
    -p CDRT \
    -s Ii \
    --chr-mapping chr_mapping.txt
```

染色体映射文件格式（空格或Tab分隔）：

```
chr1    Chr1
chr2    Chr2
scaffold_1    Scf1
```

### 跳过AGAT清洗

如果输入GFF文件已经过清洗，可以跳过AGAT步骤以节省时间：

```bash
biopytools gff-renamer \
    -i clean.gff \
    -o output.gff \
    -p CDRT \
    -s Ii \
    --skip-clean
```

### 指定命名格式

```bash
# standard格式: CDRT_Ii06g000010
biopytools gff-renamer -i input.gff -o output.gff -p CDRT -s Ii --naming-format standard

# simple格式: CDRT06G000010
biopytools gff-renamer -i input.gff -o output.gff -p CDRT -s Ii --naming-format simple

# compact格式: CDRT06g000010
biopytools gff-renamer -i input.gff -o output.gff -p CDRT -s Ii --naming-format compact
```

## 命令行参数

### 必需参数

| 参数 | 简写 | 说明 | 示例 |
|------|------|------|------|
| `--input` | `-i` | 输入GFF文件路径 | `complete.genomic.gff` |
| `--output` | `-o` | 输出GFF文件路径 | `output.gff` |
| `--prefix` | `-p` | ID前缀 | `CDRT`, `AGIS` |
| `--species` | `-s` | 物种缩写（2-4个字母） | `Ov`, `Ii`, `Os` |

### 可选参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--threads` | `12` | 并行线程数 |
| `--naming-format` | `standard` | 命名格式：`standard`/`simple`/`compact` |
| `--output-mrna-mapping` | `False` | 输出mRNA新旧ID映射文件 |
| `--mrna-mapping-file` | 自动生成 | 映射文件路径（默认：`{输出文件名}_mrna_mapping.tsv`） |
| `--chr-mapping` | - | 染色体映射文件路径 |
| `--include-utr` | `False` | 在重命名中包含UTR特征 |
| `--skip-clean` | `False` | 跳过AGAT清洗步骤 |

## 输出说明

### 1. 重命名后的GFF文件

所有feature的ID和Parent属性均被替换为规范格式：

```
# 基因
ID=CDRT_Ii06g000010;Name=CDRT_Ii06g000010

# mRNA
ID=CDRT_Ii06g000010.mRNA1;Name=CDRT_Ii06g000010.mRNA1;Parent=CDRT_Ii06g000010

# lncRNA
ID=CDRT_Ii06g000010.lncRNA1;Name=CDRT_Ii06g000010.lncRNA1;Parent=CDRT_Ii06g000010

# exon
ID=CDRT_Ii06g000010.mRNA1.exon1;Parent=CDRT_Ii06g000010.mRNA1

# CDS
ID=CDRT_Ii06g000010.mRNA1.cds1;Parent=CDRT_Ii06g000010.mRNA1

# UTR（需要--include-utr）
ID=CDRT_Ii06g000010.mRNA1.five_prime_UTR1;Parent=CDRT_Ii06g000010.mRNA1
```

### 2. mRNA映射文件

使用 `--output-mrna-mapping` 生成TSV文件，格式为：

```
新ID	旧ID	原始属性1	原始属性2	...
```

可用于追溯新旧ID对应关系。

## 支持的输入格式

| 注释流程 | gene ID示例 | mRNA ID示例 | mRNA编号提取规则 |
|----------|-------------|-------------|------------------|
| NCBI Gnomon | `gene-Iin_Chr06_001435` | `rna-Iin_Chr06_001435-R1` | `-R` 后缀 |
| EGAPx | `gene-XXX` | `mRNA-XXX-R1` | `-R` 后缀 |
| BRAKER | - | `gXXX.t1` | `.t` 后缀 |
| EviAnn | `LOC_00000001` | `LOC_00000001-mRNA-1` | `-mRNA-` 后缀 |

## 工作流程

1. **AGAT清洗**（默认）- 使用AGAT规范化GFF结构
2. **特征收集** - 单次遍历收集所有gene、transcript、exon、CDS等特征
3. **ID映射构建** - 按染色体和基因组位置排序后生成新ID
4. **属性替换** - 并行将所有feature的ID和Parent替换为新ID
5. **输出写入** - 清理冗余属性，按基因组位置排序子特征后写入

## 命名格式对照

| 格式 | gene ID示例 | mRNA ID示例 | 说明 |
|------|-------------|-------------|------|
| `standard` | `CDRT_Ii06g000010` | `CDRT_Ii06g000010.mRNA1` | 前缀_物种染色体编号，基因编号从000010起递增10 |
| `simple` | `CDRT06G000010` | `CDRT06G000010.mRNA1` | 前缀+染色体编号，无物种缩写 |
| `compact` | `CDRT06g000010` | `CDRT06g000010.mRNA1` | 与simple类似，基因编号用小写g |

## 特殊处理

- **假基因（pseudogene）**：Name属性自动添加 `.pseudogene` 后缀
- **纯lncRNA基因**：当基因下仅有lncRNA转录本（无mRNA）时，Name属性添加 `.lncRNA` 后缀
- **AGAT新增特征**：AGAT清洗可能新增UTR等feature，这些feature的Parent会被自动替换为对应的新ID
- **exon/CDS排序**：子特征按基因组位置统一排序后编号

## 故障排除

### 1. 未找到AGAT

```
命令执行失败: AGAT清洗GFF文件
```

**解决方案：**

```bash
conda install -c bioconda agat
```

或跳过清洗步骤：

```bash
biopytools gff-renamer -i input.gff -o output.gff -p CDRT -s Ii --skip-clean
```

### 2. 输入文件格式不支持

```
输入文件格式不正确: input.txt. 支持: .gff, .gff3
```

**解决方案：** 确认输入文件为GFF/GFF3格式。

### 3. 基因编号不连续

基因编号按每条染色体独立计数，从000010开始每次递增10。如果某个基因缺失，编号不会回填，这是正常行为。

## 版本历史

- **v1.0.0** (2025-12-25)
  - 初始版本发布
  - 支持NCBI Gnomon、EGAPx、BRAKER、EviAnn格式
  - 支持三种命名格式
  - 支持AGAT清洗集成
  - 支持并行处理

## 相关工具

- `biopytools chr-renamer` - 染色体名称重命名工具
- `biopytools vcf-renamer` - VCF样品名称重命名工具
- `biopytools egapx-batch` - EGAPx批量注释工具
