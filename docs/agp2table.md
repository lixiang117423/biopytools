# AGP转表格工具 | AGP to Table Converter

**专业的AGP格式转换工具 | Professional AGP Format Conversion Tool**

## 功能概述 | Overview

AGP转表格工具是一个简单易用的AGP格式文件转换工具，可以将基因组组装的AGP格式文件转换为更易读的表格格式，自动添加表头和统计信息，便于查看和分析基因组scaffolding结果。

## 主要特性 | Key Features

- **📊 格式转换**: 支持TXT、TSV、CSV多种输出格式
- **📋 自动表头**: 自动添加中英文双语表头
- **📈 统计信息**: 自动生成详细的scaffold统计信息
- **🎯 灵活分组**: 支持按scaffold分组输出
- **🔍 格式验证**: 自动检测和跳过格式错误的行
- **📝 详细日志**: 完整的处理过程日志记录

## 快速开始 | Quick Start

### 基本用法 | Basic Usage

```bash
# 转换AGP文件为表格格式
biopytools agp2table -i assembly.agp -o assembly_table.txt

# 输出为CSV格式
biopytools agp2table -i assembly.agp -o assembly.csv -f csv

# 不包含统计信息
biopytools agp2table -i assembly.agp -o output.txt --no-statistics
```

## 参数说明 | Parameters

### 必需参数 | Required Parameters

| 参数 | 描述 | 示例 |
|------|------|------|
| `-i, --input` | AGP文件路径 | `-i assembly.agp` |
| `-o, --output` | 输出表格文件路径 | `-o output.txt` |

### 可选参数 | Optional Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-f, --format` | `txt` | 输出格式 (txt/tsv/csv) |
| `--no-statistics` | `False` | 不添加统计信息 |
| `--no-headers` | `False` | 不添加表头 |
| `--no-grouping` | `False` | 不按scaffold分组 |

## 输出格式说明 | Output Format Description

### 表格列名 | Table Columns

| 列名 | 描述 |
|------|------|
| **Scaffold** | Scaffold/Chromosome名称 |
| **Start** | 起始位置 |
| **End** | 结束位置 |
| **Length** | 长度(bp) |
| **Part_Number** | 部分编号 |
| **Type** | 组件类型 (Contig/Gap) |
| **Component_ID** | 组件ID (contig名称) |
| **Component_Start** | 组件起始位置 |
| **Component_End** | 组件结束位置 |
| **Orientation** | 方向 (+/-/?) |

### 组件类型 | Component Types

| 类型 | 描述 |
|------|------|
| **W** | Contig - 已知序列 |
| **U/N** | Gap - 未知序列(通常是100bp的N) |
| **O** | Other - 其他类型 |

## 使用示例 | Usage Examples

### 示例1：基本转换

```bash
biopytools agp2table \
    -i EcA_corrected.agp \
    -o EcA_table.txt
```

输出示例：
```
Scaffold名称	起始位置	结束位置	长度(bp)	部分编号	类型	组件ID	组件起始	组件结束	方向
group1	1	190127	190127	1	Contig	ptg002519l	1	190127	+
group1	190128	190227	100	2	Gap	100	scaffold	yes	proximity_ligation
group1	190228	5516498	5326271	3	Contig	ptg000008l	1	5326271	-
```

### 示例2：CSV格式输出

```bash
biopytools agp2table \
    -i assembly.agp \
    -o assembly.csv \
    -f csv
```

### 示例3：不包含统计信息和表头

```bash
biopytools agp2table \
    -i assembly.agp \
    -o raw_data.txt \
    --no-statistics \
    --no-headers
```

## 统计信息说明 | Statistics Description

转换完成后，输出文件末尾会包含以下统计信息：

```
====================================================================================================
统计信息|Statistics
====================================================================================================

总记录数|Total records: 7792
总scaffold数|Total scaffolds: 7750
总长度|Total length: 256,789,234 bp

Contig数量|Contig count: 7784
Gap数量|Gap count: 8
Gap总长度|Gap total length: 800 bp

组件类型统计|Component type statistics:
  W: 7784
  U: 8

前20个最长的scaffold|Top 20 longest scaffolds:
   1. group1: 27,718,774 bp
   2. group2: 26,402,014 bp
   3. group3: 23,266,732 bp
   ...
```

## 输入文件格式 | Input File Format

标准AGP格式文件示例：

```
group1	1	190127	1	W	ptg002519l	1	190127	+
group1	190128	190227	2	U	100	scaffold	yes	proximity_ligation
group1	190228	5516498	3	W	ptg000008l	1	5326271	-
```

AGP格式说明：
- 第1列: object名称(scaffold/chromosome)
- 第2-3列: 起始和结束位置
- 第4列: 部分编号
- 第5列: 组件类型 (W=contig, U=gap)
- 第6-9列: 组件详细信息

## 注意事项 | Important Notes

1. **文件编码**: 输出文件使用UTF-8编码
2. **日志文件**: 会在输出目录生成`agp2table.log`日志文件
3. **格式验证**: 自动跳过格式不正确的行
4. **内存使用**: 大文件处理时需要注意内存使用

## 故障排除 | Troubleshooting

### 常见问题

**Q: 提示"文件不存在"**
```bash
# 检查文件路径是否正确
ls -lh /path/to/assembly.agp

# 使用绝对路径
biopytools agp2table -i /full/path/to/assembly.agp -o output.txt
```

**Q: 输出文件为空**
```bash
# 检查AGP文件格式是否正确
head -20 assembly.agp

# 查看日志文件
cat agp2table.log
```

## 系统要求 | System Requirements

- Python 3.7+
- click库

## 许可证 | License

MIT License
