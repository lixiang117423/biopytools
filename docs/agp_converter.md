# AGP文件转换器 | AGP File Converter

## 概述 | Overview

AGP文件转换器是一个用于将机器可读的AGP (A Golden Path) 文件转换为人类易读表格格式的工具。支持转换为TSV、CSV、Excel等格式，并提供丰富的数据过滤、排序和分析功能。

## AGP文件格式 | AGP File Format

AGP文件是一种描述基因组组装结构的文件格式，包含以下主要字段：

| 字段 | 说明 | 示例 |
|------|------|------|
| 染色体ID | 目标染色体或scaffold的标识符 | chr1 |
| 起始位置 | 组装片段的起始坐标 | 1 |
| 结束位置 | 组装片段的结束坐标 | 10000 |
| 组装顺序 | 组装片段的顺序号 | 1 |
| 组件类型 | W(Contig), U(Gap), N(Gap) | W |
| 组件ID | Contig或Gap的标识符 | contig_1 |
| 组件起始 | 组件的起始位置 | 1 |
| 组件结束 | 组件的结束位置 | 5000 |
| 方向 | 组装方向 (+或-) | + |

## 功能特点 | Features

- 📊 **多格式输出**: 支持TSV、CSV、Excel等格式
- 🔍 **智能过滤**: 按长度、染色体等条件过滤数据
- 📈 **灵活排序**: 按长度、染色体等字段排序
- 🌐 **双语支持**: 支持中英文表头和提示
- 📝 **详细日志**: 完整的处理过程记录
- ⚡ **高性能**: 优化的大文件处理能力
- 🛡️ **错误处理**: 完善的异常处理和恢复机制

## 安装要求 | Requirements

- Python 3.6+
- 无额外依赖 (纯Python实现)

## 使用方法 | Usage

### 基本用法 | Basic Usage

```bash
# 输出到屏幕
biopytools agp-converter -i genome.agp

# 输出到TSV文件
biopytools agp-converter -i genome.agp -o result.tsv

# 输出到Excel文件
biopytools agp-converter -i genome.agp -o result.xlsx
```

### 参数说明 | Parameters

#### 必需参数 | Required Parameters

| 参数 | 说明 | 示例 |
|------|------|------|
| `-i, --input` | 输入AGP文件路径 | `genome.agp` |

#### 可选参数 | Optional Parameters

| 参数 | 默认值 | 说明 | 示例 |
|------|--------|------|------|
| `-o, --output` | 屏幕输出 | 输出文件路径 | `result.tsv` |
| `--output-format` | `tsv` | 输出格式 (tsv/csv/xlsx/xls) | `xlsx` |
| `--encoding` | `utf-8` | 文件编码格式 | `gbk` |
| `--precision` | `4` | Mb长度小数位数 | `6` |
| `--include-gaps` | `False` | 是否包含Gap记录 | `--include-gaps` |
| `--sort-by` | 无 | 排序方式 | `length_desc` |
| `--filter-by` | 无 | 过滤条件 | `"length>1000"` |

#### 日志控制 | Logging Control

| 参数 | 说明 | 示例 |
|------|------|------|
| `-v, --verbose` | 详细输出模式 (-v: INFO, -vv: DEBUG) | `-vv` |
| `--quiet` | 静默模式 (只输出ERROR) | `--quiet` |
| `--log-file` | 日志文件名 | `convert.log` |
| `--log-dir` | 日志目录 | `./logs` |

## 输出格式 | Output Formats

### TSV格式 (默认)
```
染色体ID    组装顺序    原始Contig名称    方向    操作说明    长度(bp)    长度(Mb)    染色体坐标范围
chr1    1    contig_001    正向 (+)    直接拼接    5000    0.0050    1-5000
chr1    2    contig_002    反向 (-)    翻转后拼接    3000    0.0030    5001-8000
```

### CSV格式
```csv
染色体ID,组装顺序,原始Contig名称,方向,操作说明,长度(bp),长度(Mb),染色体坐标范围
chr1,1,contig_001,正向 (+),直接拼接,5000,0.0050,1-5000
chr1,2,contig_002,反向 (-),翻转后拼接,3000,0.0030,5001-8000
```

### Excel格式
- 可以直接用Microsoft Excel打开
- 自动识别列类型
- 支持中文显示

## 过滤功能 | Filtering

### 按长度过滤
```bash
# 过滤长度大于1000bp的记录
biopytools agp-converter -i genome.agp -o filtered.tsv --filter-by "length>1000"

# 过滤长度小于5000bp的记录
biopytools agp-converter -i genome.agp -o filtered.tsv --filter-by "length<5000"
```

### 按染色体过滤
```bash
# 只保留chr1的记录
biopytools agp-converter -i genome.agp -o chr1.tsv --filter-by "chr=chr1"

# 只保留chr2的记录
biopytools agp-converter -i genome.agp -o chr2.tsv --filter-by "chr=chr2"
```

## 排序功能 | Sorting

```bash
# 按长度降序排序
biopytools agp-converter -i genome.agp -o sorted.tsv --sort-by length_desc

# 按长度升序排序
biopytools agp-converter -i genome.agp -o sorted.tsv --sort-by length_asc

# 按染色体ID排序
biopytools agp-converter -i genome.agp -o sorted.tsv --sort-by chromosome
```

## 使用示例 | Examples

### 1. 基本转换
```bash
# 将AGP文件转换为易读的表格并输出到屏幕
biopytools agp-converter -i genome.agp
```

### 2. Excel格式输出
```bash
# 输出为Excel格式，方便进一步分析
biopytools agp-converter -i genome.agp -o genome_analysis.xlsx --output-format xlsx
```

### 3. 包含Gap记录
```bash
# 包含Gap记录，查看完整的组装结构
biopytools agp-converter -i genome.agp -o complete.tsv --include-gaps
```

### 4. 长片段分析
```bash
# 只保留长度大于1Mb的片段，并按长度降序排序
biopytools agp-converter -i genome.agp -o large_fragments.tsv \
    --filter-by "length>1000000" --sort-by length_desc
```

### 5. 特定染色体分析
```bash
# 分析chr1染色体，输出为CSV格式
biopytools agp-converter -i genome.agp -o chr1_analysis.csv \
    --filter-by "chr=chr1" --output-format csv
```

### 6. 详细日志
```bash
# 启用详细日志输出，了解处理过程
biopytools agp-converter -i genome.agp -o result.tsv -vv --log-file detailed.log
```

### 7. 高精度输出
```bash
# 提高Mb长度的精度到小数点后6位
biopytools agp-converter -i genome.agp -o precise.tsv --precision 6
```

## 输出字段说明 | Output Fields

| 字段名 | 说明 | 示例 |
|--------|------|------|
| 染色体ID | 目标染色体或scaffold标识 | chr1 |
| 组装顺序 | 组装片段的顺序号 | 1 |
| 原始Contig名称 | Contig或组件的标识符 | contig_001 |
| 方向 | 组装方向 | 正向 (+) / 反向 (-) |
| 操作说明 | 组装操作的描述 | 直接拼接 / 翻转后拼接 |
| 长度(bp) | 片段长度（碱基对） | 5000 |
| 长度(Mb) | 片段长度（兆碱基） | 0.0050 |
| 染色体坐标范围 | 在染色体上的位置范围 | 1-5000 |

## 性能优化 | Performance Optimization

### 大文件处理
```bash
# 处理大文件时，建议输出到文件而不是屏幕
biopytools agp-converter -i large_genome.agp -o large_result.tsv --quiet
```

### 内存优化
- 工具采用流式处理，支持处理GB级别的AGP文件
- 内存使用量与文件大小无关，适合大基因组分析

### 编码处理
- 支持多种文件编码格式
- 自动处理中文字符，避免乱码

## 错误处理 | Error Handling

### 常见错误及解决方案

#### 1. 文件编码错误
```
错误: 文件编码错误 | File encoding error
解决: 指定正确的编码格式 --encoding gbk
```

#### 2. 文件格式错误
```
错误: 文件不包含有效的AGP格式行 | File does not contain valid AGP format lines
解决: 检查文件格式，确保是标准的AGP文件
```

#### 3. 权限错误
```
错误: 文件写入错误 | File write error
解决: 检查输出目录的写入权限
```

## 日志分析 | Log Analysis

### 查看处理统计
```bash
# 查看日志文件获取详细统计信息
cat logs/agp_converter.log | grep "INFO"
```

### 错误诊断
```bash
# 查看错误信息
cat logs/agp_converter.log | grep "ERROR"
```

## 技术细节 | Technical Details

### 解析算法
- 基于正则表达式的字段解析
- 智能的类型推断和数据清洗
- 完善的错误恢复机制

### 数据验证
- AGP格式规范检查
- 字段完整性验证
- 数据一致性检查

### 格式化引擎
- 可插拔的格式化器架构
- 支持自定义输出格式
- 内存高效的表格处理

## 更新日志 | Changelog

### v1.0.0 (2024-12-22)
- ✨ 初始版本发布
- 📊 支持TSV、CSV、Excel格式输出
- 🔍 智能过滤和排序功能
- 🌐 完整的中英文支持
- 🛡️ 完善的错误处理和日志记录

## 贡献指南 | Contributing

欢迎提交Issue和Pull Request来改进这个工具。

## 许可证 | License

本项目采用MIT许可证。

## 联系方式 | Contact

如有问题或建议，请联系生信分析团队。