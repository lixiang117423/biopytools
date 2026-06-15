# RxLR效应蛋白扫描工具文档
# RxLR Effector Protein Scanner Documentation

## 功能概述 | Feature Overview

RxLR效应蛋白扫描工具用于在蛋白质序列中批量搜索RxLR和EER基序，鉴定卵菌门（如疫霉菌）的候选效应蛋白。

**主要功能|Main Features:**
- 自动读取FASTA格式的蛋白质序列文件
- 在指定窗口（默认21-120位）内搜索RxLR、QxLR、GxLR和EER基序
- 支持序列长度验证和标记
- 同时输出Excel和TSV格式结果
- 自动生成候选蛋白列表

## 快速开始 | Quick Start

### 基本用法 | Basic Usage

```bash
# 扫描蛋白质序列文件
biopytools rxlr-scanner -i proteins.fa -o rxlr_results

# 使用自定义窗口
biopytools rxlr-scanner -i proteins.fa -o rxlr_results --window-start 25 --window-end 150

# 指定输出目录
biopytools rxlr-scanner -i proteins.fa -o rxlr_results --output-dir ./my_output
```

## 参数说明 | Parameters

### 必需参数 | Required Parameters

| 参数 | 描述 | 示例 |
|------|------|------|
| `-i, --input` | 输入FASTA文件路径| `-i proteins.fa` |
| `-o, --output-prefix` | 输出文件前缀| `-o rxlr_results` |

### 搜索窗口参数 | Search Window Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--window-start` | `20` | 窗口起始位置（20对应第21位，0-based索引）|
| `--window-end` | `120` | 窗口结束位置（不包含）|
| `--min-length` | `120` | 最小序列长度，短于此长度的序列会被标记 |

**注意**: 搜索窗口为21-120位，对应Python索引[20:120]

### 输出选项 | Output Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--output-dir` | `./rxlr_scanner_output` | 输出目录路径 |
| `--no-excel` | `False` | 不生成Excel输出 |
| `--no-tsv` | `False` | 不生成TSV输出 |

### 日志选项 | Logging Options

| 参数 | 描述 |
|------|------|
| `-v, --verbose` | 详细输出模式（显示DEBUG级别日志）|
| `--log-file` | 生成日志文件 |

## 输出文件 | Output Files

### 主要输出文件 | Main Output Files

1. **{prefix}.xlsx** - Excel格式结果文件
2. **{prefix}.tsv** - TSV格式结果文件
3. **{prefix}_candidates_only.tsv** - 仅包含候选RxLR蛋白
4. **{prefix}.log** - 日志文件（如果使用--log-file）

### 结果字段说明 | Result Fields Description

| 字段名 | 描述 |
|--------|------|
| `Sequence_ID` | 序列ID |
| `Sequence_Length` | 序列长度 |
| `Valid_Length_(≥120)` | 序列长度是否≥120（Yes/No）|
| `Has_RxLR_motif` | 是否检测到RxLR相关基序（Yes/No）|
| `Has_EER_motif` | 是否检测到EER基序（Yes/No）|
| `RxLR_Candidate` | 是否为候选RxLR蛋白（Yes/No）|
| `Motif_Types` | 检测到的基序类型，用分号分隔 |
| `Motif_Positions` | 基序在原序列中的位置（1-based），用分号分隔 |
| `Total_Motif_Count` | 检测到的基序总数 |

### 基序类型 | Motif Types

- **RxLR**: 精氨酸-任意氨基酸-亮氨酸-精氨酸
- **QxLR**: 谷氨酰胺-任意氨基酸-亮氨酸-精氨酸
- **GxLR**: 甘氨酸-任意氨基酸-亮氨酸-精氨酸
- **EER**: 谷氨酸-谷氨酸-精氨酸（连续）

### 判定标准 | Criteria

**候选RxLR蛋白**: 检测到RxLR、QxLR、GxLR **或** EER任一基序

## 使用示例 | Usage Examples

### 示例1：基本扫描 | Example 1: Basic Scanning

```bash
# 扫描蛋白质文件
biopytools rxlr-scanner -i proteins.fa -o sample_results

# 输出文件：
# - sample_results.xlsx
# - sample_results.tsv
# - sample_results_candidates_only.tsv
```

### 示例2：自定义搜索窗口 | Example 2: Custom Search Window

```bash
# 扩大搜索窗口到25-150位
biopytools rxlr-scanner \
    -i proteins.fa \
    -o rxlr_results \
    --window-start 24 \
    --window-end 150
```

### 示例3：严格筛选 | Example 3: Strict Filtering

```bash
# 提高最小长度要求到150aa
biopytools rxlr-scanner \
    -i proteins.fa \
    -o rxlr_results \
    --min-length 150
```

### 示例4：仅TSV输出 | Example 4: TSV Only Output

```bash
# 不生成Excel文件，只输出TSV
biopytools rxlr-scanner \
    -i proteins.fa \
    -o rxlr_results \
    --no-excel
```

## 输入文件格式 | Input File Format

标准FASTA格式蛋白质序列文件：

```fasta
>protein1
MKTIIALSYIFCLVFA...
>protein2
MTYKLILALGAAATA...
```

**要求**:
- 标准FASTA格式
- 序列为氨基酸序列（标准20种氨基酸）
- 序列ID在`>`后的第一个词（空格分隔）

## 技术细节 | Technical Details

### 搜索算法 | Search Algorithm

1. **序列读取**: 逐条解析FASTA文件
2. **长度验证**: 检查序列长度是否≥120aa
3. **窗口提取**: 提取21-120位氨基酸片段
4. **基序搜索**:
   - 使用正则表达式搜索RxLR、QxLR、GxLR模式
   - 搜索EER连续三联体
5. **结果判定**: RxLR或EER任一匹配即为候选

### 正则表达式模式 | Regex Patterns

- `R.LR` - RxLR模式
- `Q.LR` - QxLR模式
- `G.LR` - GxLR模式
- `EER` - EER模式

### 边界处理 | Boundary Handling

- 序列长度<窗口起始：无法提取窗口，标记为无效长度
- 序列长度在窗口内：提取可用部分
- 所有序列都会被记录，但短序列会被标记

## 结果解读 | Result Interpretation

### 筛选候选蛋白 | Filtering Candidate Proteins

```bash
# 查看候选蛋白列表
cat rxlr_scanner_output/rxlr_results_candidates_only.tsv

# 或在Excel中筛选"RxLR_Candidate"列为"Yes"的行
```

### 统计信息 | Statistics

工具会在运行结束时输出：
- 总序列数
- 有效序列数（≥120aa）
- 无效序列数（<120aa）
- 候选RxLR蛋白数
- 候选率

## 常见问题 | FAQ

**Q: 为什么有些序列长度<120但仍然有结果？**
A: 工具会标记所有短序列，但仍然尝试搜索基序。短序列的"Valid_Length"字段会显示"No"。

**Q: 如何提高检测特异性？**
A: 可以要求同时具有RxLR和EER基序，通过筛选结果实现：
```bash
# 使用awk筛选同时具有两种基序的序列
awk -F'\t' '$5=="Yes" && $6=="Yes"' rxlr_results.tsv
```

**Q: 搜索窗口可以调整吗？**
A: 可以。使用--window-start和--window-end参数调整。注意参数是0-based索引。

**Q: 是否支持其他基序模式？**
A: 当前版本仅支持RxLR、QxLR、GxLR和EER。如需添加其他模式，请联系开发者。

## 系统要求 | System Requirements

### 依赖软件 | Dependencies

- **Python**: 3.7+
- **Python包**:
  - `pandas` - 数据处理
  - `openpyxl` - Excel输出
  - `click` - 命令行界面

### 安装依赖 | Installing Dependencies

```bash
# 使用pip安装
pip install pandas openpyxl click
```

## 注意事项 | Important Notes

1. **序列质量**: 确保输入序列为完整蛋白质序列，不含片段
2. **窗口选择**: 默认窗口适用于大多数情况，但可根据信号肽长度调整
3. **结果验证**: 候选蛋白建议进一步验证（如结构预测、功能注释）
4. **性能**: 大型蛋白组（>10000条序列）可能需要较长处理时间

## 引用 | Citation

如果在研究中使用此工具，请引用相关文献：

```
Win, J., Kamoun, S., & Jones, A. M. (2007).
RxLR effectors: weapons of plant pathogen infection.
Trends in Microbiology, 15(8), 341-348.
```

## 版本历史 | Version History

| 版本 | 日期 | 主要变更 |
|------|------|----------|
| 1.0.0 | 2026-02-05 | 初始版本 - 基本RxLR/EER扫描功能 |

---

**开发者**: Xiang LI
**最后更新**: 2026-02-05
