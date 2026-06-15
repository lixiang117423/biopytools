# Primer3 引物设计分析模块

**专业的PCR引物批量设计工具 | Professional Batch PCR Primer Design Tool**

## 功能概述 | Overview

Primer3 引物设计分析模块是一个强大的PCR引物批量设计工具，基于Primer3软件构建，支持从FASTA序列文件批量设计引物。提供灵活的参数配置（引物长度、退火温度、产物大小等），支持批量处理和完整的引物质量评估，适用于各种分子生物学研究。

## 主要特性 | Key Features

- **批量引物设计**: 支持从FASTA文件批量设计引物，一次处理多条序列
- **灵活参数配置**: 支持自定义引物长度范围（默认20-22 bp）、退火温度范围（默认58±5°C，即53-63°C）
- **完整质量评估**: 输出包含Tm值、GC含量、二聚体等所有Primer3质量指标
- **多格式输出**: 支持CSV、TSV、Excel格式输出，方便后续分析
- **Conda环境集成**: 自动检测并使用conda环境中的Primer3
- **详细日志记录**: 完整的处理过程日志和错误追踪

## 快速开始 | Quick Start

### 基本用法 | Basic Usage

```bash
# 基本引物设计（默认：method=all，引物在序列两端，自动调整产物大小）
biopytools primer3 -i sequences.fasta -o primer3_output

# 随机设计引物（method=random，在序列任意位置设计，覆盖率更高）
biopytools primer3 -i sequences.fasta -o primer3_output -m random

# 自定义引物长度和退火温度
biopytools primer3 -i sequences.fasta -o primer3_output \
    --primer-min-size 18 \
    --primer-max-size 24 \
    --primer-opt-tm 60.0

# 禁用自动产物大小调整，使用固定范围
biopytools primer3 -i sequences.fasta -o output \
    --no-auto-product-size \
    --product-min-size 500 \
    --product-max-size 2000

# 自定义产物大小比例
biopytools primer3 -i sequences.fasta -o output \
    --product-size-min-ratio 0.3 \
    --product-size-max-ratio 0.9
```

### 高级用法 | Advanced Usage

```bash
# 自定义所有参数的完整引物设计
biopytools primer3 -i sequences.fasta -o primer3_output \
    --primer-min-size 20 \
    --primer-opt-size 22 \
    --primer-max-size 24 \
    --primer-min-tm 55.0 \
    --primer-opt-tm 58.0 \
    --primer-max-tm 62.0 \
    --product-min-size 150 \
    --product-max-size 400 \
    --primer-num-return 10 \
    --output-format xlsx
```

## 参数说明 | Parameters

### 必需参数 | Required Parameters

| 参数 | 描述 | 示例 |
|------|------|------|
| `-i, --input-fasta` | 输入FASTA文件路径 | `-i sequences.fasta` |
| `-o, --output-dir` | 输出目录路径 | `-o primer3_output` |

### 引物长度参数 | Primer Size Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--primer-min-size` | `20` | 最小引物长度(bp) |
| `--primer-opt-size` | `20` | 最优引物长度(bp) |
| `--primer-max-size` | `22` | 最大引物长度(bp) |

### 退火温度参数 | Annealing Temperature Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--primer-min-tm` | `53.0` | 最小退火温度(°C) |
| `--primer-opt-tm` | `58.0` | 最优退火温度(°C) |
| `--primer-max-tm` | `63.0` | 最大退火温度(°C) |

**注**: 默认退火温度范围为58±5°C，即53-63°C

### 产物大小参数 | Product Size Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--product-min-size` | `100` | 最小产物大小(bp) |
| `--product-max-size` | `300` | 最大产物大小(bp) |

### 其他参数 | Other Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--method`, `-m` | `all` | 引物设计策略(all:覆盖头尾, random:随机设计) |
| `--primer-num-return` | `5` | 每条序列返回的引物对数量 |
| `--output-format` | `csv` | 输出格式(csv/tsv/xlsx) |
| `--output-header-lang` | `zh` | 输出表头语言(zh:中文, en:英文) |
| `--primer-end-margin` | `200` | 两端允许的引物位置范围bp(仅用于method=all) |
| `--auto-product-size` | `True` | 自动根据序列长度设置产物大小范围(**默认开启**) |
| `--product-size-min-ratio` | `0.5` | 产物最小长度占序列长度的比例 |
| `--product-size-max-ratio` | `1.0` | 产物最大长度占序列长度的比例 |

### 引物设计策略参数说明 | Primer Design Strategy Parameters

**`--method` / `-m`**: 控制引物设计策略

**`--method all`** (默认): 覆盖头尾模式
- **行为**: 正向引物设计在序列起始端，反向引物设计在序列终止端，确保扩增产物覆盖整个序列
- **优点**: 适用于需要扩增完整序列或验证全长的应用（如基因克隆、全长扩增）
- **缺点**: 序列覆盖率较低（约30%），引物对数较少
- **配合参数**: `--primer-end-margin`控制引物在两端的搜索范围

**`--method random`**: 随机设计模式
- **行为**: Primer3在序列任意位置设计引物，不受头尾限制
- **优点**: 序列覆盖率极高（约99%），引物对数多，引物质量更好
- **缺点**: 引物可能不在序列两端，扩增产物可能不覆盖完整序列
- **适用场景**: PCR检测、基因分型、定点扩增等不需要全长的应用

**性能对比示例**（基于941条测试序列）:
| 模式 | 引物对数 | 序列覆盖率 | 成功率 | 适用场景 |
|------|----------|-----------|--------|---------|
| method=all | 1114 | 280条 | 29.8% | 需要扩增完整序列 |
| method=random | 3764 | 938条 | 99.7% | 序列任意位置设计引物 |

**`--primer-end-margin`**: 控制method=all模式下引物在两端的搜索范围
- 值越小，引物越靠近序列两端
- 建议值：100-500 bp，根据序列长度调整
- 对于短序列(<500bp)，建议设置为序列长度的20-30%
- **注**: 此参数仅在method=all时生效

**`--auto-product-size`** (默认开启): 自动根据序列长度设置产物大小范围
- **默认行为**: 产物大小根据序列长度自动调整
- 产物最小长度 = 序列长度 × `--product-size-min-ratio` (默认0.5)
- 产物最大长度 = 序列长度 × `--product-size-max-ratio` (默认1.0)
- 不会低于全局最小值或超过全局最大值

**`--product-size-min-ratio`**: 产物最小长度占序列长度的比例
- 默认0.5，即产物最小为序列长度的50%
- 范围：0.1-1.0

**`--product-size-max-ratio`**: 产物最大长度占序列长度的比例
- 默认1.0，即产物最大为序列全长
- 范围：0.1-1.0

## 输入文件格式 | Input File Format

### FASTA序列文件 | FASTA Sequence File

标准FASTA格式的序列文件：

```fasta
>sequence1
ATCGATCGATCGATCGATCGATCGATCG...
>sequence2
GCTAGCTAGCTAGCTAGCTAGCTAGCTA...
>gene_X_target
ATGGCCATGGCGCCATGCGATACGCTAGCT...
```

**文件要求**:
- 标准FASTA格式
- 序列使用大写或小写字母（会自动转换为大写）
- 序列ID不能包含特殊字符（会自动转换为下划线）

## 输出结果 | Output Results

### 表头语言选项 | Header Language Options

**中文表头（默认）| Chinese Headers (Default)**:
```
序列ID, 输入序列, 引物对编号, 正向引物, 反向引物, 正向引物长度, 反向引物长度,
正向引物Tm值, 反向引物Tm值, 退火温度, 正向GC含量(%), 反向GC含量(%),
产物大小, 产物Tm值, 惩罚值, ...
```

**英文表头|English Headers**:
```
Sequence_ID, Input_Sequence, Pair_Number, F_Primer, R_Primer,
F_Primer_Length, R_Primer_Length, F_Primer_TM, R_Primer_TM,
Annealing_Temp, F_GC_Percent, R_GC_Percent, Product_Size_bp,
Product_TM, Penalty, ...
```

使用`--output-header-lang`参数切换表头语言。

输出目录包含以下文件：

```
primer3_output/
├── primers_result.csv          # 引物设计结果表格（默认中文表头）
└── primer3_design.log          # 运行日志
```

**注**: 默认使用中文表头，可通过`--output-header-lang en`参数切换为英文表头。

### 结果表格列说明 | Result Table Columns

| 列名 | 描述 |
|------|------|
| `ID` | 序列ID |
| `Input_Seq` | 输入序列 |
| `Pair_Num` | 引物对编号 |
| `F_Primer` | 正向引物序列 |
| `R_Primer` | 反向引物序列 |
| `F_Primer_Length` | 正向引物长度 |
| `R_Primer_Length` | 反向引物长度 |
| `F_Primer_TM` | 正向引物Tm值 |
| `R_Primer_TM` | 反向引物Tm值 |
| `Annealing_Temp` | 推荐退火温度 |
| `F_GC_Percent` | 正向引物GC含量(%) |
| `R_GC_Percent` | 反向引物GC含量(%) |
| `Product_Size` | 产物大小(bp) |
| `Product_TM` | 产物Tm值 |
| `Penalty` | 引物惩罚值(越小越好) |
| `F_Self_Any_TH` | 正向引物自身二聚体(any) |
| `R_Self_Any_TH` | 反向引物自身二聚体(any) |
| `Cross_Dimer` | 引物间二聚体(any) |
| `F_Self_End_TH` | 正向引物自身二聚体(end) |
| `R_Self_End_TH` | 反向引物自身二聚体(end) |
| `Cross_End_Dimer` | 引物间二聚体(end) |

**注**: 所有质量指标都会输出，方便用户根据需要筛选

## 使用示例 | Usage Examples

### 示例1：基本引物设计 | Example 1: Basic Primer Design

```bash
# 对基因序列进行引物设计
biopytools primer3 -i genes.fasta -o gene_primers
```

### 示例2：自定义引物长度 | Example 2: Custom Primer Size

```bash
# 设计较长的引物（22-26 bp）
biopytools primer3 -i sequences.fasta -o output \
    --primer-min-size 22 \
    --primer-opt-size 24 \
    --primer-max-size 26
```

### 示例3：调整退火温度范围 | Example 3: Adjust Temperature Range

```bash
# 设计高Tm值的引物（60-65°C）
biopytools primer3 -i sequences.fasta -o output \
    --primer-min-tm 60.0 \
    --primer-opt-tm 62.0 \
    --primer-max-tm 65.0
```

### 示例4：大产物引物设计 | Example 4: Large Product Primer Design

```bash
# 设计扩增大片段的引物（500-1000 bp）
biopytools primer3 -i sequences.fasta -o output \
    --product-min-size 500 \
    --product-max-size 1000
```

### 示例5：获得更多候选引物 | Example 5: Get More Candidate Primers

```bash
# 每条序列返回10对引物
biopytools primer3 -i sequences.fasta -o output \
    --primer-num-return 10 \
    --output-format xlsx
```

### 示例6：输出英文表头 | Example 6: Output English Headers

```bash
# 使用英文表头输出结果
biopytools primer3 -i sequences.fasta -o output \
    --output-header-lang en
```

### 示例7：随机设计引物（不限制头尾）| Example 7: Random Primer Design

```bash
# 使用method=random模式，在序列任意位置设计引物（高覆盖率）
biopytools primer3 -i sequences.fasta -o output \
    -m random
```

### 示例8：覆盖头尾设计（扩增完整序列）| Example 8: Cover Entire Sequence

```bash
# 使用method=all模式，强制引物在序列两端（扩增完整序列）
biopytools primer3 -i sequences.fasta -o output \
    --method all \
    --primer-end-margin 150
```

## 系统要求 | System Requirements

### 依赖软件 | Dependencies

- **Primer3** (版本 2.6.1 或更新)
  - conda环境: `primer3_v.2.6.1`
  - 可通过以下命令安装:
    ```bash
    conda create -n primer3_v.2.6.1 -c bioconda primer3
    ```
- **Python** (版本 3.7+)
- **Python包**:
  - `click` - 命令行界面
  - `biopython` - FASTA文件解析
  - `pandas` - 结果表格处理
  - `openpyxl` - Excel输出（可选）

### 安装依赖软件 | Installing Dependencies

```bash
# 创建conda环境并安装Primer3
conda create -n primer3_v.2.6.1 -c bioconda primer3

# 激活环境
conda activate primer3_v.2.6.1

# 验证安装
primer3_core --version
```

### 安装Python包 | Installing Python Packages

```bash
# 安装必需的Python包
pip install click biopython pandas

# 如需Excel输出，安装openpyxl
pip install openpyxl
```

## 注意事项 | Important Notes

1. **序列质量**: 确保输入序列质量良好，避免大量N或低复杂度区域
2. **参数合理性**: 引物长度、温度等参数应在合理范围内
3. **输出筛选**: 结果表格包含所有质量指标，可根据需要筛选
4. **路径规范**: 使用相对路径或~路径，避免硬编码绝对路径

## 故障排除 | Troubleshooting

### 常见问题 | Common Issues

**Q: "primer3_core: command not found" 错误**
```bash
# 检查Primer3是否安装
which primer3_core

# 如未找到，安装Primer3
conda create -n primer3_v.2.6.1 -c bioconda primer3
```

**Q: 引物设计失败，没有输出结果**
```bash
# 检查输入序列格式
head -20 sequences.fasta

# 尝试放宽参数限制
biopytools primer3 -i sequences.fasta -o output \
    --primer-min-size 18 \
    --primer-max-size 26 \
    --product-min-size 80 \
    --product-max-size 500
```

**Q: 某些序列没有设计出引物**
- 这是正常情况，可能序列特征不满足参数要求
- 可以尝试调整参数或检查序列质量
- 查看日志文件了解详细信息

## 结果解读指南 | Result Interpretation Guide

### 引物质量评估标准 | Primer Quality Assessment Standards

**高质量引物特征**:
1. **Penalty值**: 越低越好，通常<1.0为优质
2. **Tm值**: 正反向引物Tm值接近（差值<2°C）
3. **GC含量**: 40-60%为最佳
4. **二聚体**: Self_Any_TH和Cross_Dimer值越低越好
5. **产物大小**: 符合预期范围

### 引物筛选建议 | Primer Filtering Suggestions

```python
# 使用pandas筛选高质量引物
import pandas as pd

df = pd.read_csv('primers_result.csv')

# 筛选条件
high_quality = df[
    (df['Penalty'] < 1.0) &
    (df['F_GC_Percent'] >= 40) &
    (df['F_GC_Percent'] <= 60) &
    (df['R_GC_Percent'] >= 40) &
    (df['R_GC_Percent'] <= 60) &
    (abs(df['F_Primer_TM'] - df['R_Primer_TM']) < 2.0)
]

# 保存筛选结果
high_quality.to_csv('high_quality_primers.csv', index=False)
```

## 相关资源 | Related Resources

- [Primer3官方文档](http://primer3.org/manual.html)
- [Primer3 GitHub仓库](https://github.com/primer3-org/primer3)
- [Primer3参数详解](http://primer3.org/manual.html#globalTags)
- [引物设计最佳实践](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3484727/)

## 引用信息 | Citation

如果在学术研究中使用此工具，请引用Primer3相关文献：

```
Untergasser A, Cutcutache I, Koressaar T, Ye J, Faircloth BC, Remm M and Rozen SG.
Primer3--new capabilities and interfaces.
Nucleic Acids Res. 2012 Aug 1;40(15):e115.
```

以及Primer3核心算法文献：

```
Koressaar T and Remm M.
Enhancements and modifications of primer design program Primer3.
Bioinformatics 2007;23(10):1289-1291.
```
