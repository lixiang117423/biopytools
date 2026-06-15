# 染色体重命名工具

基于minimap2全基因组比对的染色体重命名工具

## 功能概述

染色体重命名工具通过minimap2全基因组比对，自动识别group命名方式的染色体与参考基因组标准Chr命名的对应关系，实现染色体名称的批量转换。适用于基因组挂载后的命名标准化，为后续比较基因组分析提供便利。

## 主要特性

- **迭代映射策略**: 从高覆盖度阈值到低阈值自动迭代，优先使用高质量映射
  - 第1轮 (>=90%): 最高质量序列优先映射
  - 第2轮 (>=80%): 高质量序列
  - 第3轮 (>=70%): 较高质量序列
  - 第4轮 (>=60%): 中等质量序列
  - 第5轮 (>=50%): 中低质量序列
  - 第6轮 (>=40%): 低质量序列
  - 第7轮 (>=30%): 更低质量序列
  - 第8轮 (>=20%): 最低质量序列
  - **已映射的序列不参与后续轮次**，确保高质量映射优先
- 基于minimap2的全基因组比对（支持asm5/asm10/asm20模式）
- 智能映射关系建立（处理一对一、一对多、多对一等复杂情况）
- 详细的映射表和汇总报告输出（包含每轮映射统计）
- 完整的日志记录和错误追踪
- 自动检测并提示嵌合/错误组装和碎片化组装

## 快速开始

### 基本用法

```bash
# 使用默认参数进行染色体重命名
biopytools chr_rename \
    -r reference.fa \
    -q query.fa \
    -o output_dir
```

### 高级用法

```bash
# 自定义比对参数和过滤阈值
biopytools chr_rename \
    -r reference.fa \
    -q query.fa \
    -x asm10 \
    -c 0.8 \
    -i 0.95 \
    -t 24 \
    -o output_dir
```

## 参数说明

### 必需参数

| 参数 | 描述 | 示例 |
|------|------|------|
| `-r, --ref` | 参考基因组FASTA文件路径 | `-r reference.fa` |
| `-q, --query` | 待重命名的基因组FASTA文件路径 | `-q query.fa` |

### 可选参数

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-o, --output-dir` | `./chr_rename_output` | 输出目录路径 |
| `-a, --minimap2-path` | `minimap2` | minimap2软件路径 |
| `-x, --preset` | `asm5` | minimap2预设模式 (asm5/asm10/asm20) |
| `-t, --threads` | `12` | 线程数 |
| `-i, --min-identity` | `0.9` | 最小序列一致性阈值 (0-1) |

**注意**: 工具使用自动迭代映射策略，会从高覆盖度（90%）到低覆盖度（20%）自动尝试映射，无需手动指定覆盖度阈值。

### 比对模式说明

| 模式 | 适用场景 | 序列分歧度 |
|------|----------|-----------|
| `asm5` | 近缘物种/同一物种不同品系 | ~0.1% |
| `asm10` | 同属物种 | ~1% |
| `asm20` | 远缘物种 | ~5% |

## 输入文件格式

### 参考基因组文件

标准FASTA格式，染色体名使用标准命名（如Chr1, Chr2等）：

```fasta
>Chr1
ATCGATCGATCGATCGATCGATCGATCG...
>Chr2
GCTAGCTAGCTAGCTAGCTAGCTAGCTA...
```

### 待重命名基因组文件

标准FASTA格式，染色体名使用group命名（如group1, group2等）：

```fasta
>group1
ATCGATCGATCGATCGATCGATCGATCG...
>group2
GCTAGCTAGCTAGCTAGCTAGCTAGCTA...
```

## 输出结果

### 输出目录结构

```
chr_rename_output/
├── alignment.paf              # minimap2比对结果
├── chromosome_mapping.tsv     # 染色体映射关系表
├── rename_summary.txt         # 重命名汇总报告
├── renamed_genome.fa          # 重命名后的基因组文件
└── chr_rename.log             # 运行日志
```

### 输出文件说明

#### chromosome_mapping.tsv

染色体映射关系表，包含以下列：

- Query_Name: 原始染色体名（如group1）
- Target_Name: 目标染色体名（如Chr1）
- Query_Length: query序列长度
- Target_Length: target序列长度
- Coverage: 覆盖度（比对长度/query长度）
- Identity: 序列一致性（匹配数/比对长度）

示例：

```
Query_Name	Target_Name	Query_Length	Target_Length	Coverage	Identity
group1	Chr1	50000000	50000000	95.20%	99.85%
group2	Chr2	45000000	45000000	88.50%	99.80%
```

#### rename_summary.txt

重命名汇总报告，包含：

1. 基本信息：输入文件、比对参数
2. 统计信息：成功映射数量、一对多映射数量
3. 一对多映射详情：展示碎片化组装情况
4. 警告信息：嵌合组装、未映射序列等

#### renamed_genome.fa

重命名后的基因组FASTA文件，染色体名已转换为标准命名。

## 使用示例

### 示例1：基本重命名流程

```bash
# 对基因组挂载结果进行染色体重命名（使用自动迭代映射策略）
biopytools chr_rename \
    -r reference_genome.fa \
    -q assembled_genome.fa \
    -o renamed_output
```

### 示例2：远缘物种比对

```bash
# 对远缘物种使用asm20模式
biopytools chr_rename \
    -r reference.fa \
    -q query.fa \
    -x asm20 \
    -i 0.85 \
    -o distant_rename
```

### 示例3：指定minimap2路径

```bash
# 使用自定义minimap2路径
biopytools chr_rename \
    -r reference.fa \
    -q query.fa \
    -a /usr/local/bin/minimap2 \
    -o output
```

## 迭代映射策略详解

工具采用智能迭代映射策略，从高覆盖度阈值到低阈值自动尝试映射，确保优先使用高质量比对结果。

### 工作原理

1. **第1轮（覆盖度 >= 90%）**: 只映射高质量序列
2. **第2轮（覆盖度 >= 80%）**: 映射剩余未映射的高质量序列
3. **第3轮（覆盖度 >= 70%）**: 映射剩余未映射的较高质量序列
4. **第4轮（覆盖度 >= 60%）**: 映射剩余未映射的中等质量序列
5. **第5轮（覆盖度 >= 50%）**: 映射剩余未映射的中低质量序列
6. **第6轮（覆盖度 >= 40%）**: 映射剩余未映射的低质量序列
7. **第7轮（覆盖度 >= 30%）**: 映射剩余未映射的更低质量序列
8. **第8轮（覆盖度 >= 20%）**: 映射剩余未映射的最低质量序列

### 核心优势

- **质量优先**: 高覆盖度的序列优先映射，避免被低质量比对覆盖
- **自动化**: 无需手动指定覆盖度阈值，工具自动寻找最佳阈值
- **全面性**: 即使是低覆盖度序列（如部分contig）也能成功映射
- **透明性**: 汇总报告中详细记录每轮的映射情况

### 示例输出

```
第4轮|Round 4: 阈值|Threshold >= 60%, 映射|Mapped 3 个序列|sequences
  (group1 -> Chr1: 64.18%, group3 -> Chr7: 66.92%, group8 -> Chr8: 60.25%)

第5轮|Round 5: 阈值|Threshold >= 50%, 映射|Mapped 3 个序列|sequences
  (group4 -> Chr6: 55.59%, group5 -> Chr4: 58.33%, group7 -> Chr3: 57.41%)

第6轮|Round 6: 阈值|Threshold >= 40%, 映射|Mapped 2 个序列|sequences
  (group2 -> Chr5: 44.06%, group6 -> Chr2: 43.81%)
```

### 与固定阈值策略的对比

| 策略 | 优点 | 缺点 |
|------|------|------|
| **迭代策略（默认）** | 自动优化、质量优先、全面覆盖 | 需要多轮迭代（速度略慢） |
| 固定高阈值（如0.8） | 速度快、只保留高质量 | 可能遗漏低质量序列 |
| 固定低阈值（如0.3） | 能映射所有序列 | 高质量序列可能被低质量比对覆盖 |

## 映射关系处理

### 一对多映射（多个group映射到同一Chr）

**含义**：可能表示基因组碎片化组装，一个参考染色体被组装成了多个contig/scaffold

**处理方式**：保留所有映射关系，全部映射到同一目标染色体

**输出提示**：
```
INFO: 多个group映射到|Multiple groups mapped to Chr1:
  group1, group2, group5
```

### 多对一映射（一个group映射到多个Chr）

**含义**：可能表示嵌合组装或错误组装

**处理方式**：选择最佳比对（覆盖度最高、比对长度最长）作为目标

**输出警告**：
```
WARNING: group5 可能是嵌合/错误组装|may be chimeric/misassembled
  最佳比对|Best alignment: Chr3 (覆盖度|coverage: 85%)
  所有比对目标|All targets: Chr3, Chr7, Chr12
```

### 未映射序列

**含义**：在参考基因组中找不到对应的序列

**处理方式**：保留原始名称，不进行重命名

**输出警告**：
```
WARNING: 未映射|Unmapped: group15，保留原名|keeping original name
```

## 系统要求

### 依赖软件

- **minimap2** (版本 2.24 或更新)
  - 下载地址: https://github.com/lh3/minimap2
  - 安装: `conda install -c bioconda minimap2` 或从源码编译
- **Python** (版本 3.7+)
- **Python包**:
  - `click` - 命令行界面（已包含在biopytools中）

### 安装依赖软件

```bash
# 使用conda安装minimap2
conda install -c bioconda minimap2

# 或从源码编译
git clone https://github.com/lh3/minimap2.git
cd minimap2 && make
sudo cp minimap2 /usr/local/bin/
```

### 硬件建议

- **CPU**: 多核处理器（推荐4核以上）
- **RAM**: 最少4GB（大基因组推荐16GB以上）
- **存储**: 预留基因组文件大小3倍的磁盘空间

## 注意事项

1. **参考基因组质量**：参考基因组应使用高质量、染色体级别的组装
2. **命名规范**：参考基因组使用标准命名（如Chr1, Chr2），query使用group命名
3. **比对模式选择**：根据物种间遗传距离选择合适的preset模式
4. **覆盖度阈值**：默认0.5可根据实际情况调整，过低可能导致错误映射
5. **一致性阈值**：默认0.9适用于近缘物种，远缘物种可适当降低

## 故障排除

### 常见问题

**Q: "minimap2: command not found" 错误**

```bash
# 安装minimap2
conda install -c bioconda minimap2

# 或指定minimap2路径
biopytools chr_rename ... --minimap2-path /path/to/minimap2
```

**Q: 大量序列未映射**

```bash
# 检查并降低覆盖度阈值
biopytools chr_rename ... --min-coverage 0.3

# 或降低一致性阈值
biopytools chr_rename ... --min-identity 0.85

# 或更改比对模式
biopytools chr_rename ... --preset asm10
```

**Q: 出现大量一对多映射**

- 这是正常现象，可能表示query基因组较为碎片化
- 查看rename_summary.txt了解详细信息
- 如果是高质量基因组组装，考虑检查覆盖度阈值是否过低

**Q: 映射关系不正确**

```bash
# 提高过滤阈值
biopytools chr_rename ... \
    --min-coverage 0.8 \
    --min-identity 0.98

# 手动检查PAF文件
cat chr_rename_output/alignment.paf | less
```

## 结果解读

### 检查映射质量

1. 查看chromosome_mapping.tsv中的Coverage和Identity列
2. 高质量映射：Coverage > 80%, Identity > 95%
3. 中等质量映射：Coverage 50-80%, Identity 90-95%
4. 低质量映射：Coverage < 50% 或 Identity < 90%（需要人工检查）

### 验证重命名结果

```bash
# 统计重命名前后的染色体数量
grep "^>" query.fa | wc -l  # 原始数量
grep "^>" chr_rename_output/renamed_genome.fa | wc -l  # 重命名后数量

# 检查是否所有序列都被重命名
grep "^>" chr_rename_output/renamed_genome.fa | grep "group"
# 如果有输出，说明这些序列未找到映射，保留了原名
```

## 相关资源

- [minimap2官方文档](https://github.com/lh3/minimap2)
- [PAF格式规范](https://github.com/lh3/minimap2/blob/master/PAF.md)
- [比较基因组学最佳实践](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5411875/)

## 许可证

本项目采用MIT许可证 - 详见 [LICENSE](LICENSE) 文件

## 引用信息

如果在学术研究中使用此工具，请引用minimap2相关文献：

```
Li, H. (2018).
Minimap2: pairwise alignment for nucleotide sequences.
Bioinformatics, 34(18), 3094-3100.
```
