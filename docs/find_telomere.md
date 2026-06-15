# 🧬 端粒识别分析模块

**专业的基因组端粒重复序列识别和可视化工具 | Professional Telomeric Repeat Identification and Visualization Tool**

## 📖 功能概述 | Overview

端粒识别分析模块是一个强大的端粒重复序列识别和可视化工具，基于 tidk (Telomere Identification toolKit) 构建。该模块提供从端粒重复序列探索、查找、搜索到可视化的完整流程，支持多种分析模式和灵活的参数配置，适用于各种基因组端粒研究。

## ✨ 主要特性 | Key Features

- **🔬 自动探索模式**: 自动识别基因组中的端粒重复序列单元
- **🎯 智能查找模式**: 基于分类群的已知端粒重复序列进行查找
- **🔍 自定义搜索**: 使用自定义序列搜索端粒重复
- **📊 可视化绘图**: 将端粒分布结果绘制为 SVG 矢量图
- **🌐 广泛物种支持**: 涵盖昆虫、脊椎动物、植物等多个类群
- **⚙️ 高度可配置**: 灵活的参数设置和输出格式选择
- **📝 详细日志记录**: 完整的处理过程日志和错误追踪
- **🚀 高效处理**: 优化的处理流程，支持大规模基因组数据

## 🚀 快速开始 | Quick Start

### 基本用法 | Basic Usage

```bash
# 探索模式 - 自动识别端粒重复序列
biopytools find_telomere \
    -i genome.fa \
    -m explore \
    -o explore_results

# 查找模式 - 根据分类群查找端粒
biopytools find_telomere \
    -i genome.fa \
    -m find \
    -c Mammalia \
    -o find_results

# 搜索模式 - 使用自定义序列搜索
biopytools find_telomere \
    -i genome.fa \
    -m search \
    -s TTAGGG \
    -o search_results

# 绘图模式 - 可视化端粒分布
biopytools find_telomere \
    -m plot \
    -t telomere_telomeric_repeat_windows.tsv \
    -o plot_results
```

### 高级用法 | Advanced Usage

```bash
# 探索模式 - 自定义参数
biopytools find_telomere \
    -i genome.fa \
    -m explore \
    --explore-min 6 \
    --explore-max 15 \
    --explore-threshold 150 \
    --explore-distance 0.02 \
    -o custom_explore

# 查找模式 - 调整窗口大小
biopytools find_telomere \
    -i genome.fa \
    -m find \
    -c Coleoptera \
    --window 5000 \
    -o coleoptera_telomeres

# 搜索模式 - 输出为 bedgraph 格式
biopytools find_telomere \
    -i genome.fa \
    -m search \
    -s TTAGGG \
    --format bedgraph \
    -o bedgraph_results

# 绘图模式 - 自定义图像尺寸
biopytools find_telomere \
    -m plot \
    -t results/telomere_windows.tsv \
    --plot-width 1500 \
    --plot-height 300 \
    --plot-fontsize 14 \
    -o custom_plot
```

## 📋 参数说明 | Parameters

### 必需参数 | Required Parameters

| 参数 | 描述 | 示例 |
|------|------|------|
| `-i, --genome` | 基因组序列文件路径 (FASTA 格式) | `-i genome.fa` |

### 分析模式选择 | Analysis Mode Selection

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-m, --mode` | `find` | 🎯 分析模式: explore(探索) / find(查找) / search(搜索) / plot(绘图) |

### 输出参数 | Output Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-o, --output-dir` | `./telomere_output` | 📁 输出目录路径 |
| `-p, --prefix` | `telomere` | 📝 输出文件前缀 |

### Explore 模式参数 | Explore Mode Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--explore-min` | `5` | 🔬 最小重复序列长度 |
| `--explore-max` | `12` | 🔬 最大重复序列长度 |
| `--explore-threshold` | `100` | 🎯 重复序列出现阈值 |
| `--explore-distance` | `0.01` | 📍 染色体末端搜索距离比例 (0-0.5) |

### Find 模式参数 | Find Mode Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-c, --clade` | `必需` | 🧬 分类群名称 (如: Mammalia, Coleoptera, Poales) |
| `-w, --window` | `10000` | 📏 计算端粒重复数量的窗口大小 |

### Search 模式参数 | Search Mode Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-s, --search-string` | `必需` | 🔍 搜索的 DNA 字符串 (如: TTAGGG) |
| `--format` | `tsv` | 📊 输出格式: tsv / bedgraph |

### Plot 模式参数 | Plot Mode Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-t, --tsv` | `必需` | 📈 输入的 TSV 文件路径 |
| `--plot-height` | `200` | 📐 子图高度 (像素) |
| `--plot-width` | `1000` | 📐 图像总宽度 (像素) |
| `--plot-fontsize` | `12` | 🔤 字体大小 |
| `--plot-strokewidth` | `2` | ✏️ 线条宽度 |

### 软件配置 | Software Configuration

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--tidk-path` | `/share/org/YZWL/yzwl_lixg/miniforge3/envs/tidk_v.0.2.65/bin/tidk` | 🛠️ tidk 软件安装路径 |

### 日志选项 | Logging Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-v, --verbose` | `False` | 📢 详细输出模式 |
| `--log-file` | `None` | 📝 日志文件路径 |

### 其他选项 | Other Options

| 参数 | 描述 |
|------|------|
| `--print-clades` | 📋 打印所有支持的分类群列表 |

## 📁 输入文件格式 | Input File Formats

### 基因组序列文件 | Genome Sequence File

标准 FASTA 格式的基因组序列文件：

```fasta
>chromosome1
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>chromosome2
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
```

**文件要求**:
- 标准 FASTA 格式
- 可以是压缩或未压缩
- 支持完整基因组或部分 scaffold/contig

### TSV 文件 (用于 Plot 模式) | TSV File (for Plot Mode)

由 find 或 search 模式生成的 TSV 文件：

```tsv
chromosome	position	start	end	telomeric_repeat_count
chromosome1	0	10000	10000	45
chromosome1	10000	20000	20000	2
chromosome2	0	10000	10000	38
```

## 📂 输出文件 | Output Files

### Explore 模式输出 | Explore Mode Output

| 文件 | 描述 |
|------|------|
| `{prefix}_explore.tsv` | 端粒重复序列探索结果 |

输出示例：
```tsv
repeat	reverse_complement	count
TTAGGG	AACCCC	1523
TTAGG	CCTAA	845
```

### Find 模式输出 | Find Mode Output

| 文件 | 描述 |
|------|------|
| `{prefix}_telomeric_repeat_windows.tsv` | 端粒重复序列窗口统计 |
| `{prefix}_telomeric_repeat_counts.tsv` | 端粒重复序列总数统计 |

### Search 模式输出 | Search Mode Output

| 文件 | 描述 |
|------|------|
| `{prefix}_search_windows.tsv` | 搜索序列窗口统计 |
| `{prefix}_search_counts.tsv` | 搜索序列总数统计 |

如果选择 bedgraph 格式：
```bedgraph
chromosome1	0	10000	45
chromosome1	10000	20000	2
```

### Plot 模式输出 | Plot Mode Output

| 文件 | 描述 |
|------|------|
| `{prefix}_plot.svg` | SVG 矢量格式的端粒分布图 |

## 🌐 支持的分类群 | Supported Clades

### 脊椎动物 | Vertebrates

| 分类群 | 端粒重复序列 | 描述 |
|--------|------------|------|
| `Mammalia` | TTAGGG | 哺乳动物 |
| `Carnivora` | TTAGGG | 食肉目 |
| `Rodentia` | TTAGGG | 啮齿目 |
| `Chiroptera` | TTAGGG | 翼手目 (蝙蝠) |
| `Aves` | TTAGGG | 鸟类 |
| `Accipitriformes` | TTAGGG | 鹰形目 |
| `Anura` | TTAGGG | 无尾目 (青蛙) |
| `Caprimulgiformes` | TTAGGG | 夜鹰目 |
| `Actinopterygii` | TTAGGG | 辐鳍鱼纲 |
| `Perciformes` | TTAGGG | 鲈形目 |
| `Salmoniformes` | TTAGGG | 鲑形目 |
| `Cypriniformes` | TTAGGG | 鲤形目 |
| `Labriformes` | TTAGGG | 隆头鱼目 |
| `Syngnathiformes` | TTAGGG | 海龙鱼目 |
| `Carangiformes` | TTAGGG | 鲹形目 |
| `Pleuronectiformes` | TTAGGG | 鲽形目 |
| `Carcharhiniformes` | TTAGGG | 真鲨目 |

### 昆虫 | Insects

| 分类群 | 端粒重复序列 | 描述 |
|--------|------------|------|
| `Coleoptera` | TTAGG | 鞘翅目 |
| `Hymenoptera` | TTAGG | 膜翅目 (蚂蚁、蜜蜂) |
| `Lepidoptera` | TTAGG | 鳞翅目 (蝴蝶、蛾) |
| `Diptera` | TTAGG | 双翅目 (苍蝇、蚊子) |
| `Hemiptera` | TTAGG | 半翅目 |
| `Orthoptera` | TTAGG, TTAGGG | 直翅目 |
| `Odonata` | TTAGG | 蜻蜓目 |
| `Plecoptera` | TTAGG | 襀翅目 |
| `Trichoptera` | TTAGG | 毛翅目 |
| `Symphypleona` | TTAGG | 圆跳虫科 |

### 植物 | Plants

| 分类群 | 端粒重复序列 | 描述 |
|--------|------------|------|
| `Arabidopsis` | TTTAGGG | 拟南芥属 |
| `Poales` | TTTAGGG | 禾本目 (包括禾本科) |
| `Rosales` | TTTAGGG | 蔷薇目 |
| `Fabales` | TTTAGGG | 豆目 |
| `Malpighiales` | TTTAGGG | 金虎尾目 |
| `Myrtales` | TTTAGGG | 桃金娘目 |
| `Sapindales` | TTTAGGG | 无患子目 |
| `Caryophyllales` | TTTAGGG | 石竹目 |
| `Asterales` | TTTAGGG | 菊目 |
| `Lamiales` | TTTAGGG | 唇形目 |
| `Solanales` | TTTAGGG | 茄目 |
| `Apiales` | TTTAGGG | 伞形目 |
| `Fagales` | TTTAGGG | 壳斗目 |
| `Buxales` | TTTAGGG | 黄杨目 |
| `Hypnales` | TTTAGGG | 灰藓目 |
| `Chlamydomonadales` | TTTAGGG | 衣藻目 |

### 其他生物 | Other Organisms

| 分类群 | 端粒重复序列 | 描述 |
|--------|------------|------|
| `Nematoda` | TTAGGC | 线虫门 |
| `Arachnida` | TTAGGG | 蛛形纲 |
| `Crustacea` | TTAGGG | 甲壳类 |
| `Fungi` | TTAGGG | 真菌 |

**查看所有支持的分类群**:
```bash
biopytools find_telomere --print-clades
```

## 🔬 分析模式详解 | Analysis Mode Details

### 1. Explore 模式 | 探索模式

**功能** | Purpose: 自动识别基因组中的端粒重复序列单元

**适用场景** | Use cases:
- 未知端粒重复序列的物种
- 想要验证端粒重复序列
- 新物种基因组分析

**工作原理** | How it works:
1. 在染色体末端 (默认 1% 区域) 搜索重复序列
2. 测试不同长度的 kmer (默认 5-12 bp)
3. 统计连续出现的重复序列
4. 输出最可能的端粒重复序列及其反向互补序列

**示例输出** | Example output:
```
repeat	reverse_complement	count
TTAGGG	AACCCC	1523
```
表示端粒重复序列为 TTAGGG，反向互补为 AACCCC，出现 1523 次。

### 2. Find 模式 | 查找模式

**功能** | Purpose: 基于分类群的已知端粒重复序列查找端粒位置

**适用场景** | Use cases:
- 已知物种的分类归属
- 快速定位端粒位置
- 染色体级别基因组组装质量评估

**工作原理** | How it works:
1. 根据分类群选择对应的端粒重复序列
2. 在基因组中以滑动窗口方式搜索
3. 统计每个窗口内的重复序列数量
4. 生成端粒分布的统计表

**输出文件说明** | Output files:
- `*_telomeric_repeat_windows.tsv`: 每个窗口的详细统计
- `*_telomeric_repeat_counts.tsv`: 每条染色体的总统计

### 3. Search 模式 | 搜索模式

**功能** | Purpose: 使用自定义序列搜索端粒

**适用场景** | Use cases:
- 已知端粒重复序列
- 分类群不在数据库中
- 研究非典型端粒序列

**工作原理** | How it works:
1. 使用用户提供的序列搜索基因组
2. 统计窗口内的匹配数量
3. 支持 TSV 和 bedgraph 输出格式

**与 Find 模式的区别** | Difference from Find mode:
- Find: 使用内置数据库的分类群-序列映射
- Search: 使用用户自定义的搜索序列

### 4. Plot 模式 | 绘图模式

**功能** | Purpose: 将端粒分布结果可视化

**适用场景** | Use cases:
- 直观展示端粒分布
- 识别染色体末端
- 评估基因组组装质量
- 生成发表级别的图表

**工作原理** | How it works:
1. 读取 find 或 search 生成的 TSV 文件
2. 为每条染色体绘制子图
3. X 轴为染色体位置，Y 轴为端粒重复序列数量
4. 输出高分辨率 SVG 矢量图

**图像特点** | Image features:
- 矢量格式，可无损缩放
- 每条染色体独立子图
- 端粒区域显示明显的峰值
- 适合论文发表

## 📊 应用场景 | Applications

### 1. 基因组组装质量评估 | Genome Assembly Quality Assessment

端粒是染色体末端的保护结构，完整的基因组组装应该在染色体末端检测到端粒重复序列的富集。

```bash
# 查找端粒并绘图
biopytools find_telomere \
    -i genome.fa \
    -m find \
    -c Mammalia \
    -o assembly_quality_check

biopytools find_telomere \
    -m plot \
    -t assembly_quality_check/telomere_telomeric_repeat_windows.tsv \
    -o assembly_quality_check
```

**质量评估标准**:
- ✅ 优质组装: 每条染色体两端都有明显的端粒峰
- ⚠️ 中等组装: 大部分染色体有端粒峰
- ❌ 低质组装: 端粒峰缺失或不明显

### 2. 新物种端粒序列鉴定 | Telomere Sequence Identification in New Species

对于新测序的物种，可能不知道其端粒重复序列类型，可以使用 explore 模式自动识别。

```bash
biopytools find_telomere \
    -i new_species_genome.fa \
    -m explore \
    --explore-min 5 \
    --explore-max 12 \
    -o telomere_identification
```

### 3. 比较基因组学 | Comparative Genomics

比较不同物种或品种的端粒特征：

```bash
# 物种 A
biopytools find_telomere -i species_a.fa -m find -c Mammalia -o species_a
biopytools find_telomere -m plot -t species_a/telomere_windows.tsv -o species_a_plot

# 物种 B
biopytools find_telomere -i species_b.fa -m find -c Mammalia -o species_b
biopytools find_telomere -m plot -t species_b/telomere_windows.tsv -o species_b_plot
```

### 4. 染色体进化研究 | Chromosome Evolution Studies

通过分析端粒分布模式，研究染色体进化事件：

```bash
biopytools find_telomere \
    -i genome.fa \
    -m find \
    -c Coleoptera \
    --window 5000 \
    -o chromosome_evolution
```

### 5. 端粒长度变异分析 | Telomere Length Variation Analysis

对于不同组织、发育阶段或处理条件的样本：

```bash
# 对照组
biopytools find_telomere -i control_genome.fa -m find -c Mammalia -o control
# 处理组
biopytools find_telomere -i treated_genome.fa -m find -c Mammalia -o treated
```

## ⚠️ 注意事项 | Notes

### 1. 输入文件要求 | Input File Requirements

- **格式**: 必须是标准 FASTA 格式
- **质量**: 建议使用染色体级别的基因组组装
- **完整性**: Contig/N scaffolding 可能影响端粒检测

### 2. 参数选择建议 | Parameter Selection Recommendations

#### Explore 模式
- `--explore-distance`: 染色体级别用 0.01，scaffold 级别可适当增大
- `--explore-threshold`: 基因组大时增大阈值 (避免假阳性)

#### Find/Search 模式
- `--window`: 大基因组用大窗口 (10-50kb)，小基因组用小窗口 (1-5kb)
- 分类群选择: 尽量选择最精确的分类群

#### Plot 模式
- `--plot-width`: 根据染色体数量调整，数量多则增大宽度
- `--plot-height`: 一般 200-300 像素即可

### 3. 结果解读 | Result Interpretation

**正常端粒分布**:
- 染色体两端有明显的峰
- 峰值从末端向内部逐渐降低
- 峰宽通常在 100-500kb

**异常端粒分布**:
- 染色体中间出现峰 (可能是染色体融合)
- 末端无峰 (组装不完整或端粒丢失)
- 全基因组低水平分布 (假阳性或污染)

### 4. 性能优化 | Performance Optimization

- **内存需求**: 主要取决于基因组大小，一般 10-50 GB 足够
- **运行时间**: 染色体级别基因组 10-30 分钟，scaffold 级别更长
- **并行处理**: 目前不支持多线程，可通过分割基因组加速

### 5. 常见问题 | Common Issues

#### Q1: 探索模式没有找到端粒序列？
**可能原因**:
- 染色体末端序列质量差
- `--explore-distance` 设置太小
- 基因组是 scaffold/contig 水平，缺少完整末端

**解决方案**:
- 尝试增大 `--explore-distance` 到 0.5 (处理完整序列)
- 使用 search 模式手动测试已知序列

#### Q2: 找到多个可能的端粒序列？
**可能原因**:
- 存在多种端粒重复序列类型
- 基因组中存在其他串联重复序列

**解决方案**:
- 选择计数最多的序列
- 查看文献确认该物种的端粒类型
- 检查是否为染色体末端富集

#### Q3: 绘图时某些染色体没有端粒峰？
**可能原因**:
- 该染色体组装不完整
- 该染色体实际上是 scaffold (缺少末端)
- 端粒序列在该物种中退化

**解决方案**:
- 检查基因组组装报告
- 查看原始序列确认是否有端粒序列
- 考虑使用 explore 模式重新鉴定

## 📚 参考文献 | References

1. **tidk 原始文献** | Original tidk paper:
   Brown MR, de La Rosa PG, Blaxter M. tidk: a toolkit to rapidly identify telomeric repeats from genomic datasets. Bioinformatics. 2025;btaf049. doi: 10.1093/bioinformatics/btaf049

2. **端粒数据库** | Telomere repeat database:
   https://github.com/tolkit/a-telomeric-repeat-database

3. **软件主页** | Software homepage:
   https://github.com/tolkit/telomeric-identifier

4. **Bioconda 安装** | Bioconda installation:
   https://bioconda.github.io/recipes/tidk/README.html

## 📝 版本信息 | Version Information

- **模块版本 | Module Version**: 1.0.0
- **tidk 版本 | tidk Version**: 0.2.65
- **最后更新 | Last Updated**: 2026-01-02

## 🤝 技术支持 | Technical Support

如有问题或建议，请联系：
For questions or suggestions, please contact:

- **GitHub Issues**: https://github.com/anthropics/biopytools/issues
- **Email**: yzwl_lixg@org.YZWL.share

---

**生成时间 | Generated**: 2026-01-02
**文档版本 | Doc Version**: 1.0.0
