# GenomeScope2 基因组特征评估工具 | GenomeScope2 Genome Characteristics Evaluation Tool

## 概述 | Overview

GenomeScope2是基于k-mer分析的基因组特征评估工具，能够估算基因组大小、杂合度、重复序列比例等关键参数，为基因组组装策略制定提供重要参考。

## 功能特点 | Features

- 🧬 **基因组大小估算**: 基于k-mer频率分布估算基因组总长度
- 🧩 **杂合度分析**: 评估基因组杂合程度，判断组装难度
- 🔁 **重复序列分析**: 量化重复序列含量和比例
- 📊 **自动报告生成**: 生成详细的技术报告和通俗解读
- 🚀 **高性能处理**: 支持多线程和大数据集处理
- 📈 **结果可视化**: 自动生成k-mer分布图和模型拟合图

## 安装依赖 | Dependencies

### 必需软件 | Required Software

1. **Jellyfish**: k-mer计数工具
   ```bash
   # Ubuntu/Debian
   sudo apt-get install jellyfish

   # CentOS/RHEL
   sudo yum install jellyfish
   ```

2. **R 环境**: 运行GenomeScope R脚本
   ```bash
   # Ubuntu/Debian
   sudo apt-get install r-base

   # CentOS/RHEL
   sudo yum install R
   ```

3. **R 包依赖**: 安装GenomeScope依赖包
   ```r
   # 在R中执行
   install.packages(c("ggplot2", "dplyr", "tidyr"))
   ```

### 脚本依赖 | Script Dependencies

确保GenomeScope R脚本存在：
```bash
# 默认路径
/share/org/YZWL/yzwl_lixg/software/scripts/genomescope.R
```

## 使用方法 | Usage

### 基本用法 | Basic Usage

```bash
# 最简单的使用方式
biopytools genomescope \
    -i fastq_files/ \
    -o sample_analysis \
    -l 150
```

### 参数说明 | Parameters

#### 必需参数 | Required Parameters

| 参数 | 说明 | 示例 |
|------|------|------|
| `-i, --input-dir` | 包含FASTQ文件的输入目录 | `fastq_files/` |
| `-o, --output-prefix` | 输出文件前缀 | `sample_analysis` |
| `-l, --read-length` | 测序数据平均读长 | `150` |

#### 可选参数 | Optional Parameters

| 参数 | 默认值 | 说明 | 示例 |
|------|--------|------|------|
| `-t, --threads` | `64` | 使用的线程数 | `32` |
| `-s, --hash-size` | `10G` | Jellyfish哈希表大小 | `20G` |
| `-k, --kmer-size` | `21` | k-mer大小 | `31` |
| `-c, --max-kmer-cov` | `1000` | 最大k-mer覆盖度 | `2000` |
| `--genomescope-r-script` | 见下方 | GenomeScope R脚本路径 | `/path/to/genomescope.R` |

#### 日志控制 | Logging Control

| 参数 | 说明 | 示例 |
|------|------|------|
| `-v, --verbose` | 详细输出模式(-v: INFO, -vv: DEBUG) | `-vv` |
| `--quiet` | 静默模式(只输出ERROR) | `--quiet` |
| `--log-file` | 日志文件名 | `analysis.log` |
| `--log-dir` | 日志目录 | `./logs` |
| `--dry-run` | 模拟运行(不实际执行) | `--dry-run` |

### 使用示例 | Examples

#### 1. 小基因组分析 (<100Mb)
```bash
biopytools genomescope \
    -i small_genome_reads/ \
    -o small_genome_analysis \
    -l 150 \
    -t 16 \
    -s 5G \
    -k 21
```

#### 2. 中等基因组分析 (100-500Mb)
```bash
biopytools genomescope \
    -i medium_genome_reads/ \
    -o medium_genome_analysis \
    -l 150 \
    -t 32 \
    -s 15G \
    -k 31
```

#### 3. 大基因组分析 (>500Mb)
```bash
biopytools genomescope \
    -i large_genome_reads/ \
    -o large_genome_analysis \
    -l 150 \
    -t 64 \
    -s 50G \
    -k 31 \
    -c 2000
```

#### 4. 高杂合度基因组
```bash
biopytools genomescope \
    -i heterozygous_genome/ \
    -o het_genome_analysis \
    -l 150 \
    -k 21 \
    -c 3000 \
    -vv
```

#### 5. 详细日志和调试
```bash
biopytools genomescope \
    -i reads/ \
    -o debug_analysis \
    -l 150 \
    -vv \
    --log-file debug.log \
    --log-dir ./debug_logs
```

#### 6. 模拟运行测试
```bash
biopytools genomescope \
    -i test_reads/ \
    -o test_analysis \
    -l 150 \
    --dry-run
```

## 输出文件 | Output Files

GenomeScope2分析会生成以下输出文件：

### 主要输出文件

1. **Jellyfish输出**
   - `{prefix}.jf`: Jellyfish k-mer数据库
   - `{prefix}.histo`: k-mer频率直方图

2. **GenomeScope输出目录**: `{prefix}_genomescope_output/`
   - `summary.txt`: 分析结果摘要
   - `model.png`: 拟合模型图
   - `log2_transformed_model.png`: 对数变换模型图
   - `linear_transformed_model.png`: 线性变换模型图
   - `coverage_distribution.png`: 覆盖度分布图

3. **易读报告**
   - `{prefix}_report_readable.txt`: 通俗易懂的分析报告

### 日志文件

- `{log_dir}/{log_file}`: 详细执行日志

## 结果解读 | Results Interpretation

### 关键参数说明

#### 1. 基因组大小 (Genome Haploid Length)
- **含义**: 基因组单倍体长度估算值
- **单位**: bp (碱基对)
- **参考值**:
  - 细菌: 1-10 Mb
  - 酵母: 10-20 Mb
  - 果蝇: ~140 Mb
  - 人类: ~3,000 Mb

#### 2. 杂合度 (Heterozygosity)
- **含义**: 基因组中杂合位点的比例
- **单位**: 百分比 (%)
- **解读指南**:
  - < 0.5%: 低杂合度，适合单倍型组装
  - 0.5-1%: 较低杂合度，常规组装
  - 1-2%: 中等杂合度
  - > 2%: 高杂合度，需要特殊策略

#### 3. 重复序列比例 (Repeat Content)
- **含义**: 重复序列在基因组中的占比
- **单位**: 百分比 (%)
- **影响**: 重复序列越多，组装难度越大

#### 4. 模型拟合度 (Model Fit)
- **含义**: 数据与理论模型的吻合程度
- **单位**: 百分比 (%)
- **质量评估**:
  - > 90%: 数据质量优秀
  - 70-90%: 数据质量良好
  - < 70%: 数据质量需检查

#### 5. 测序错误率 (Read Error Rate)
- **含义**: 测序过程中碱基错误率
- **单位**: 百分比 (%)
- **质量标准**:
  - < 0.1%: 测序质量优秀
  - 0.1-0.5%: 测序质量良好
  - > 0.5%: 测序质量需关注

### 下一步建议

#### 基于杂合度的组装策略

1. **低杂合度 (<0.5%)**
   - 推荐工具: Hifiasm, Flye
   - 策略: 单倍型组装
   - 测序深度: HiFi 30-50x, ONT 50-80x

2. **中等杂合度 (0.5-2%)**
   - 推荐工具: Hifiasm (默认模式), Canu
   - 策略: 常规二倍体组装
   - 测序深度: HiFi 30-50x, ONT 60-100x

3. **高杂合度 (>2%)**
   - 推荐工具: Hifiasm (三倍体模式), Trio-binning
   - 策略: 分型组装
   - 建议数据: HiFi + Hi-C 或 亲本数据

## 性能优化 | Performance Optimization

### 内存配置建议

| 基因组大小 | 测序深度 | 哈希表大小 | 内存需求 |
|------------|----------|------------|----------|
| < 100 Mb | 50x | 5G | 8-16 GB |
| 100-500 Mb | 50x | 10-20G | 16-32 GB |
| 500 Mb-1 Gb | 50x | 30-50G | 32-64 GB |
| > 1 Gb | 50x | 50-100G | 64-128 GB |

### 线程数建议

- **小数据集**: 8-16 线程
- **中等数据集**: 16-32 线程
- **大数据集**: 32-64 线程
- **注意**: 过多线程可能导致I/O瓶颈

## 故障排除 | Troubleshooting

### 常见问题

#### 1. 内存不足
**错误**: `jellyfish: hash table size too large`
**解决方案**:
- 减少哈希表大小 (`-s` 参数)
- 使用更少的线程
- 增加系统内存

#### 2. 模型拟合度低
**可能原因**:
- 测序数据质量问题
- k-mer大小选择不当
- 样本污染
**解决方案**:
- 检查原始数据质量
- 尝试不同的k-mer大小 (15, 21, 31, 51)
- 检查是否有污染物

#### 3. 杂合度估算异常
**可能原因**:
- 二倍体样本实际上是多倍体
- 存在样本污染
- 测序覆盖度过低
**解决方案**:
- 检查样本倍性
- 验证样本纯度
- 增加测序覆盖度

#### 4. 分析运行缓慢
**优化方法**:
- 使用SSD存储
- 增加线程数
- 优化哈希表大小
- 使用压缩格式减少I/O

### 日志分析

查看详细日志了解问题原因：
```bash
# 查看执行日志
tail -f logs/genomescope.log

# 查看错误信息
grep ERROR logs/genomescope.log

# 查看警告信息
grep WARNING logs/genomescope.log
```

## 技术细节 | Technical Details

### 分析流程

1. **k-mer计数**: 使用Jellyfish统计所有k-mer的出现频率
2. **频率分布**: 生成k-mer频率直方图
3. **模型拟合**: GenomeScope拟合理论模型到实际数据
4. **参数提取**: 从拟合结果提取基因组特征参数
5. **结果验证**: 评估模型拟合度和数据质量

### 算法原理

GenomeScope基于以下假设：
- 大部分k-mer来自基因组唯一序列
- 杂合位点的k-mer覆盖度约为纯合位点的一半
- 重复序列导致k-mer覆盖度增加
- 测序错误产生低频k-mer

## 更新日志 | Changelog

### v1.0.0 (2024-12-22)
- ✨ 初始版本发布
- 🧬 完整的GenomeScope2分析流程
- 📊 自动报告生成
- 🚀 高性能处理支持
- 📝 详细日志记录

## 引用 | Citation

如果使用本工具，请引用原始GenomeScope论文：

```
Vurture GW, Sedlazeck FJ, Nattestad M, et al.
GenomeScope: fast reference-free genome profiling from short reads.
Bioinformatics. 2017;33(14):2202-2208.
```

## 联系方式 | Contact

如有问题或建议，请联系生信分析团队或提交Issue。