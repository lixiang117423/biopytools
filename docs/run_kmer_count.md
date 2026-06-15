# 🧬 K-mer丰度分析工具 | K-mer Abundance Analysis Tool

[![Python](https://img.shields.io/badge/Python-3.7%2B-blue.svg)](https://python.org)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Jellyfish](https://img.shields.io/badge/Jellyfish-2.2%2B-orange.svg)](https://github.com/gmarcais/Jellyfish)

使用Jellyfish进行高效的k-mer计数和丰度分析的Python工具包，支持大规模基因组数据处理。

## 🚀 功能特点 | Features

- **🔍 自动文件识别**: 智能识别双末端FASTQ文件（支持.fq.gz、.fastq等格式）
- **📊 批量样本处理**: 同时处理多个样本，自动合并结果矩阵
- **🧬 灵活输出格式**: 支持BED文件注释或k-mer ID格式输出
- **🪟 滑动窗口分析**: 基于基因组区间的k-mer存在比例分析
- **🔢 二进制矩阵**: 生成0/1存在缺失矩阵，便于下游分析
- **⚡ 高性能处理**: 基于Jellyfish的多线程k-mer计数
- **📝 详细日志**: 完整的分析过程记录和统计报告

## 📦 安装依赖 | Dependencies

### 必需软件

```bash
# 安装Jellyfish
# Ubuntu/Debian
sudo apt-get install jellyfish

# 或使用Conda
conda install -c bioconda jellyfish

# 或从源码编译
wget https://github.com/gmarcais/Jellyfish/releases/download/v2.3.0/jellyfish-2.3.0.tar.gz
tar -xzf jellyfish-2.3.0.tar.gz
cd jellyfish-2.3.0
./configure && make && sudo make install
```

### Python依赖

```bash
pip install pandas numpy pathlib
```

## 🖥️ 使用方法 | Usage

### 基本命令格式

```bash
run_kmer_count -i <input_dir> -p <pattern> -k <kmer_lib> -o <output_dir> [options]
```

### 快速开始

```bash
# 基础k-mer丰度分析
run_kmer_count -i ./fastq_files -p "*_1.clean.fq.gz" -k kmers.fasta -o ./results

# 带BED注释的完整分析
run_kmer_count -i ./fastq_files -p "*_1.clean.fq.gz" -k kmers.fasta -b kmers.bed -o ./results

# 高性能分析（多线程 + 大内存）
run_kmer_count -i ./fastq_files -p "*_1.clean.fq.gz" -k kmers.fasta -o ./results -t 32 -s 10000M

# 滑动窗口分析
run_kmer_count -i ./fastq_files -p "*_1.clean.fq.gz" -k kmers.fasta -b kmers.bed -w 500000 --step-size 100000 -o ./results
```

## 📋 参数说明 | Parameters

### 必需参数

| 参数 | 描述 | 示例 |
|------|------|------|
| `-i, --input` | 📁 FASTQ文件输入目录 | `./fastq_files` |
| `-p, --pattern` | 🔍 文件模式（支持通配符*） | `"*_1.clean.fq.gz"` |
| `-k, --kmer-lib` | 🧬 k-mer库文件(FASTA格式) | `kmers.fasta` |
| `-o, --output` | 📂 输出目录 | `./results` |

### 可选参数

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-b, --bed-file` | None | 📋 BED注释文件路径 |
| `-m, --kmer-size` | 51 | 📏 k-mer长度 |
| `-s, --hash-size` | 1000M | 🗂️ Jellyfish哈希表大小 |
| `-t, --threads` | 8 | 🧵 线程数 |
| `-w, --window-size` | 500000 | 🪟 滑动窗口大小(bp) |
| `--step-size` | window-size/5 | 👣 滑动窗口步长(bp) |
| `-C, --canonical` | False | 🔄 统计正向和反向互补链 |
| `--keep-temp` | False | 💾 保留临时文件 |
| `--keep-binary` | False | 🔢 保留0/1存在缺失矩阵 |
| `--jellyfish-path` | jellyfish | 🐙 Jellyfish程序路径 |
| `-v, --verbose` | False | 📝 详细输出模式 |

## 📁 输入文件格式 | Input File Formats

### 1. FASTQ文件命名规范

支持的文件命名模式：
```
Sample1_1.clean.fq.gz, Sample1_2.clean.fq.gz
Sample2_R1.fastq, Sample2_R2.fastq
OV8-105_1.fq.gz, OV8-105_2.fq.gz
```

使用通配符指定模式：
```bash
--pattern "*_1.clean.fq.gz"  # 匹配 Sample*_1.clean.fq.gz 和 Sample*_2.clean.fq.gz
--pattern "*_R1.fastq"       # 匹配 Sample*_R1.fastq 和 Sample*_R2.fastq
```

### 2. K-mer库文件 (FASTA格式)

```fasta
>OV12_1_51
ctctaaaccctaaacccaaccctaaaccctaaaccctaaaccccaaacctc
>OV12_2_52
tctaaaccctaaacccaaccctaaaccctaaaccctaaaccccaaacctcc
>OV12_3_53
ctaaaccctaaacccaaccctaaaccctaaaccctaaaccccaaacctcct
```

### 3. BED注释文件 (可选)

```bed
OV12    0       51      ctctaaaccctaaacccaaccctaaaccctaaaccctaaaccccaaacctc     0       .
OV12    1       52      tctaaaccctaaacccaaccctaaaccctaaaccctaaaccccaaacctcc     0       .
OV12    2       53      ctaaaccctaaacccaaccctaaaccctaaaccctaaaccccaaacctcct     0       .
```

## 📊 输出结果 | Output Results

### 主要输出文件

| 文件名 | 描述 | 生成条件 |
|--------|------|----------|
| `kmer_abundance_matrix.tsv` | 📈 k-mer丰度矩阵 | 总是生成 |
| `sliding_window_analysis.tsv` | 🪟 滑动窗口分析 | 提供BED文件时 |
| `kmer_binary_matrix.tsv` | 🔢 0/1存在缺失矩阵 | `--keep-binary`时 |
| `analysis_summary.txt` | 📋 统计摘要报告 | 总是生成 |
| `kmer_count.log` | 📝 详细运行日志 | 总是生成 |

### 输出格式示例

#### 不使用BED文件时：
```tsv
kmer_id                 kmer                                                    Sample1    Sample2    Sample3
OV12_1_51    CTCTAAACCCTAAACCCAACCCTAAACCCTAAACCCTAAACCCCAAACCTC    15         8          0
OV12_2_52    TCTAAACCCTAAACCCAACCCTAAACCCTAAACCCTAAACCCCAAACCTCC    12         5          3
```

#### 使用BED文件时：
```tsv
chr    start    end    kmer                                                    tmp_1    tmp_2    Sample1    Sample2    Sample3
OV12   0        51     CTCTAAACCCTAAACCCAACCCTAAACCCTAAACCCTAAACCCCAAACCTC  0        .        15         8          0
OV12   1        52     TCTAAACCCTAAACCCAACCCTAAACCCTAAACCCTAAACCCCAAACCTCC  0        .        12         5          3
```

#### 滑动窗口分析结果：
```tsv
chr    start     end       kmer_count    Sample1    Sample2    Sample3
OV12   0         500000    1234          0.8567     0.7234     0.5678
OV12   100000    600000    1456          0.9012     0.8345     0.6789
```

## 🏗️ 模块结构 | Module Structure

```
biopytools/kmer_count/
├── __init__.py           # 🏠 包初始化
├── config.py             # ⚙️ 配置管理
├── utils.py              # 🔧 工具函数和日志
├── file_processor.py     # 📁 文件识别和处理
├── jellyfish_processor.py # 🐙 Jellyfish操作
├── data_processor.py     # 📊 数据解析和合并
├── window_analyzer.py    # 🪟 滑动窗口分析
├── main.py               # 🚀 主程序逻辑
└── cli.py                # 🖥️ 命令行接口
```

## 🔧 高级用法 | Advanced Usage

### 大规模数据处理

```bash
# 处理TB级数据，使用大内存和多线程
run_kmer_count \
    -i /data/large_dataset \
    -p "*_R1.fastq.gz" \
    -k large_kmer_library.fasta \
    -o /results/large_analysis \
    -t 64 \
    -s 50000M \
    --verbose
```

### 精细化窗口分析

```bash
# 小窗口高分辨率分析
run_kmer_count \
    -i ./samples \
    -p "*_1.clean.fq.gz" \
    -k kmers.fasta \
    -b regions.bed \
    -w 50000 \
    --step-size 10000 \
    -C \
    -o fine_scale_analysis
```

### 调试模式

```bash
# 保留所有中间文件用于调试
run_kmer_count \
    -i ./samples \
    -p "*_1.fq.gz" \
    -k kmers.fasta \
    -o debug_results \
    --keep-temp \
    --keep-binary \
    --verbose
```

## ⚡ 性能优化建议 | Performance Optimization

### 内存设置

根据数据大小调整`--hash-size`：
- 🔸 小数据集（<10GB）: `1000M-5000M`
- 🔸 中等数据集（10-100GB）: `5000M-20000M`
- 🔸 大数据集（>100GB）: `20000M-100000M`

### 线程配置

```bash
# 建议设置为CPU核心数的80%
-t $(nproc --all | awk '{print int($1*0.8)}')
```

### 磁盘空间

确保有足够的磁盘空间：
- 📦 解压缩FASTQ文件需要额外空间
- 🗂️ Jellyfish数据库文件可能很大
- 📊 临时文件和结果文件

## 🐛 故障排除 | Troubleshooting

### 常见问题

**Q: 为什么找不到FASTQ文件？**
```bash
# 检查文件路径和模式
ls -la /path/to/fastq/files/*_1.clean.fq.gz
```

**A: 确保：**
- 输入目录路径正确
- 文件模式匹配实际文件名
- 文件权限允许读取

**Q: 内存不足错误？**

**A: 解决方案：**
- 减少`--hash-size`参数
- 分批处理大文件
- 使用更多内存的机器

**Q: Jellyfish命令未找到？**

**A: 解决方案：**
```bash
# 检查Jellyfish是否安装
which jellyfish

# 或指定完整路径
--jellyfish-path /path/to/jellyfish
```

### 调试技巧

1. **使用详细模式**: `--verbose`
2. **保留临时文件**: `--keep-temp`
3. **检查日志文件**: `kmer_count.log`
4. **从小数据集开始测试**

## 📈 输出解读 | Result Interpretation

### 丰度值含义

- **数值**: k-mer在样本中的出现次数
- **0**: k-mer在该样本中不存在
- **>0**: k-mer存在，数值越大表示丰度越高

### 滑动窗口比例

- **值范围**: 0.0 - 1.0
- **含义**: 窗口内存在的k-mer占总k-mer的比例
- **1.0**: 窗口内所有k-mer都存在
- **0.0**: 窗口内没有k-mer存在

## 📄 引用 | Citation

如果这个工具对您的研究有帮助，请引用：

```
K-mer Abundance Analysis Tool (biopytools.kmer_count)
Author: biopytools team
Version: 1.0.0
```

## 🤝 贡献指南 | Contributing

欢迎提交Issue和Pull Request！

1. Fork 项目
2. 创建特性分支 (`git checkout -b feature/AmazingFeature`)
3. 提交更改 (`git commit -m 'Add some AmazingFeature'`)
4. 推送到分支 (`git push origin feature/AmazingFeature`)
5. 打开Pull Request

## 📜 许可证 | License

本项目采用MIT许可证 - 查看 [LICENSE](LICENSE) 文件了解详情。

## 📧 联系方式 | Contact

- 📧 Email: biopytools@example.com
- 🐛 Issues: [GitHub Issues](https://github.com/biopytools/issues)
- 💬 Discussions: [GitHub Discussions](https://github.com/biopytools/discussions)

---

**💡 提示**: 首次使用建议从小数据集开始测试，熟悉各种参数的效果后再处理大规模数据。