# GTX Joint Calling 命令生成工具

## 功能概述

`biopytools gtx-joint` 是一个专门用于生成GTX Joint Calling命令脚本的工具。它可以根据参考基因组的染色体信息，自动生成按染色体或按区间拆分的GTX joint calling命令，支持大规模多样品联合变异检测任务的并行化处理。

## 核心功能

- 📂 **自动扫描GVCF文件** - 自动发现并验证输入的GVCF文件及其索引
- 🧬 **染色体智能处理** - 从参考基因组索引读取染色体信息
- 🔄 **灵活的拆分策略** - 支持按染色体或按固定区间拆分任务
- 📝 **生成可执行脚本** - 生成可以直接运行的shell命令脚本
- 🔍 **染色体过滤** - 支持使用正则表达式过滤需要处理的染色体
- ⏰ **faketime支持** - 可选的faketime支持，解决时间相关问题

## 安装依赖

使用本工具需要预先安装以下软件：

- **GTX** - 用于Joint Calling的核心软件
- **samtools** - 用于参考基因组索引(faidx)
- **faketime** (可选) - 用于解决软件时间限制问题

### 安装示例

```bash
# 安装samtools
conda install -c bioconda samtools

# 为参考基因组创建索引（必需）
samtools faidx genome.fa
```

## 使用方法

### 基本用法 - 按染色体生成

```bash
biopytools gtx-joint \
    -g /share/org/YZWL/yzwl_lixg/software/gtx/bin/gtx \
    -r reference.fa \
    -i ./gvcf_files \
    -o ./joint_output
```

这将生成一个名为 `run_gtx_joint.sh` 的脚本，每个染色体一条命令。

### 按10M区间拆分（推荐用于大基因组）

```bash
biopytools gtx-joint \
    -g ./gtx \
    -r genome.fa \
    -i ./gvcf \
    -o ./output \
    -w 10000000
```

使用区间模式可以有效避免内存溢出问题，特别适合大型基因组项目。

### 只处理主染色体

```bash
biopytools gtx-joint \
    -g ./gtx \
    -r genome.fa \
    -i ./gvcf \
    -o ./output \
    -p "^Chr[0-9]+$"
```

使用正则表达式过滤，只处理匹配模式的染色体。

### 完整参数示例

```bash
biopytools gtx-joint \
    -g /path/to/gtx \
    -r reference.fa \
    -i ./gvcf_files \
    -o ./joint_output \
    -t 24 \
    -T /tmp \
    -s my_joint_commands.sh \
    -f "2020-01-01 00:00:00" \
    -p "^Chr" \
    -w 10000000
```

## 命令行参数

### 必需参数

| 参数 | 简写 | 说明 | 示例 |
|------|------|------|------|
| `--gtx` | `-g` | GTX可执行文件路径 | `/share/org/YZWL/yzwl_lixg/software/gtx/bin/gtx` |
| `--ref` | `-r` | 参考基因组文件路径（需要.fai索引） | `genome.fa` |
| `--input` | `-i` | GVCF文件所在目录 | `./gvcf_files` |
| `--output` | `-o` | 输出结果目录 | `./joint_output` |

### 可选参数

| 参数 | 简写 | 默认值 | 说明 |
|------|------|--------|------|
| `--threads` | `-t` | 88 | 每个任务的线程数 |
| `--tmp-dir` | `-T` | `./tmp` | 临时文件目录 |
| `--script` | `-s` | `run_gtx_joint.sh` | 输出脚本文件名 |
| `--faketime` | `-f` | `2020-10-20 00:00:00` | faketime时间字符串 |
| `--pattern` | `-p` | None | 染色体过滤正则表达式 |
| `--window` | `-w` | None | 区间大小(bp)，设置后按区间拆分 |

## 输出说明

### 生成的脚本格式

生成的脚本包含所有GTX joint calling命令，每条命令处理一个染色体或区间。

**按染色体模式示例：**
```bash
faketime '2020-10-20 00:00:00' /path/to/gtx joint \
    -r genome.fa \
    -o output/Chr01.joint.vcf.gz \
    -L Chr01 \
    --tmp-dir ./tmp \
    -t 88 \
    -v sample1.g.vcf.gz -v sample2.g.vcf.gz -v sample3.g.vcf.gz
```

**按区间模式示例：**
```bash
faketime '2020-10-20 00:00:00' /path/to/gtx joint \
    -r genome.fa \
    -o output/Chr01_1-10000000.joint.vcf.gz \
    -L Chr01:1-10000000 \
    --tmp-dir ./tmp \
    -t 88 \
    -v sample1.g.vcf.gz -v sample2.g.vcf.gz -v sample3.g.vcf.gz
```

### 输出文件

- **run_gtx_joint.sh** - 生成的命令脚本
- **gtx_joint.log** - 程序运行日志
- **.joint.vcf.gz** - 每个染色体/区间的联合VCF文件

## 执行方式

### 串行执行

```bash
bash run_gtx_joint.sh
```

适合小规模项目或单机环境。

### GNU Parallel并行执行

```bash
# 并行执行4个任务
parallel -j 4 < run_gtx_joint.sh

# 或者指定日志
parallel -j 4 --joblog parallel.log < run_gtx_joint.sh
```

适合多核CPU环境，可以显著提升速度。

### 使用集群批量提交

```bash
# 使用你的批量提交脚本
your_batch_submit_script.sh run_gtx_joint.sh
```

适合HPC集群环境。

## 区间模式的后续处理

如果使用`-w`参数按区间拆分，运行完成后需要合并同一染色体的所有区间VCF文件：

```bash
# 合并Chr01的所有区间
bcftools concat -o Chr01.merged.vcf.gz output/Chr01_*.joint.vcf.gz

# 或循环处理所有染色体
for chr in Chr01 Chr02 Chr03; do
    bcftools concat -o ${chr}.merged.vcf.gz output/${chr}_*.joint.vcf.gz
done
```

## 资源建议

### 硬件配置

- **内存**：建议≥32GB/任务，大基因组项目需≥64GB
- **磁盘**：预留输入数据的2-3倍空间
- **CPU**：根据内存情况调整并行任务数

### 性能优化

- 💾 使用SSD作为临时目录可以显著提升I/O性能
- 🧵 根据可用内存调整并行任务数和线程数
- 📊 监控系统资源使用情况，避免过度订阅

## 故障排除

### 1. GTX可执行文件不存在

```
错误: GTX可执行文件不存在: /path/to/gtx
```

**解决方案：**
- 使用`-g`参数指定正确的GTX路径
- 确保GTX文件具有可执行权限：`chmod +x gtx`

### 2. 未找到任何*.g.vcf.gz文件

```
错误: 未找到任何 *.g.vcf.gz 文件
```

**解决方案：**
- 检查`-i`指定的目录是否正确
- 确认GVCF文件以`.g.vcf.gz`结尾
- 检查目录权限

### 3. 参考基因组索引不存在

```
错误: 参考基因组索引不存在: genome.fa.fai
```

**解决方案：**
```bash
samtools faidx genome.fa
```

### 4. faketime未找到

```
警告: faketime 未找到，将不使用 faketime
```

这是警告信息，不影响核心功能。如果需要使用faketime：

```bash
# conda安装
conda install -c conda-forge faketime

# 或Ubuntu/Debian安装
sudo apt-get install faketime
```

## 应用场景

- 🧬 **群体遗传学研究** - 大规模样本的联合变异检测
- 📊 **GWAS分析预处理** - 为关联分析提供高质量的基因型数据
- 🌱 **育种项目** - 品种资源的基因型鉴定和比较
- 🔬 **进化研究** - 物种基因组变异分析

## 技术细节

### 染色体过滤

使用Python正则表达式进行匹配，常见的模式示例：

| 模式 | 说明 | 匹配示例 |
|------|------|----------|
| `^Chr[0-9]+$` | 只匹配主染色体 | Chr01, Chr02, ..., Chr99 |
| `^chr[0-9]+$` | 小写开头的染色体 | chr1, chr2, ..., chrX |
| `^[0-9]+$` | 纯数字染色体 | 1, 2, 3, ..., X, Y |
| `^(Chr|chr)[0-9]+$` | Chr或chr开头 | Chr01, chr01 |

### 内存估算

区间大小与内存使用的关系：

| 区间大小 | 预估内存 | 适用场景 |
|----------|----------|----------|
| 全染色体 | >64GB | 小基因组，充足内存 |
| 50M | ~32GB | 中等基因组 |
| 10M | ~16GB | 大基因组（推荐） |
| 5M | ~8GB | 内存受限环境 |

## 版本历史

- **v1.0.0** (2025-12-25)
  - 初始版本发布
  - 支持按染色体和区间拆分
  - 支持染色体过滤和faketime

## 相关工具

- `biopytools gtx` - GTX WGS单样品分析工具
- `biopytools fastq2vcf-gtx` - FastQ到VCF(GTX)全流程
- `biopytools gatk-joint` - GATK Joint Genotyping工具

## 参考文献

如果使用本工具，请引用GTX软件和相关方法学文献。
