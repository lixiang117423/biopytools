# 🧬 GTX Joint Calling 命令生成模块

**按染色体/区间生成GTX joint calling分析脚本的专用工具 | Specialized Tool for Generating GTX Joint Calling Analysis Scripts by Chromosome/Window**

## 📖 功能概述 | Overview

GTX Joint Calling 命令生成模块是一个专门用于生成GTX joint calling分析脚本的工具。通过按染色体或固定区间拆分任务，支持大规模样本的并行joint calling分析，有效避免内存溢出问题，提高分析效率。

## ✨ 主要特性 | Key Features

- **🔄 自动化脚本生成**: 扫描GVCF文件并自动生成执行脚本
- **📊 灵活拆分模式**: 支持按染色体或按固定区间拆分
- **🎯 染色体过滤**: 支持正则表达式过滤特定染色体
- **🛡️ 环境检查**: 完整的依赖检查和验证机制
- **📝 详细统计**: 完整的统计信息和执行建议
- **🚀 并行优化**: 生成的脚本支持GNU Parallel并行执行
- **🧠 内存优化**: 区间拆分模式避免大染色体导致的内存溢出
- **📋 自动备份**: 自动备份已存在的脚本文件

## 🚀 快速开始 | Quick Start

### 基本用法 | Basic Usage

```bash
# 按染色体拆分（最常用）
biopytools gtx-joint \
    -g /path/to/gtx \
    -r genome.fa \
    -i ./gvcf_dir \
    -o ./output

# 按10M区间拆分（防止内存溢出）
biopytools gtx-joint \
    -g /path/to/gtx \
    -r genome.fa \
    -i ./gvcf_dir \
    -o ./output \
    -w 10000000
```

### 高级用法 | Advanced Usage

```bash
# 自定义线程和临时目录
biopytools gtx-joint \
    -g /path/to/gtx \
    -r genome.fa \
    -i ./gvcf_dir \
    -o ./output \
    -t 24 \
    -T /tmp

# 只处理主染色体（使用正则过滤）
biopytools gtx-joint \
    -g /path/to/gtx \
    -r genome.fa \
    -i ./gvcf_dir \
    -o ./output \
    -p "^Chr[0-9]{2}$"
```

## 📋 参数说明 | Parameters

### 必需参数 | Required Parameters

| 参数 | 描述 | 示例 |
|------|------|------|
| `-g, --gtx` | GTX可执行文件路径 | `-g /path/to/gtx` |
| `-r, --ref` | 参考基因组文件(.fa) | `-r genome.fa` |
| `-i, --input` | GVCF文件所在目录 | `-i ./gvcf_dir` |
| `-o, --output` | 输出结果目录 | `-o ./output` |

### 可选参数 | Optional Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-t, --threads` | `88` | 每个任务的线程数 |
| `-T, --tmp-dir` | `./tmp` | 临时文件目录 |
| `-s, --script` | `run_gtx_joint.sh` | 输出脚本文件名 |
| `-f, --faketime` | `2020-10-20 00:00:00` | faketime时间设置 |
| `-p, --pattern` | `""` | 染色体过滤正则表达式 |
| `-w, --window` | `None` | 区间大小(bp)，设置后按区间拆分 |

### 日志参数 | Logging Parameters

| 参数 | 描述 |
|------|------|
| `-v, --verbose` | 详细输出模式 (-v: INFO, -vv: DEBUG) |
| `--quiet` | 静默模式 (只输出ERROR) |
| `--log-file` | 日志文件路径 |

## 💡 使用示例 | Usage Examples

### 示例1：基本染色体拆分 | Example 1: Basic Chromosome Splitting

```bash
# 最常用的场景：按染色体拆分
biopytools gtx-joint \
    -g /share/org/YZWL/yzwl_lixg/software/gtx/bin/gtx \
    -r Arabidopsis_thaliana.fa \
    -i ./gvcf_files \
    -o ./joint_results
```

**生成的文件**:
- `./joint_results/run_gtx_joint.sh` - 执行脚本
- `./joint_results/Chr1.joint.vcf.gz` - 染色体1的结果
- `./joint_results/Chr2.joint.vcf.gz` - 染色体2的结果
- ...

### 示例2：按10M区间拆分（推荐用于大基因组）| Example 2: Split by 10M Windows (Recommended for Large Genomes)

```bash
# 大基因组或大样本量时，按区间拆分避免内存溢出
biopytools gtx-joint \
    -g /path/to/gtx \
    -r wheatgenome.fa \
    -i ./gvcf_dir \
    -o ./joint_results \
    -w 10000000
```

**⚠️ 重要提示 - 区间拆分后的合并步骤**:

使用区间拆分（`-w`参数）后，每个染色体会被分割成多个区间文件。**分析完成后必须合并同一染色体的所有区间文件**：

```bash
# 进入输出目录
cd ./joint_results

# 合并Chr1的所有区间
bcftools concat -o Chr1.merged.vcf.gz Chr1_*.joint.vcf.gz
bcftools index -t Chr1.merged.vcf.gz

# 合并Chr2的所有区间
bcftools concat -o Chr2.merged.vcf.gz Chr2_*.joint.vcf.gz
bcftools index -t Chr2.merged.vcf.gz

# 对所有染色体执行相同操作
for chr in Chr[0-9]*; do
    echo "Merging $chr..."
    bcftools concat -o ${chr}.merged.vcf.gz ${chr}_*.joint.vcf.gz
    bcftools index -t ${chr}.merged.vcf.gz
done

# 最后合并所有染色体（如果需要全基因组单文件）
bcftools concat -o genome_merged.vcf.gz Chr*.merged.vcf.gz
bcftools index -t genome_merged.vcf.gz
```

**批量合并脚本示例**:
```bash
#!/bin/bash
# merge_intervals.sh - 区间结果合并脚本

OUTPUT_DIR="./joint_results"
cd "$OUTPUT_DIR"

# 按染色体合并区间文件
for chr_prefix in $(ls | grep -oP 'Chr\w+(?=_\d+-\d+\.joint\.vcf\.gz)' | sort -u); do
    echo "处理染色体 | Processing chromosome: $chr_prefix"

    # 查找该染色体的所有区间文件
    files=(${chr_prefix}_*.joint.vcf.gz)

    if [ ${#files[@]} -gt 1 ]; then
        # 合并
        bcftools concat -Oz -o ${chr_prefix}.merged.vcf.gz ${chr_prefix}_*.joint.vcf.gz
        bcftools index -t ${chr_prefix}.merged.vcf.gz
        echo "  ✓ 已生成 | Created: ${chr_prefix}.merged.vcf.gz"
    else
        # 只有一个文件，直接重命名
        mv ${files[0]} ${chr_prefix}.merged.vcf.gz
        echo "  ✓ 已重命名 | Renamed: ${chr_prefix}.merged.vcf.gz"
    fi
done

echo ""
echo "✅ 所有染色体合并完成！| All chromosomes merged!"
echo "📊 合并后的文件 | Merged files:"
ls -lh *.merged.vcf.gz*
```

### 示例3：只处理主染色体 | Example 3: Process Main Chromosomes Only

```bash
# 只处理Chr01-Chr09，跳过scaffold和contig
biopytools gtx-joint \
    -g /path/to/gtx \
    -r genome.fa \
    -i ./gvcf_dir \
    -o ./output \
    -p "^Chr[0-9]{2}$"
```

**常用正则表达式**:
- `^Chr[0-9]+$` - 只处理Chr1, Chr2...（单数染色体）
- `^Chr[0-9]{2}$` - 只处理Chr01-Chr99（两位数染色体）
- `^(chr)?[0-9]+$` - 匹配Chr1和chr1格式

### 示例4：完整参数配置 | Example 4: Full Parameter Configuration

```bash
biopytools gtx-joint \
    -g /opt/gtx/bin/gtx \
    -r /data/reference/human_g1k_v37.fa \
    -i /data/gvcf/1000genomes \
    -o /data/joint_output \
    -t 32 \
    -T /scratch/tmp \
    -s human_g1k_joint.sh \
    -w 50000000 \
    -v
```

### 示例5：批量提交到集群 | Example 5: Batch Submit to Cluster

```bash
# 1. 生成执行脚本
biopytools gtx-joint \
    -g /path/to/gtx \
    -r genome.fa \
    -i ./gvcf_dir \
    -o ./output \
    -w 10000000

# 2. 使用GNU Parallel并行执行（推荐）
cd ./output
parallel -j 10 < run_gtx_joint.sh

# 3. 或使用自定义的批量提交脚本
your_cluster_submit_script.sh ./output/run_gtx_joint.sh

# 4. 如果使用了区间拆分，执行合并脚本
bash merge_intervals.sh
```

## 📁 输入文件格式 | Input File Formats

### GVCF文件 | GVCF Files

标准的GATK GVCF格式文件（`.g.vcf.gz`）：

```bash
gvcf_dir/
├── sample1.g.vcf.gz
├── sample1.g.vcf.gz.tbi
├── sample2.g.vcf.gz
├── sample2.g.vcf.gz.tbi
└── ...
```

**要求**:
- 必须是 `.g.vcf.gz` 格式
- 必须有对应的 `.tbi` 索引文件
- 所有样本应使用相同的参考基因组

### 参考基因组 | Reference Genome

标准FASTA格式，必须已建立索引：

```bash
genome.fa
genome.fa.fai  # 必需
```

如果缺少索引，需要先创建：
```bash
samtools faidx genome.fa
```

## 📊 输出结果 | Output Results

### 染色体拆分模式 | Chromosome Splitting Mode

```
output_dir/
├── run_gtx_joint.sh          # 执行脚本
├── Chr1.joint.vcf.gz         # Chr1的joint calling结果
├── Chr1.joint.vcf.gz.tbi     # 索引
├── Chr2.joint.vcf.gz
├── Chr2.joint.vcf.gz.tbi
└── ...
```

### 区间拆分模式 | Window Splitting Mode

```
output_dir/
├── run_gtx_joint.sh              # 执行脚本
├── Chr1_1-10000000.joint.vcf.gz  # Chr1的区间1
├── Chr1_10000001-20000000.joint.vcf.gz  # Chr1的区间2
├── Chr1_20000001-27642331.joint.vcf.gz  # Chr1的区间3
├── Chr2_1-10000000.joint.vcf.gz
├── Chr2_10000001-18053412.joint.vcf.gz
└── ...
```

**⚠️ 区间模式必须合并同一染色体的文件** (见示例2)

## 🔧 系统要求 | System Requirements

### 依赖软件 | Dependencies

- **GTX** (版本 1.0+)
  - GTX joint calling工具
- **bcftools** (版本 1.10+)
  - 用于合并区间拆分的结果
- **samtools** (版本 1.10+)
  - 用于参考基因组索引
- **Python** (版本 3.7+)
- **Python包**:
  - `click` - 命令行界面

### 硬件建议 | Hardware Recommendations

| 参数 | 小规模 | 中等规模 | 大规模 |
|------|--------|----------|--------|
| **样本数** | < 100 | 100-1000 | > 1000 |
| **内存/任务** | 16GB | 32GB | 64GB |
| **磁盘** | 100GB | 500GB | 1TB+ |
| **线程/任务** | 8-16 | 24-48 | 88+ |

## ⚙️ 运行模式对比 | Running Mode Comparison

| 模式 | 优点 | 缺点 | 适用场景 |
|------|------|------|----------|
| **按染色体拆分**<br>(默认) | • 简单直接<br>• 结果无需合并<br>• 易于管理 | • 大染色体可能内存溢出 | • 小基因组<br>• 样本数< 500<br>• 染色体< 20条 |
| **按区间拆分**<br>(-w参数) | • 避免内存溢出<br>• 适合大基因组<br>• 任务更均衡 | • 需要合并结果<br>• 任务数更多 | • 大基因组<br>• 样本数> 500<br>• 染色体> 20条 |

## 🚀 执行流程 | Execution Workflow

### 完整工作流程（区间拆分模式）| Complete Workflow (Window Splitting Mode)

```bash
# 步骤1: 生成命令脚本
biopytools gtx-joint \
    -g /path/to/gtx \
    -r genome.fa \
    -i ./gvcf_dir \
    -o ./joint_output \
    -w 10000000

# 步骤2: 执行joint calling
cd ./joint_output
parallel -j 10 < run_gtx_joint.sh

# 步骤3: 合并区间结果（必需！）
# 方法1: 手动合并每个染色体
bcftools concat -Oz -o Chr1.merged.vcf.gz Chr1_*.joint.vcf.gz
bcftools index -t Chr1.merged.vcf.gz

# 方法2: 使用批量合并脚本
bash merge_intervals.sh  # 见示例2

# 步骤4: 合并所有染色体（可选）
bcftools concat -Oz -o genome_full.vcf.gz Chr*.merged.vcf.gz
bcftools index -t genome_full.vcf.gz
```

### 完整工作流程（染色体拆分模式）| Complete Workflow (Chromosome Splitting Mode)

```bash
# 步骤1: 生成命令脚本
biopytools gtx-joint \
    -g /path/to/gtx \
    -r genome.fa \
    -i ./gvcf_dir \
    -o ./joint_output

# 步骤2: 执行joint calling
cd ./joint_output
parallel -j 10 < run_gtx_joint.sh

# 步骤3: 合并所有染色体（可选）
bcftools concat -Oz -o genome_full.vcf.gz Chr*.joint.vcf.gz
bcftools index -t genome_full.vcf.gz
```

## ⚠️ 注意事项 | Important Notes

1. **区间拆分必须合并**: 使用`-w`参数后，必须合并同一染色体的区间文件
2. **索引文件**: GVCF文件必须有`.tbi`索引，否则会警告
3. **内存管理**: 大样本量建议使用区间拆分模式
4. **磁盘空间**: 确保有足够的磁盘空间（输入数据的2-3倍）
5. **参考基因组**: 必须有`.fai`索引文件
6. **线程设置**: 建议根据服务器配置合理设置线程数
7. **时间设置**: faketime用于解决GTX的许可证时间限制问题

## 🐛 故障排除 | Troubleshooting

### 常见问题 | Common Issues

**Q: 提示缺少`.fai`索引文件**
```bash
# 创建参考基因组索引
samtools faidx genome.fa
```

**Q: 内存不足错误**
```bash
# 使用区间拆分模式
biopytools gtx-joint ... -w 10000000

# 或减少每个任务的线程数
biopytools gtx-joint ... -t 16
```

**Q: 如何确定区间大小？**
```bash
# 根据染色体长度和可用内存决定
# 小基因组(< 500M): -w 10000000  (10M)
# 中基因组(500M-2G): -w 50000000  (50M)
# 大基因组(> 2G): -w 100000000  (100M)
```

**Q: 合并时出现顺序错误**
```bash
# 使用bcftools concat的默认顺序（按位置）
bcftools concat -o merged.vcf.gz Chr_*.joint.vcf.gz

# 或明确指定顺序
ls -v Chr_*.joint.vcf.gz | xargs bcftools concat -o merged.vcf.gz
```

## 📚 应用场景 | Application Scenarios

1. **大规模群体测序**: 1000+样本的joint calling
2. **植物泛基因组**: 大基因组的变异检测
3. **临床外显子组**: 大批样本的exome seq数据分析
4. **种质资源评估**: 自然群体遗传变异分析
5. **进化基因组学**: 物种间比较分析

## 📄 许可证 | License

本项目采用MIT许可证 - 详见 [LICENSE](LICENSE) 文件

---

## 🔬 引用信息 | Citation

如果在学术研究中使用此工具，请引用相关文献：

```
[您的项目名称] (2025).
GTX Joint Calling Command Generator: Automated script generation for large-scale joint variant calling.
https://github.com/yourusername/biopytools
```
