# KMERIA 输入文件准备指南

## 快速生成示例文件

KMERIA提供了一个辅助脚本，帮助您快速创建输入文件模板。

### 方法1: 生成随机示例数据

用于测试流程：

```bash
python3 /share/org/YZWL/yzwl_lixg/software/biopytools/biopytools/kmeria/scripts/create_example_inputs.py \
    example \
    -n 10 \
    -o ./kmeria_input_files
```

**参数说明**：
- `-n 10`: 生成10个样本
- `-o ./kmeria_input_files`: 输出到指定目录

**生成的文件**：
- `samples.txt` - 样本列表
- `depth.txt` - 深度文件（随机值30-60X）
- `pheno.txt` - 表型文件（随机值-2到2）

---

### 方法2: 从FASTQ目录生成

如果您已有FASTQ文件，可以自动提取样本名：

```bash
python3 /share/org/YZWL/yzwl_lixg/software/biopytools/biopytools/kmeria/scripts/create_example_inputs.py \
    from-fastq \
    -i /path/to/your/fastq_files \
    -o ./kmeria_input_files
```

**生成的文件**：
- `samples.txt` - 从FASTQ文件名提取的样本列表
- `depth.txt` - 深度文件（默认30X，需要手动更新）
- `pheno.txt` - 表型文件（默认0.0，需要手动更新）

**⚠️ 重要**：使用此方法后，**必须手动编辑** `depth.txt` 和 `pheno.txt` 填入真实数据！

---

## 手动创建输入文件

如果您想手动创建文件，请遵循以下格式：

### 1. samples.txt

纯文本文件，每行一个样本名：

```text
sample1
sample2
sample3
```

**要求**：
- 样本名必须与FASTQ文件的前缀一致
- 例如：`sample1_R1.fq.gz` → 样本名为 `sample1`

---

### 2. depth.txt

Tab分隔文件（TSV）：

```text
sample1	45.2
sample2	52.8
sample3	38.9
```

**要求**：
- 第1列：样本名（与samples.txt一致）
- 第2列：测序深度（浮点数）
- 使用Tab分隔（`\t`）

**如何获取深度值**：

如果有BAM文件：
```bash
# 使用samtools计算平均深度
samtools depth -a sample.bam | awk '{sum+=$3} END {print "sample\t"sum/NR}'
```

或使用覆盖度统计：
```bash
# 计算覆盖度
samtools flagstat sample.bam
bedtools genomecov -ibam sample.bam
```

---

### 3. pheno.txt

Tab或空格分隔：

```text
sample1	1.5
sample2	2.3
sample3	1.8
```

**要求**：
- 第1列：样本名（与samples.txt一致）
- 第2列：表型值（数值型）
  - 连续型性状：如产量、株高（小数）
  - 分类性状：0/1或其他编码

**表型数据预处理建议**：

1. **标准化表型值**（推荐）：
```r
# R语言示例
phenotype <- scale(phenotype)  # mean=0, sd=1
```

2. **去除异常值**：
```r
# 保留3个标准差范围内的值
phenotype <- phenotype[abs(phenotype) < 3]
```

3. **分类性状编码**：
```text
# 二分类
resistant	0
susceptible	1

# 多分类
category_A	1
category_B	2
category_C	3
```

---

## 验证输入文件

运行分析前，建议验证文件格式：

### 检查样本名一致性

```bash
# 1. 提取所有文件的样本名
cut -f1 samples.txt > samples_names.txt
cut -f1 depth.txt > depth_names.txt
cut -f1 pheno.txt > pheno_names.txt

# 2. 检查是否一致
diff samples_names.txt depth_names.txt
diff samples_names.txt pheno_names.txt

# 3. 或者使用一行命令
awk 'NR==FNR{a[$1]=1; next} !($1 in a){print "Missing in depth:", $1}' samples.txt depth.txt
```

### 统计信息

```bash
# 样本数量
wc -l samples.txt depth.txt pheno.txt

# 深度范围
awk '{print $2}' depth.txt | sort -n | head -1  # 最小深度
awk '{print $2}' depth.txt | sort -n | tail -1  # 最大深度
awk '{sum+=$2; count++} END {print "Mean:", sum/count}' depth.txt  # 平均深度

# 表型值范围
awk '{print $2}' pheno.txt | sort -n | head -1  # 最小表型
awk '{print $2}' pheno.txt | sort -n | tail -1  # 最大表型
awk '{sum+=$2; count++} END {print "Mean:", sum/count}' pheno.txt  # 平均表型
```

---

## 完整工作流程示例

### 示例1: 使用真实数据

```bash
# 1. 准备FASTQ文件目录
mkdir -p my_project/fastq
# 将您的FASTQ文件放到此目录

# 2. 生成输入文件模板
python3 /share/org/YZWL/yzwl_lixg/software/biopytools/biopytools/kmeria/scripts/create_example_inputs.py \
    from-fastq \
    -i my_project/fastq \
    -o my_project/input_files

# 3. 编辑depth.txt和pheno.txt，填入真实数据
vi my_project/input_files/depth.txt
vi my_project/input_files/pheno.txt

# 4. 运行KMERIA分析
biopytools kmeria pipeline \
    -i my_project/fastq \
    --samples my_project/input_files/samples.txt \
    -d my_project/input_files/depth.txt \
    -p my_project/input_files/pheno.txt \
    -o my_project/kmeria_results \
    -t 24
```

### 示例2: 测试流程

```bash
# 生成10个样本的测试数据
python3 /share/org/YZWL/yzwl_lixg/software/biopytools/biopytools/kmeria/scripts/create_example_inputs.py \
    example \
    -n 10 \
    -o ./test_input

# 查看生成的文件
ls -lh test_input/
cat test_input/samples.txt
cat test_input/depth.txt
cat test_input/pheno.txt

# 使用测试数据运行流程
biopytools kmeria pipeline \
    -i /path/to/test/fastq \
    --samples test_input/samples.txt \
    -d test_input/depth.txt \
    -p test_input/pheno.txt \
    -o test_results \
    -t 8
```

---

## 常见问题

### Q1: 深度文件如何获取？

**A**: 如果您有BAM文件：
```bash
# 方法1: 使用samtools depth
samtools depth sample.bam | awk '{sum+=$3} END {print sum/NR}'

# 方法2: 使用mosdepth
mosdepth sample sample.bam
cat sample.mosdepth.global.dist.txt
```

如果没有BAM文件，可以根据测序数据量估算：
```
深度 = (总reads数 × read长度) / 基因组大小
```

### Q2: 表型值应该如何预处理？

**A**:
1. **连续型性状**：建议标准化（mean=0, sd=1）
2. **分类性状**：使用数值编码（0,1,2...）
3. **去除异常值**：建议删除超过3个标准差的值
4. **检查分布**：确保没有严重的偏态

### Q3: FASTQ文件命名有什么要求？

**A**:
- 支持双端：`sample_R1.fq.gz` + `sample_R2.fq.gz`
- 或：`sample_1.fq.gz` + `sample_2.fq.gz`
- 自动识别配对关系
- 支持压缩格式：`.gz`

### Q4: 可以有重复的样本名吗？

**A**: 不可以。样本名必须唯一且在所有文件中保持一致。

---

## 参考资源

- **完整文档**: `/share/org/YZWL/yzwl_lixg/software/biopytools/docs/kmeria.md`
- **官方示例**: `/share/org/YZWL/yzwl_lixg/tmp/test_kmeria/KMERIA-main/examples/`
- **KMERIA Wiki**: https://github.com/Sh1ne111/KMERIA/wiki
