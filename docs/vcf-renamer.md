# VCF样品名称重命名工具

## 功能概述

`biopytools vcf-renamer` 是一个专门用于重命名VCF文件中样品名称的工具。它可以将长样品名称重命名为简洁的格式（如 S1, S2, S3...），防止样品名称被GATK、GTX等软件截断的问题。

## 核心功能

- 📋 **自动提取样品名称** - 从VCF文件中自动读取所有样品名称
- 🏷️ **批量重命名** - 按顺序重命名样品（S1, S2, S3...或自定义前缀）
- 📝 **生成映射文件** - 保存新旧样品名称对应关系
- 🔄 **bcftools集成** - 使用bcftools reheader进行安全的重命名操作
- 📇 **自动索引** - 重命名后自动创建VCF索引文件
- 🧹 **可选清理** - 可选择保留或删除映射文件

## 安装依赖

使用本工具需要预先安装bcftools：

```bash
# 使用conda安装
conda install -c bioconda bcftools

# 或使用apt安装 (Ubuntu/Debian)
sudo apt update && sudo apt install bcftools
```

## 使用方法

### 基本用法（使用默认前缀S）

```bash
biopytools vcf-renamer \
    -i variation.filtered.snp.vcf.gz \
    -o variation.renamed.vcf.gz
```

这会将样品重命名为 S1, S2, S3...

### 自定义样品名前缀

```bash
biopytools vcf-renamer \
    -i input.vcf.gz \
    -o output.vcf.gz \
    -p Sample
```

样品将被重命名为 Sample1, Sample2, Sample3...

### 指定映射文件路径

```bash
biopytools vcf-renamer \
    -i input.vcf.gz \
    -o output.vcf.gz \
    -m sample_mapping.txt
```

### 不保留映射文件

```bash
biopytools vcf-renamer \
    -i input.vcf.gz \
    -o output.vcf.gz \
    --no-mapping
```

## 命令行参数

### 必需参数

| 参数 | 简写 | 说明 | 示例 |
|------|------|------|------|
| `--input` | `-i` | 输入VCF文件路径 | `input.vcf.gz` |
| `--output` | `-o` | 输出VCF文件路径 | `output.vcf.gz` |

### 可选参数

| 参数 | 简写 | 默认值 | 说明 |
|------|------|--------|------|
| `--prefix` | `-p` | `S` | 新样品名前缀 |
| `--mapping` | `-m` | 自动生成 | 映射文件路径 |
| `--no-mapping` | | False | 不保留映射文件 |

## 输出说明

### 1. 重命名后的VCF文件

样品名称已更新为新的格式。可以使用bcftools验证：

```bash
bcftools query -l output.vcf.gz
```

输出示例：
```
S1
S2
S3
```

### 2. VCF索引文件

自动生成`.tbi`索引文件：`output.vcf.gz.tbi`

### 3. 映射文件

保存新旧样品名称对应关系（除非使用`--no-mapping`）

格式示例：
```
Original_Sample_Name_1    S1
Original_Sample_Name_2    S2
Original_Sample_Name_3    S3
```

## 工作流程

1. **步骤1**: 提取原始样品名称
2. **步骤2**: 生成新旧ID映射表
3. **步骤3**: 使用bcftools reheader重命名
4. **步骤4**: 创建VCF索引

输出示例：
```
[2025-12-25 10:30:15] INFO: ============================================================
[2025-12-25 10:30:15] INFO: 🧬 开始VCF样品名称重命名 | Starting VCF Sample Name Renaming
[2025-12-25 10:30:15] INFO: ============================================================
[2025-12-25 10:30:15] INFO: 📋 步骤 1/4: 提取样品名称 | Step 1/4: Extracting sample names
[2025-12-25 10:30:15] INFO: ✅ 检测到 10 个样品 | Found 10 samples
[2025-12-25 10:30:15] INFO: 📝 步骤 2/4: 生成新旧ID映射表 | Step 2/4: Generating mapping
[2025-12-25 10:30:16] INFO: 🔄 步骤 3/4: 使用bcftools重命名样品 | Step 3/4: Renaming samples
[2025-12-25 10:30:45] INFO: 📇 步骤 4/4: 创建VCF索引 | Step 4/4: Creating index
[2025-12-25 10:30:46] INFO: 🎉 重命名完成！| Renaming completed!
```

## 应用场景

### 1. 防止样品名截断

某些软件（如GATK、GTX）对样品名称长度有限制，长名称会被截断：

```bash
# 原始样品名称（可能被截断）
VeryLongSampleNameFromPopulationA_2024
VeryLongSampleNameFromPopulationB_2024

# 重命名后（不会被截断）
S1
S2
```

### 2. 简化样品名称

简化复杂或难以处理的样品名称：

```bash
# 原始名称
Sample@2024#01_WGS_R1
Sample@2024#01_WGS_R2

# 重命名后
S1
S2
```

### 3. 批量标准化

为多个项目统一样品命名格式：

```bash
# 前缀"Pop"表示群体研究
biopytools vcf-renamer -i input.vcf.gz -o output.vcf.gz -p Pop

# 结果: Pop1, Pop2, Pop3...
```

## 故障排除

### 1. 未找到bcftools命令

```
错误: 未找到 'bcftools' 命令
```

**解决方案：**
```bash
conda install -c bioconda bcftools
```

### 2. 输入VCF文件不存在

```
错误: 输入VCF文件不存在: input.vcf.gz
```

**解决方案：**
- 检查`-i`参数指定的文件路径
- 确认文件格式为`.vcf.gz`或`.vcf`
- 使用`ls -lh`检查文件是否存在

### 3. 提取样品名称失败

```
错误: 提取样品名称失败
```

**解决方案：**
```bash
# 检查VCF文件格式
bcftools view -h input.vcf.gz | head -20

# 验证VCF文件完整性
bcftools index -i input.vcf.gz
```

### 4. 重命名失败

```
错误: 重命名失败
```

**可能原因：**
- 输出目录无写权限
- 磁盘空间不足
- 映射文件格式错误

**解决方案：**
```bash
# 检查输出目录权限
ls -ld output_dir/

# 检查磁盘空间
df -h

# 手动测试bcftools
bcftools reheader -s mapping.txt -o test.vcf.gz input.vcf.gz
```

## 最佳实践

### 1. 保留映射文件

建议保留映射文件，方便后续追溯原始样品名：

```bash
# 默认保留映射文件
biopytools vcf-renamer -i input.vcf.gz -o output.vcf.gz
```

映射文件可以用于：
- 恢复原始样品名称
- 记录样品对应关系
- 文档化和结果追溯

### 2. 选择有意义的前缀

根据项目选择有意义的样品名前缀：

| 场景 | 前缀建议 | 结果示例 |
|------|----------|----------|
| 群体研究 | `Pop` | Pop1, Pop2, Pop3 |
| 品系研究 | `Strain` | Strain1, Strain2 |
| 处理组 | `Treat` | Treat1, Treat2 |
| 对照组 | `Ctrl` | Ctrl1, Ctrl2 |
| 时间序列 | `Time` | Time1, Time2 |

### 3. 验证结果

重命名后务必验证样品名称：

```bash
# 查看重命名后的样品名称
bcftools query -l output.vcf.gz

# 检查VCF文件完整性
bcftools stats output.vcf.gz | less
```

### 4. 批量处理

结合shell脚本批量处理多个VCF文件：

```bash
#!/bin/bash
# 批量重命名VCF文件

for vcf in *.vcf.gz; do
    output="renamed_${vcf}"
    biopytools vcf-renamer -i "$vcf" -o "$output" -p Sample
done
```

## 恢复原始样品名称

如果需要恢复原始样品名称，可以使用映射文件：

```bash
# 1. 读取映射文件
cat vcf_name_mapping.txt
# 输出: OriginalName    S1

# 2. 创建反向映射文件
cat > reverse_mapping.txt << EOF
S1    OriginalName
S2    OriginalName2
...

# 3. 使用bcftools恢复
bcftools reheader -s reverse_mapping.txt -o restored.vcf.gz renamed.vcf.gz
```

## 技术细节

### bcftools reheader工作原理

`bcftools reheader`命令的工作流程：
1. 读取VCF文件的头部信息
2. 根据映射文件替换样品名称
3. 保持VCF文件的其他内容不变
4. 输出新的VCF文件

### 映射文件格式要求

映射文件必须是两列格式：
- 第一列：原始样品名称
- 第二列：新的样品名称
- 分隔符：空格或制表符

示例：
```
Sample_A_2024    S1
Sample_B_2024    S2
Sample_C_2024    S3
```

## 版本历史

- **v1.0.0** (2025-12-25)
  - 初始版本发布
  - 支持样品名称重命名
  - 支持自定义前缀
  - 支持映射文件生成
  - 自动索引创建

## 相关工具

- `biopytools gtx-joint` - GTX Joint Calling命令生成工具
- `biopytools vcf-filter` - VCF文件过滤工具
- `biopytools annovar` - ANNOVAR变异注释工具
