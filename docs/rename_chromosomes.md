# Rename Chromosomes - 染色体重命名工具

## 功能简介

将FASTA文件的序列进行标准化重命名，常用于基因组组装后的染色体编号规范：

- **前N条序列**：命名为 `Chr01`, `Chr02`, ..., `ChrNN`
- **剩余序列**：命名为 `HiC_scaffold_01`, `HiC_scaffold_02`, ...

## 系统要求

- Python >= 3.8
- AWK工具（系统自带）

## 安装方法

```bash
# 从源码安装
git clone https://github.com/lixiang117423/biopytools.git
cd biopytools
pip install -e .
```

## 使用方法

### 查看帮助

```bash
biopytools rename-chromosomes -h
```

### 基本用法

```bash
# 重命名基因组文件，前20条为染色体
biopytools rename-chromosomes \
    -i scaffolds_final.fa \
    -o renamed_genome.fa \
    -n 20
```

### 参数说明

| 短参数 | 长参数 | 类型 | 必需 | 说明 |
|--------|--------|------|------|------|
| `-i` | `--input` | str | 是 | 输入FASTA文件路径 |
| `-o` | `--output` | str | 是 | 输出FASTA文件路径 |
| `-n` | `--number` | int | 是 | 染色体数量 |

### 输出说明

工具会生成以下文件：

1. **renamed_genome.fa** - 重命名后的FASTA文件
2. **rename_chromosomes_YYYYMMDD_HHMMSS.log** - 运行日志

## 命名规则

### 染色体命名
- 格式：`ChrNN`（两位数，不足补零）
- 示例：`Chr01`, `Chr02`, `Chr03`, ..., `Chr20`

### Scaffold命名
- 格式：`HiC_scaffold_NN`（两位数，不足补零）
- 示例：`HiC_scaffold_01`, `HiC_scaffold_02`, ...

### 示例输出

```fasta
>Chr01
ATCGATCGATCGATCG...
>Chr02
GCTAGCTAGCTAGCTA...
...
>Chr20
TTAATTAATTAATTAA...
>HiC_scaffold_01
CCGGCCGGCCGGCCGG...
>HiC_scaffold_02
AATTAATTAATTAATT...
```

## 使用场景

### 1. 基因组组装后处理

```bash
# Hi-C scaffolding后的染色体重命名
biopytools rename-chromosomes \
    -i assembly_scaffolds.fa \
    -o final_chromosomes.fa \
    -n 24
```

### 2. 批量处理

```bash
# 处理多个样本
for sample in sample1 sample2 sample3; do
    biopytools rename-chromosomes \
        -i ${sample}_scaffolds.fa \
        -o ${sample}_renamed.fa \
        -n 20
done
```

## 技术细节

### 实现原理

使用AWK进行流式处理，避免大文件内存问题：

1. 读取FASTA文件的header行（以`>`开头）
2. 维护两个计数器：`chr_count`和`scaf_count`
3. 根据计数器位置决定命名格式
4. 序列内容保持不变

### 性能特点

- **内存效率**：流式处理，内存占用低
- **速度**：AWK原生处理，速度快
- **可扩展性**：支持超大FASTA文件

## 注意事项

1. **输入格式**：必须是标准的FASTA格式
2. **文件覆盖**：输出文件存在时会直接覆盖
3. **编号格式**：使用两位数格式（01-99），超过99条需调整代码
4. **日志管理**：每次运行生成新的日志文件

## 示例工作流

### 完整的基因组组装流程

```bash
# 1. 基因组组装
# hifiasm -o assembly -t 32 reads.fastq.gz

# 2. Hi-C scaffolding
# biopytools haphic -a assembly.fa -b hic.bam -c 24

# 3. 染色体重命名
biopytools rename-chromosomes \
    -i hap hic_scaffolds.fa \
    -o final_genome.fa \
    -n 24

# 4. 质量评估
# busco -i final_genome.fa -l embryophyta -m genome
```

## 常见问题

### Q: 如何处理超过99条染色体？
A: 需要修改代码中的格式化字符串，将`%02d`改为`%03d`（三位数）。

### Q: 可以自定义命名格式吗？
A: 可以修改`utils.py`中的AWK脚本部分。

### Q: 支持压缩的FASTA文件吗？
A: 当前版本不支持，需要先解压。

## 许可证

MIT License

## 作者信息

**李详 (Xiang Li)**
- Email: lixiang117423@gmail.com
- GitHub: [@lixiang117423](https://github.com/lixiang117423)

## 相关工具

- [HapHiC](./haphic.md) - Hi-C基因组scaffolding
- [BUSCO](./busco.md) - 基因组质量评估
- [EGAPx](./egapx-batch.md) - 基因预测
