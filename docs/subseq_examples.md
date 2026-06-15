# 🧬 Subseq模块使用示例

本文档提供了subseq模块的实际使用示例和常见场景应用。

## 📋 基础示例

### 示例1：从基因组中提取特定染色体

```bash
# 创建染色体ID列表
echo -e "chr1\nchr2\nchr3\nchr4\nchr5" > chromosomes.txt

# 提取指定染色体
biopytools subseq -i human_genome.fasta -l chromosomes.txt -o selected_chromosomes.fasta
```

### 示例2：提取特定基因家族成员

```bash
# 提取所有以"ABC_"开头的基因
biopytools subseq -i annotated_genome.fasta -p "ABC_" -o ABC_gene_family.fasta --pattern-type startswith

# 提取所有包含"kinase"的序列
biopytools subseq -i protein_sequences.fasta -p "kinase" -o kinase_sequences.fasta
```

### 示例3：基于长度的序列筛选

```bash
# 提取长度大于1000bp的序列
biopytools subseq -i all_sequences.fasta -o long_sequences.fasta --length-only --min-length 1000

# 提取长度在500-2000bp之间的序列
biopytools subseq -i transcripts.fasta -o medium_transcripts.fasta --length-only --min-length 500 --max-length 2000
```

## 🔬 生物信息学应用场景

### 场景1：转录组分析

```bash
# 从转录组中提取完整转录本（>500bp）
biopytools subseq -i raw_transcripts.fasta -o complete_transcripts.fasta --length-only --min-length 500

# 提取特定基因的所有转录本
echo -e "gene1\n gene2\n gene3" > target_genes.txt
biopytools subseq -i transcriptome.fasta -l target_genes.txt -o target_transcripts.fasta
```

### 场景2：基因组组装质量控制

```bash
# 移除过短的contigs（<500bp）
biopytools subseq -i assembly.fasta -o filtered_assembly.fasta --length-only --min-length 500

# 统计序列长度分布
grep "^>" filtered_assembly.fasta | wc -l  # 统计过滤后序列数
```

### 场景3：标记基因提取

```bash
# 使用正则表达式提取COX基因
biopytools subseq -i mitochondrial_genome.fasta -p "COX[1-3]" -o COX_genes.fasta --pattern-type regex

# 提取rRNA基因
biopytools subseq -i bacterial_genome.fasta -p "16S\|23S\|5S" -o rRNA_genes.fasta --pattern-type regex
```

### 场景4：蛋白质序列处理

```bash
# 提取特定结构域的蛋白质
biopytools subseq -i proteome.fasta -p "kinase_domain" -o kinases.fasta --ignore-case

# 提取特定长度的蛋白质（50-500氨基酸）
biopytools subseq -i proteins.fasta -o medium_proteins.fasta --length-only --min-length 50 --max-length 500
```

## 📊 批量处理脚本

### 批量提取脚本示例

```bash
#!/bin/bash
# batch_extract.sh - 批量序列提取脚本

INPUT_DIR="input_fastas"
OUTPUT_DIR="output_fastas"
ID_LISTS_DIR="id_lists"

# 创建输出目录
mkdir -p "$OUTPUT_DIR"

# 处理每个FASTA文件
for fasta_file in "$INPUT_DIR"/*.fasta; do
    base_name=$(basename "$fasta_file" .fasta)
    id_file="$ID_LISTS_DIR/${base_name}_ids.txt"

    if [[ -f "$id_file" ]]; then
        echo "处理: $base_name"
        biopytools subseq \
            -i "$fasta_file" \
            -l "$id_file" \
            -o "$OUTPUT_DIR/${base_name}_extracted.fasta"

        echo "完成: $base_name -> ${base_name}_extracted.fasta"
    else
        echo "警告: 未找到ID文件 $id_file"
    fi
done
```

### 多条件筛选脚本

```bash
#!/bin/bash
# multi_filter.sh - 多条件序列筛选

input_file="$1"
output_dir="$2"

mkdir -p "$output_dir"

echo "开始多条件筛选..."

# 条件1: 提取基因序列
biopytools subseq \
    -i "$input_file" \
    -p "gene" \
    -o "$output_dir/gene_sequences.fasta" \
    --ignore-case

# 条件2: 提取长度>1000的序列
biopytools subseq \
    -i "$input_file" \
    -o "$output_dir/long_sequences.fasta" \
    --length-only \
    --min-length 1000

# 条件3: 提取特定模式的序列
biopytools subseq \
    -i "$input_file" \
    -p "chr[1-9]" \
    -o "$output_dir/chromosome_sequences.fasta" \
    --pattern-type regex

echo "筛选完成！结果保存在 $output_dir"
```

## 🧪 测试数据准备

### 创建测试FASTA文件

```bash
# 创建测试FASTA文件
cat > test_sequences.fasta << 'EOF'
>chr1_gene1
ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>chr2_gene2
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT
>chr1_gene3
TTTTAAAAAAAAAAAAAAAAGGGGGCCCCCCCCCCCCCCCCCCCCCCCCCC
>protein_kinase
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
>contig_00001
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
>contig_10001
ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>short_seq
ATGC
>gene_XYZ
GATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
EOF
```

### 创建测试ID列表

```bash
# 创建测试ID列表
cat > test_ids.txt << 'EOF'
chr1_gene1
chr1_gene3
gene_XYZ
nonexistent_id
EOF
```

## 🔍 结果验证

### 检查提取结果

```bash
# 检查输出文件
echo "提取的序列数量:"
grep "^>" output.fasta | wc -l

echo "提取的序列ID:"
grep "^>" output.fasta

# 检查序列长度
awk '/^>/ {if (seq) print length(seq); seq=""; print; next} {seq=seq $0} END {print length(seq)}' output.fasta
```

### 统计信息验证

```bash
# 验证日志中的统计信息
tail -20 subseq_extraction.log | grep "统计信息"
```

## ⚠️ 常见问题和解决方案

### 问题1：ID不匹配

**现象**: 提取序列数量为0或少于预期

**解决方案**:
```bash
# 检查ID是否存在于FASTA文件中
for id in $(cat ids.txt); do
    grep -q "^>$id" sequences.fasta && echo "找到: $id" || echo "未找到: $id"
done

# 使用模糊匹配
biopytools subseq -i sequences.fasta -p "partial_id" -o output.fasta --ignore-case
```

### 问题2：内存不足

**现象**: 处理大文件时内存不足

**解决方案**:
```bash
# 检查内存使用
free -h

# 分批处理
split -l 10000 large_ids.txt small_ids_
for file in small_ids_*; do
    biopytools subseq -i large.fasta -l "$file" -o "output_$file.fasta"
    cat "output_$file.fasta" >> final_output.fasta
done
```

### 问题3：正则表达式错误

**现象**: 模式匹配无结果

**解决方案**:
```bash
# 测试正则表达式
echo "gene_001_test" | grep -E "^gene_[0-9]{3}.*$"

# 使用简单模式验证
biopytools subseq -i test.fasta -p "gene" -o test_output.fasta
```

## 📈 性能优化建议

### 提高处理速度

1. **使用简单匹配模式**: 尽量避免复杂的正则表达式
2. **批量处理**: 将多个小文件合并后一次性处理
3. **过滤预处理**: 先使用长度筛选减少序列数量
4. **SSD存储**: 使用SSD存储提高I/O性能

### 内存优化

1. **分批处理**: 大ID列表分批处理
2. **清理中间文件**: 及时删除临时文件
3. **监控内存**: 使用`htop`监控内存使用

这些示例涵盖了subseq模块的主要使用场景和最佳实践。根据具体需求选择合适的提取方式和参数配置。