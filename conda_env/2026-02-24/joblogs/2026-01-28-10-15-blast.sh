#!/bin/bash

# 设置变量
INPUT_FASTA="41.unused.fa"
OUTPUT_DIR="./blast"
BATCH_SIZE=100
THREADS=12

# 创建输出目录
mkdir -p $OUTPUT_DIR/tmp_fasta

# 分割FASTA文件，每个文件包含BATCH_SIZE条序列
echo "正在分割FASTA文件..."

# 使用awk正确分割FASTA文件
awk -v size=$BATCH_SIZE 'BEGIN {
    n = 0
    file_num = 0
}
/^>/ {
    if (n % size == 0 && n > 0) {
        close(file)
        file_num++
    }
    if (n % size == 0) {
        file = sprintf("'$OUTPUT_DIR'/tmp_fasta/batch_%04d.fa", file_num)
    }
    print > file
    n++
    next
}
{
    print > file
}
END {
    close(file)
    print "总共分割成 " (file_num + 1) " 个批次文件"
}' $INPUT_FASTA

# 统计分割后的文件数量
BATCH_NUM=$(ls -1 $OUTPUT_DIR/tmp_fasta/*.fa 2>/dev/null | wc -l)
echo "共分割成 $BATCH_NUM 个批次文件"

# 检查分割文件是否有内容
echo "检查分割文件..."
for f in $OUTPUT_DIR/tmp_fasta/*.fa; do
    seq_count=$(grep -c "^>" "$f")
    size=$(ls -lh "$f" | awk '{print $5}')
    echo "$(basename $f): $seq_count 条序列, 大小: $size"
done

# 环境检查
echo "=== 环境检查 ==="
echo "blastn 路径: $(which blastn)"
echo "blastn 版本: $(blastn -version 2>&1 | head -n 1)"
echo "当前目录: $(pwd)"
echo "HOME: $HOME"

# 批量运行BLAST
echo "开始批量BLAST搜索..."
for batch_file in $OUTPUT_DIR/tmp_fasta/*.fa; do
    batch_name=$(basename $batch_file .fa)
    echo "========================================"
    echo "正在处理: $batch_name"
    echo "命令: blastn -query $batch_file -db ~/database/nt/nt -out $OUTPUT_DIR/${batch_name}.blast.txt -outfmt 6 -evalue 1e-5 -num_threads $THREADS"

    blastn -query "$batch_file" \
           -db ~/database/nt/nt \
           -out "$OUTPUT_DIR/${batch_name}.blast.txt" \
           -outfmt 6 \
           -evalue 1e-5 \
           -num_threads $THREADS

    # 检查BLAST命令是否成功
    blast_exit_code=$?
    echo "BLAST 退出码: $blast_exit_code"

    # 检查输出文件
    if [ -s "$OUTPUT_DIR/${batch_name}.blast.txt" ]; then
        result_lines=$(wc -l < "$OUTPUT_DIR/${batch_name}.blast.txt")
        echo "✓ 完成: $batch_name, 结果行数: $result_lines"
    else
        echo "✗ 警告: $batch_name 输出文件为空！"
        echo "输出文件路径: $OUTPUT_DIR/${batch_name}.blast.txt"
    fi
done

# 合并所有结果
echo "正在合并所有BLAST结果..."
cat $OUTPUT_DIR/batch_*.blast.txt > $OUTPUT_DIR/41.unused.nt.blast.txt

# 清理临时文件
echo "清理临时文件..."
rm -rf $OUTPUT_DIR/tmp_fasta
# 保留批次结果文件，如果不需要可以删除
# rm -f $OUTPUT_DIR/batch_*.blast.txt

echo "全部完成！"
echo "最终结果文件: $OUTPUT_DIR/41.unused.nt.blast.txt"
