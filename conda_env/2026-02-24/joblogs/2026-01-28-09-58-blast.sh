#!/bin/bash

# 设置变量
INPUT_FASTA="41.unused.fa"
OUTPUT_DIR="./blast"
DB="~/database/nt/nt"
BATCH_SIZE=100
THREADS=12

# 创建输出目录
mkdir -p $OUTPUT_DIR/tmp_fasta

# 分割FASTA文件，每个文件包含BATCH_SIZE条序列
echo "正在分割FASTA文件..."
python3 - <<EOF
from Bio import SeqIO
import os

batch_size = $BATCH_SIZE
output_dir = "$OUTPUT_DIR/tmp_fasta"
input_file = "$INPUT_FASTA"

batch_num = 0
record_count = 0
current_batch = []

for record in SeqIO.parse(input_file, "fasta"):
    current_batch.append(record)
    record_count += 1

    if len(current_batch) >= batch_size:
        output_file = os.path.join(output_dir, f"batch_{batch_num:04d}.fa")
        with open(output_file, "w") as f:
            SeqIO.write(current_batch, f, "fasta")
        print(f"写入批次 {batch_num}: {len(current_batch)} 条序列")
        batch_num += 1
        current_batch = []

# 写入剩余的序列
if current_batch:
    output_file = os.path.join(output_dir, f"batch_{batch_num:04d}.fa")
    with open(output_file, "w") as f:
        SeqIO.write(current_batch, f, "fasta")
    print(f"写入批次 {batch_num}: {len(current_batch)} 条序列")
    batch_num += 1

print(f"总共分割成 {batch_num} 个批次文件")
EOF

# 统计分割后的文件数量
BATCH_NUM=$(ls -1 $OUTPUT_DIR/tmp_fasta/*.fa 2>/dev/null | wc -l)
echo "共分割成 $BATCH_NUM 个批次文件"

# 批量运行BLAST
echo "开始批量BLAST搜索..."
for batch_file in $OUTPUT_DIR/tmp_fasta/*.fa; do
    batch_name=$(basename $batch_file .fa)
    echo "正在处理: $batch_name"

    blastn -query $batch_file \
           -db $DB \
           -out $OUTPUT_DIR/${batch_name}.blast.txt \
           -outfmt 6 \
           -evalue 1e-5 \
           -num_threads $THREADS

    echo "完成: $batch_name"
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
