#!/bin/bash

# 设置输入目录和输出目录
FASTA_DIR=~/database/brassicaceae/fasta
OUTPUT_DIR=$(pwd)

# 遍历所有 .fa 文件
for fa_file in "$FASTA_DIR"/*.fa; do
    # 获取文件名（不带路径和扩展名）
    basename=$(basename "$fa_file" .fa)

    echo "Processing $basename..."

    # 执行 wgsim 命令
    wgsim "$fa_file" -N 50000000 -1 150 -2 150 \
        "${OUTPUT_DIR}/${basename}_1.fq.gz" \
        "${OUTPUT_DIR}/${basename}_2.fq.gz"

    # 检查命令是否成功
    if [ $? -eq 0 ]; then
        echo "Successfully processed $basename"
    else
        echo "Error processing $basename" >&2
    fi
done

echo "All files processed!"
