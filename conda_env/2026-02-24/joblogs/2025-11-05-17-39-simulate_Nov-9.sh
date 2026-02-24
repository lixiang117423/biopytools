#!/bin/bash

# --- 设置参数 ---

# 1. 参考基因组文件路径
GENOME_FILE="./GCA_009729435.1_ASM972943v1_genomic.fna"

# 2. 输出的 Reads 1 文件名
OUTPUT_R1="Nov-9_1.clean.fq"

# 3. 输出的 Reads 2 文件名
OUTPUT_R2="Nov-9_2.clean.fq"

# 4. 模拟的深度
DEPTH=100

# 5. Reads 长度
READ_LENGTH=250

# 6. 双末端测序的插入片段平均长度
FRAGMENT_SIZE=500

# --- 检查文件是否存在 ---
if [ ! -f "$GENOME_FILE" ]; then
    echo "错误: 基因组文件未找到: $GENOME_FILE"
    exit 1
fi

# --- 计算基因组大小 ---
echo "正在计算基因组大小..."
GENOME_SIZE=$(grep -v ">" "$GENOME_FILE" | wc -c | awk '{print $1}')
echo "基因组大小为: $GENOME_SIZE bp"

# --- 计算需要生成的 reads 数量 ---
# N = (基因组大小 * 深度) / (Reads 长度 * 2)
NUM_READS=$(echo "($GENOME_SIZE * $DEPTH) / ($READ_LENGTH * 2)" | bc)
echo "为达到 ${DEPTH}x 深度，需要生成 $NUM_READS 对 reads"

# --- 运行 wgsim 命令 ---
echo "开始使用 wgsim 模拟重测序 reads..."
echo "-------------------------------------"

wgsim \
    -N "$NUM_READS" \
    -1 "$READ_LENGTH" \
    -2 "$READ_LENGTH" \
    -d "$FRAGMENT_SIZE" \
    -e 0.001 \
    -r 0.0001 \
    "$GENOME_FILE" \
    "$OUTPUT_R1" \
    "$OUTPUT_R2"

echo "-------------------------------------"
echo "模拟完成！"
echo "生成的文件: $OUTPUT_R1 和 $OUTPUT_R2"
