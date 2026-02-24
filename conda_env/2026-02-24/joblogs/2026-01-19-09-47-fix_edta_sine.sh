#!/bin/bash
# 手动修复EDTA SINE缺失问题的脚本|Manual fix script for EDTA SINE missing issue

set -e

GENOME="/share/org/YZWL/yzwl_lixg/tmp/test_edta/OV53_Chr.fa"
OUTPUT_DIR="/share/org/YZWL/yzwl_lixg/tmp/test_edta/output"
THREADS=64

cd "$OUTPUT_DIR"

# 获取基因组文件名（不含路径和扩展名）
GENOME_NAME=$(basename "$GENOME" .fa)

echo "========================================="
echo "修复EDTA SINE缺失问题|Fixing EDTA SINE missing issue"
echo "========================================="
echo "基因组|Genome: $GENOME_NAME"
echo "输出目录|Output: $OUTPUT_DIR"
echo ""

# 创建空的SINE文件
echo "步骤1: 创建空SINE文件|Step 1: Create empty SINE file"
SINE_FILE="$OUTPUT_DIR/${GENOME_NAME}.mod.EDTA.raw/${GENOME_NAME}.mod.SINE.raw.fa"
mkdir -p "$(dirname "$SINE_FILE")"
touch "$SINE_FILE"
echo "已创建|Created: $SINE_FILE"
echo ""

# 继续运行EDTA filter步骤
echo "步骤2: 使用--step filter继续运行EDTA|Step 2: Continue EDTA with --step filter"
echo "这可能需要几个小时|This may take several hours"
echo ""

perl /share/org/YZWL/yzwl_lixg/miniforge3/envs/EDTA_v.2.2.2/share/EDTA/EDTA.pl \
    --genome "$GENOME" \
    --species others \
    --step filter \
    --overwrite 0 \
    --maxdiv 40 \
    --u 1.3e-08 \
    --threads $THREADS \
    --debug 0 \
    --anno 1

echo ""
echo "========================================="
echo "EDTA分析完成|EDTA analysis completed"
echo "========================================="
