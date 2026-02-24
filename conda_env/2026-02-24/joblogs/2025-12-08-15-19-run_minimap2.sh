#!/bin/bash

# ================= 配置区域 =================
# 参考基因组 (旧版本, X轴)
REF="genome_chr.fa"
# 你的组装结果 (Hap1, Y轴)
QRY="OV53_YaHS_Chr.fa"
# 输出前缀
PREFIX="OV53_vs_Ref"
# 线程数
THREADS=64
# ===========================================

# 1. 检查输入文件
if [[ ! -f "$REF" || ! -f "$QRY" ]]; then
    echo "Error: 输入文件不存在，请检查文件名！"
    echo "Ref: $REF"
    echo "Query: $QRY"
    exit 1
fi

echo "=========================================="
echo "Step 1: Running Minimap2 Alignment..."
echo "=========================================="

# 运行 minimap2
# -x asm5: 适用于同物种比较 (divergence < 5%)
# -t: 线程
# 输出 .paf 格式
if [[ ! -f "${PREFIX}.paf" ]]; then
    minimap2 -x asm5 -t $THREADS "$REF" "$QRY" > "${PREFIX}.paf"
    echo "Alignment finished: ${PREFIX}.paf generated."
else
    echo "PAF file already exists, skipping alignment."
fi
