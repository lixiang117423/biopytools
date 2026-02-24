#!/bin/bash

# 脚本用途：使用 RagTag 进行基因组错误校正
# 作者：Claude Code
# 日期：2025-12-23

set -euo pipefail

# 工作目录
WORK_DIR="/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/32.ragtag校正"
cd "$WORK_DIR"

# 输入文件
REFERENCE="${WORK_DIR}/ref.fa"
QUERY="${WORK_DIR}/OV53_1.primary_corrected.chr.fa"

# 输出目录
OUTPUT_DIR="${WORK_DIR}/ragtag_correct_output"
mkdir -p "$OUTPUT_DIR"

# 日志文件
LOG_FILE="${WORK_DIR}/ragtag_correct.log"
exec > >(tee -a "$LOG_FILE") 2>&1

echo "=========================================="
echo "RagTag 基因组校正"
echo "开始时间: $(date '+%Y-%m-%d %H:%M:%S')"
echo "=========================================="
echo ""
echo "参考序列: $REFERENCE"
echo "待校正序列: $QUERY"
echo "输出目录: $OUTPUT_DIR"
echo ""

# 检查输入文件
if [ ! -f "$REFERENCE" ]; then
    echo "错误: 参考序列文件不存在: $REFERENCE"
    exit 1
fi

if [ ! -f "$QUERY" ]; then
    echo "错误: 待校正序列文件不存在: $QUERY"
    exit 1
fi

# 统计输入序列信息
echo "参考序列统计:"
seqkit stats "$REFERENCE"
echo ""
echo "待校正序列统计:"
seqkit stats "$QUERY"
echo ""

# 运行 RagTag correct
echo "=========================================="
echo "开始运行 RagTag correct..."
echo "=========================================="

ragtag.py correct \
    -t 16 \
    -o "$OUTPUT_DIR" \
    -u \
    -w \
    "$REFERENCE" \
    "$QUERY"

echo ""
echo "=========================================="
echo "RagTag correct 完成"
echo "结束时间: $(date '+%Y-%m-%d %H:%M:%S')"
echo "=========================================="
echo ""

# 统计输出结果
if [ -f "${OUTPUT_DIR}/query.corrected.fa" ]; then
    echo "校正后的序列统计:"
    seqkit stats "${OUTPUT_DIR}/query.corrected.fa"
    echo ""
    echo "校正序列文件: ${OUTPUT_DIR}/query.corrected.fa"
fi

# 比较校正前后的序列数量
echo "=========================================="
echo "校正前后对比:"
echo "=========================================="
ref_count=$(grep -c "^>" "$REFERENCE" || echo 0)
query_count=$(grep -c "^>" "$QUERY" || echo 0)
corrected_count=$(grep -c "^>" "${OUTPUT_DIR}/query.corrected.fa" 2>/dev/null || echo 0)

printf "参考序列数:    %d\n" "$ref_count"
printf "原始序列数:    %d\n" "$query_count"
printf "校正后序列数:  %d\n" "$corrected_count"
printf "新增断裂数:    %d\n" $((corrected_count - query_count))
echo ""
