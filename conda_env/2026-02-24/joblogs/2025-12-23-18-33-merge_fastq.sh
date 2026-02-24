#!/bin/bash

# 脚本用途：合并双端测序数据
# 作者：Claude Code
# 日期：2025-12-23
# 说明：将每个样品的所有 _1.fq.gz 和 _2.fq.gz 文件分别合并

set -euo pipefail

# 当前脚本所在目录
SCRIPT_DIR="/share/org/YZWL/yzwl_lixg/project/21.野生大豆群体/01.data/raw"
cd "$SCRIPT_DIR"

# 输出目录
OUTPUT_DIR="${SCRIPT_DIR}/merged"
mkdir -p "$OUTPUT_DIR"

# 日志文件
LOG_FILE="${SCRIPT_DIR}/merge_fastq.log"
exec > >(tee -a "$LOG_FILE") 2>&1

echo "=========================================="
echo "开始合并FASTQ文件"
echo "开始时间: $(date '+%Y-%m-%d %H:%M:%S')"
echo "=========================================="

# 遍历upload目录下的所有批次目录
for batch_dir in upload/F25A040008383_*; do
    if [ ! -d "$batch_dir" ]; then
        continue
    fi

    echo ""
    echo "处理批次目录: $batch_dir"

    # 遍历每个样品目录
    for sample_dir in "$batch_dir"/*; do
        # 跳过非目录文件（如Clean_stat.xls）
        if [ ! -d "$sample_dir" ]; then
            continue
        fi

        # 获取样品名称
        sample_name=$(basename "$sample_dir")

        # 检查是否有fastq.gz文件
        if ! ls "$sample_dir"/*_1.fq.gz &>/dev/null; then
            echo "  跳过 $sample_name: 没有找到 _1.fq.gz 文件"
            continue
        fi

        echo "  处理样品: $sample_name"

        # 合并 _1.fq.gz 文件
        output_1="${OUTPUT_DIR}/${sample_name}_1.fq.gz"
        if [ -f "$output_1" ]; then
            echo "    警告: $output_1 已存在，将覆盖"
            rm -f "$output_1"
        fi

        echo "    合并 _1.fq.gz 文件..."
        cat "$sample_dir"/*_1.fq.gz > "$output_1"
        original_size=$(du -ch "$sample_dir"/*_1.fq.gz | tail -1 | cut -f1)
        merged_size=$(du -h "$output_1" | cut -f1)
        echo "    原始大小: $original_size -> 合并后: $merged_size"

        # 合并 _2.fq.gz 文件
        output_2="${OUTPUT_DIR}/${sample_name}_2.fq.gz"
        if [ -f "$output_2" ]; then
            echo "    警告: $output_2 已存在，将覆盖"
            rm -f "$output_2"
        fi

        echo "    合并 _2.fq.gz 文件..."
        cat "$sample_dir"/*_2.fq.gz > "$output_2"
        original_size=$(du -ch "$sample_dir"/*_2.fq.gz | tail -1 | cut -f1)
        merged_size=$(du -h "$output_2" | cut -f1)
        echo "    原始大小: $original_size -> 合并后: $merged_size"
    done
done

echo ""
echo "=========================================="
echo "合并完成"
echo "结束时间: $(date '+%Y-%m-%d %H:%M:%S')"
echo "输出目录: $OUTPUT_DIR"
echo "=========================================="

# 统计合并后的文件数量
total_files=$(ls -1 "$OUTPUT_DIR"/*.fq.gz 2>/dev/null | wc -l)
echo "共生成 $total_files 个FASTQ文件"

# 显示合并后的文件列表
echo ""
echo "合并后的文件列表:"
ls -lh "$OUTPUT_DIR"/*.fq.gz | awk '{print $9, $5}'
