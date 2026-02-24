#!/bin/bash
#
# 统计Chr12:135000000到末尾的碱基覆盖度
# 输出格式: 第一列是位置，后面每列是一个样品的覆盖度
#

set -e

# 定义路径
BAM_DIR="/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/31.重测序数据比对到挂载的结果上/02.mapping/bam"
OUTPUT_DIR="/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/31.重测序数据比对到挂载的结果上/02.mapping"
REGION="Chr12:135000000-"
OUTPUT="${OUTPUT_DIR}/Chr12_135M_to_end_coverage.txt"

# 临时目录
TEMP_DIR="${OUTPUT_DIR}/temp_coverage"
mkdir -p "$TEMP_DIR"

echo "=========================================="
echo "BAM文件覆盖度统计"
echo "=========================================="
echo "BAM文件夹: $BAM_DIR"
echo "目标区域: $REGION"
echo "输出文件: $OUTPUT"
echo "------------------------------------------"

# 检查bam文件夹是否存在
if [ ! -d "$BAM_DIR" ]; then
    echo "错误: BAM文件夹不存在: $BAM_DIR"
    exit 1
fi

# 检查samtools是否可用
if ! command -v samtools &> /dev/null; then
    echo "错误: 未找到samtools，请先安装samtools"
    exit 1
fi

# 获取所有bam文件（包括.bam和.cram）
echo "正在查找BAM文件..."
BAM_FILES=()
for bam in "$BAM_DIR"/*.bam "$BAM_DIR"/*.cram; do
    if [ -f "$bam" ]; then
        BAM_FILES+=("$bam")
        echo "  找到: $(basename "$bam")"
    fi
done

if [ ${#BAM_FILES[@]} -eq 0 ]; then
    echo "错误: 在 $BAM_DIR 中未找到BAM/CRAM文件"
    exit 1
fi

echo "------------------------------------------"
echo "共找到 ${#BAM_FILES[@]} 个BAM文件"
echo "开始提取覆盖度信息..."
echo "------------------------------------------"

# 提取每个样品的覆盖度
SAMPLE_NAMES=()
INDEX=0

for bam in "${BAM_FILES[@]}"; do
    BASENAME=$(basename "$bam")
    # 去掉.bam或.cram后缀作为样品名
    SAMPLE_NAME="${BASENAME%.bam}"
    SAMPLE_NAME="${SAMPLE_NAME%.cram}"
    SAMPLE_NAMES+=("$SAMPLE_NAME")

    TEMP_FILE="${TEMP_DIR}/${SAMPLE_NAME}_depth.txt"

    echo "[$((INDEX+1))/${#BAM_FILES[@]}] 处理: $SAMPLE_NAME"

    # 使用samtools depth获取覆盖度
    samtools depth -r "$REGION" "$bam" > "$TEMP_FILE"

    INDEX=$((INDEX+1))
done

echo "------------------------------------------"
echo "合并所有样品的覆盖度数据..."

# 获取第一个文件作为参考
FIRST_TEMP="${TEMP_DIR}/${SAMPLE_NAMES[0]}_depth.txt"

# 检查第一个文件是否有数据
if [ ! -s "$FIRST_TEMP" ]; then
    echo "警告: 目标区域 ($REGION) 没有覆盖度数据"
    echo "可能的原因:"
    echo "  1. 染色体名称不是 'Chr12' (注意大小写)"
    echo "  2. 起始位置超出染色体长度"
    echo ""
    echo "正在检查可用的染色体..."
    samtools view -H "${BAM_FILES[0]}" | grep "^@SQ" | cut -f2 | sed 's/SN://' | sort
    exit 1
fi

# 创建输出文件的header
HEADER="Position"
for sample in "${SAMPLE_NAMES[@]}"; do
    HEADER="${HEADER}\t${sample}"
done
echo -e "$HEADER" > "$OUTPUT"

# 使用paste命令合并所有文件
# 首先获取所有临时文件
TEMP_FILES=()
for sample in "${SAMPLE_NAMES[@]}"; do
    TEMP_FILES+=("${TEMP_DIR}/${sample}_depth.txt")
done

# 合并文件 (只保留位置列和覆盖度列)
# paste会按列合并，我们需要：
# 从第一个文件取第1列(位置)，从每个文件取第3列(覆盖度)
paste "${TEMP_FILES[@]}" | \
    awk -v n Samples=${#SAMPLE_NAMES[@]} '{
        # 第1列是位置(从第一个文件)
        printf "%s", $1;
        # 然后每个文件的覆盖度(第3列)
        for (i=1; i<=n; i++) {
            col = i*3 - 1;  # 计算覆盖度所在的列号
            printf "\t%s", $col;
        }
        printf "\n";
    }' >> "$OUTPUT"

echo "------------------------------------------"
echo "清理临时文件..."
rm -rf "$TEMP_DIR"

echo "------------------------------------------"
echo "完成！"
echo "输出文件: $OUTPUT"
echo ""
echo "统计信息:"
TOTAL_POS=$(tail -n +2 "$OUTPUT" | wc -l)
echo "  总碱基数: $TOTAL_POS"
echo "  样品数: ${#SAMPLE_NAMES[@]}"
echo ""
echo "前10行预览:"
head -n 10 "$OUTPUT"
echo "=========================================="
