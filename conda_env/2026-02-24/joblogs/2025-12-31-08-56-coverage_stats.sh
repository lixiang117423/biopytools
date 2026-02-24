#!/bin/bash

# 完整的Chr12覆盖度提取和合并脚本 - 修正版

# ==================== 设置参数 ====================
BAM_DIR="/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/31.重测序数据比对到挂载的结果上/02.mapping/bam"
CHR="Chr12"
START=125000000
# 使用绝对路径
OUTPUT_DIR="$(pwd)/chr12_coverage"

# ==================== 创建输出目录 ====================
echo "Creating output directories..."
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${OUTPUT_DIR}/temp"

# 验证目录创建成功
if [ ! -d "${OUTPUT_DIR}/temp" ]; then
    echo "Error: Failed to create output directory ${OUTPUT_DIR}/temp"
    exit 1
fi

echo "============================================"
echo "Chr12 Coverage Analysis Pipeline"
echo "============================================"
echo "BAM directory: ${BAM_DIR}"
echo "Output directory: ${OUTPUT_DIR}"
echo "Target region: ${CHR}:${START}-end"
echo "============================================"

# ==================== 步骤1: 获取Chr12长度 ====================
echo ""
echo "[Step 1] Getting Chr12 length..."

if [ ! -d "${BAM_DIR}" ]; then
    echo "Error: BAM directory does not exist: ${BAM_DIR}"
    exit 1
fi

cd "${BAM_DIR}" || exit 1

# 从第一个BAM文件获取染色体长度
FIRST_BAM=$(ls *.bam 2>/dev/null | head -1)
if [ -z "${FIRST_BAM}" ]; then
    echo "Error: No BAM files found in ${BAM_DIR}"
    exit 1
fi

CHR_LENGTH=$(samtools view -H ${FIRST_BAM} | grep -w "SN:${CHR}" | awk '{print $3}' | sed 's/LN://')

if [ -z "${CHR_LENGTH}" ]; then
    echo "Error: Could not find ${CHR} in BAM header"
    exit 1
fi

echo "Chr12 length: ${CHR_LENGTH} bp"
echo "Target region: ${CHR}:${START}-${CHR_LENGTH}"
REGION_SIZE=$((CHR_LENGTH - START + 1))
echo "Region size: ${REGION_SIZE} bp"

# ==================== 步骤2: 为每个BAM文件提取覆盖度 ====================
echo ""
echo "[Step 2] Extracting coverage for each BAM file..."

BAM_COUNT=$(ls *.bam 2>/dev/null | wc -l)
echo "Found ${BAM_COUNT} BAM files"

# 确保在BAM目录中
cd "${BAM_DIR}" || exit 1

CURRENT=0
SUCCESS_COUNT=0
FAIL_COUNT=0

for bam in *.bam; do
    CURRENT=$((CURRENT + 1))
    sample=$(basename ${bam} .bam)
    
    # 输出文件使用绝对路径
    OUTPUT_FILE="${OUTPUT_DIR}/temp/${sample}.chr12.depth"
    
    echo "  [${CURRENT}/${BAM_COUNT}] Processing ${sample}..."
    
    # 检查BAM索引
    if [ ! -f "${bam}.bai" ]; then
        echo "    Creating index for ${bam}..."
        samtools index "${bam}"
    fi
    
    # 提取覆盖度 - 去掉 -a 参数以节省空间（只输出有覆盖的位点）
    samtools depth -r ${CHR}:${START}-${CHR_LENGTH} -Q 0 -q 0 "${bam}" > "${OUTPUT_FILE}" 2>&1
    
    # 检查命令执行状态
    if [ $? -ne 0 ]; then
        echo "    ERROR: samtools depth failed for ${sample}"
        FAIL_COUNT=$((FAIL_COUNT + 1))
        continue
    fi
    
    # 检查输出文件
    if [ ! -f "${OUTPUT_FILE}" ]; then
        echo "    ERROR: Output file not created: ${OUTPUT_FILE}"
        FAIL_COUNT=$((FAIL_COUNT + 1))
    elif [ ! -s "${OUTPUT_FILE}" ]; then
        echo "    WARNING: Output file is empty for ${sample}"
        FAIL_COUNT=$((FAIL_COUNT + 1))
    else
        LINE_COUNT=$(wc -l < "${OUTPUT_FILE}")
        FILE_SIZE=$(du -h "${OUTPUT_FILE}" | cut -f1)
        AVG_COV=$(awk '{sum+=$3} END {if(NR>0) printf "%.2f", sum/NR; else print "0"}' "${OUTPUT_FILE}")
        echo "    SUCCESS: ${LINE_COUNT} positions, ${FILE_SIZE}, Avg: ${AVG_COV}X"
        SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
    fi
done

echo ""
echo "Extraction summary: ${SUCCESS_COUNT} succeeded, ${FAIL_COUNT} failed"

if [ ${SUCCESS_COUNT} -eq 0 ]; then
    echo "Error: No samples were successfully processed!"
    exit 1
fi

# ==================== 步骤3: 合并所有样本的覆盖度数据 ====================
echo ""
echo "[Step 3] Merging all samples..."

cd "${OUTPUT_DIR}/temp" || exit 1

# 获取所有非空的depth文件
VALID_FILES=($(ls *.chr12.depth 2>/dev/null | xargs -I{} bash -c '[ -s {} ] && echo {}'))

if [ ${#VALID_FILES[@]} -eq 0 ]; then
    echo "Error: No valid depth files to merge"
    exit 1
fi

echo "Merging ${#VALID_FILES[@]} valid files..."

# 创建表头
printf "Chr\tPos" > "${OUTPUT_DIR}/merged_chr12_coverage.txt"
for file in "${VALID_FILES[@]}"; do
    sample=$(basename ${file} .chr12.depth)
    printf "\t%s" "${sample}" >> "${OUTPUT_DIR}/merged_chr12_coverage.txt"
done
printf "\n" >> "${OUTPUT_DIR}/merged_chr12_coverage.txt"

# 使用 paste 和 awk 合并（更高效）
echo "  Merging data..."

# 方案：使用所有文件生成统一的位置列表，然后逐个添加样本数据
# 首先获取所有唯一位置
echo "  Step 3.1: Creating position list..."
cat "${VALID_FILES[@]}" | cut -f1,2 | sort -u -k2,2n > "${OUTPUT_DIR}/all_positions.txt"

TOTAL_POS=$(wc -l < "${OUTPUT_DIR}/all_positions.txt")
echo "  Total unique positions: ${TOTAL_POS}"

# 为每个样本创建完整的位置-覆盖度映射
echo "  Step 3.2: Processing each sample..."
for i in "${!VALID_FILES[@]}"; do
    file="${VALID_FILES[$i]}"
    sample=$(basename ${file} .chr12.depth)
    echo "    Processing ${sample} ($((i+1))/${#VALID_FILES[@]})..."
    
    # 使用 join 合并位置列表和该样本的覆盖度
    # 如果位置不存在，填充0
    awk 'NR==FNR {depth[$1"\t"$2]=$3; next} 
         {key=$1"\t"$2; print (key in depth) ? depth[key] : 0}' \
         "${file}" "${OUTPUT_DIR}/all_positions.txt" > "${OUTPUT_DIR}/temp_${sample}.txt"
done

echo "  Step 3.3: Combining all columns..."
# 合并所有列
paste "${OUTPUT_DIR}/all_positions.txt" ${OUTPUT_DIR}/temp_*.txt >> "${OUTPUT_DIR}/merged_chr12_coverage.txt"

# 清理临时文件
rm -f "${OUTPUT_DIR}/all_positions.txt" "${OUTPUT_DIR}/temp_*.txt"

echo "  Merge completed!"

# ==================== 步骤4: 生成统计摘要 ====================
echo ""
echo "[Step 4] Generating summary statistics..."

cd "${OUTPUT_DIR}" || exit 1

printf "Sample\tTotal_Positions\tMean_Coverage\tMedian_Coverage\tMax_Coverage\tPositions_0X\tPositions_>0X\tPositions_>10X\tPositions_>30X\tCoverage_%%_>0X\tCoverage_%%_>10X\n" > coverage_summary.txt

for file in temp/*.chr12.depth; do
    if [ ! -s "$file" ]; then
        continue
    fi
    
    sample=$(basename ${file} .chr12.depth)
    
    awk -v sample="${sample}" '
    {
        sum += $3
        count++
        depths[count] = $3
        
        if ($3 == 0) zero++
        if ($3 > 0) gt0++
        if ($3 > 10) gt10++
        if ($3 > 30) gt30++
        if ($3 > max) max = $3
    }
    END {
        mean = (count > 0) ? sum/count : 0
        
        n = asort(depths)
        if (n % 2 == 1) {
            median = depths[int(n/2) + 1]
        } else {
            median = (depths[n/2] + depths[n/2 + 1]) / 2
        }
        
        pct_gt0 = (count > 0) ? (gt0 / count) * 100 : 0
        pct_gt10 = (count > 0) ? (gt10 / count) * 100 : 0
        
        printf "%s\t%d\t%.2f\t%.2f\t%d\t%d\t%d\t%d\t%d\t%.2f\t%.2f\n", 
               sample, count, mean, median, max, zero, gt0, gt10, gt30, pct_gt0, pct_gt10
    }' "${file}" >> coverage_summary.txt
done

# ==================== 完成 ====================
MERGED_LINES=$(wc -l < merged_chr12_coverage.txt)
MERGED_SIZE=$(du -h merged_chr12_coverage.txt | cut -f1)

echo ""
echo "============================================"
echo "Analysis Complete!"
echo "============================================"
echo "Region: ${CHR}:${START}-${CHR_LENGTH} (${REGION_SIZE} bp)"
echo "Samples processed: ${SUCCESS_COUNT}/${BAM_COUNT}"
echo "Positions in merged file: $((MERGED_LINES - 1))"
echo ""
echo "Output files:"
echo "  1. ${OUTPUT_DIR}/merged_chr12_coverage.txt (${MERGED_SIZE})"
echo "  2. ${OUTPUT_DIR}/coverage_summary.txt"
echo "  3. ${OUTPUT_DIR}/temp/*.chr12.depth"
echo ""
echo "Preview:"
head -5 merged_chr12_coverage.txt | cut -f1-6
echo ""
echo "============================================"
