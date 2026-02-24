#!/bin/bash

# 完整的Chr12覆盖度提取和合并脚本
# 使用 samtools depth 方法（更可靠）

# ==================== 设置参数 ====================
BAM_DIR="/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/31.重测序数据比对到挂载的结果上/02.mapping/bam"
CHR="Chr12"
START=125000000
OUTPUT_DIR="./chr12_coverage"
THREADS=8  # 可根据服务器资源调整

# ==================== 创建输出目录 ====================
mkdir -p ${OUTPUT_DIR}
mkdir -p ${OUTPUT_DIR}/temp

echo "============================================"
echo "Chr12 Coverage Analysis Pipeline"
echo "============================================"
echo "BAM directory: ${BAM_DIR}"
echo "Target region: ${CHR}:${START}-end"
echo "Output directory: ${OUTPUT_DIR}"
echo "============================================"

# ==================== 步骤1: 获取Chr12长度 ====================
echo ""
echo "[Step 1] Getting Chr12 length..."
cd ${BAM_DIR}

# 从第一个BAM文件获取染色体长度
FIRST_BAM=$(ls *.bam | head -1)
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

BAM_COUNT=$(ls *.bam | wc -l)
echo "Found ${BAM_COUNT} BAM files"

CURRENT=0
for bam in *.bam; do
    CURRENT=$((CURRENT + 1))
    sample=$(basename ${bam} .bam)
    echo "  [${CURRENT}/${BAM_COUNT}] Processing ${sample}..."
    
    # 检查BAM索引是否存在，不存在则创建
    if [ ! -f "${bam}.bai" ]; then
        echo "    Creating index for ${bam}..."
        samtools index ${bam}
    fi
    
    # 使用 samtools depth 直接提取指定区域
    # -r 指定区域
    # -a 输出所有位点（包括0覆盖度的位点）
    # -Q 0 不过滤mapping quality
    # -q 0 不过滤base quality
    samtools depth -r ${CHR}:${START}-${CHR_LENGTH} -a -Q 0 -q 0 ${bam} \
    > ${OUTPUT_DIR}/temp/${sample}.chr12.depth
    
    # 检查输出文件
    if [ ! -s ${OUTPUT_DIR}/temp/${sample}.chr12.depth ]; then
        echo "    Warning: No data extracted for ${sample}"
    else
        LINE_COUNT=$(wc -l < ${OUTPUT_DIR}/temp/${sample}.chr12.depth)
        # 快速统计平均覆盖度
        AVG_COV=$(awk '{sum+=$3} END {printf "%.2f", sum/NR}' ${OUTPUT_DIR}/temp/${sample}.chr12.depth)
        echo "    Extracted ${LINE_COUNT} positions, Avg coverage: ${AVG_COV}X"
    fi
done

# ==================== 步骤3: 检查提取的文件 ====================
echo ""
echo "[Step 3] Checking extracted files..."

DEPTH_FILES=(${OUTPUT_DIR}/temp/*.chr12.depth)
DEPTH_COUNT=${#DEPTH_FILES[@]}

if [ ${DEPTH_COUNT} -eq 0 ]; then
    echo "Error: No depth files generated"
    exit 1
fi

echo "Successfully generated ${DEPTH_COUNT} depth files"

# 检查是否有空文件
EMPTY_COUNT=0
for file in ${OUTPUT_DIR}/temp/*.chr12.depth; do
    if [ ! -s "$file" ]; then
        EMPTY_COUNT=$((EMPTY_COUNT + 1))
    fi
done

if [ ${EMPTY_COUNT} -gt 0 ]; then
    echo "Warning: ${EMPTY_COUNT} files are empty"
fi

# ==================== 步骤4: 合并所有样本的覆盖度数据 ====================
echo ""
echo "[Step 4] Merging all samples..."
echo "This may take a few minutes for large datasets..."

cd ${OUTPUT_DIR}/temp

# 获取所有非空的depth文件
VALID_FILES=($(ls -S *.chr12.depth 2>/dev/null | xargs -I{} bash -c '[ -s {} ] && echo {}'))

if [ ${#VALID_FILES[@]} -eq 0 ]; then
    echo "Error: No valid depth files to merge"
    exit 1
fi

echo "Merging ${#VALID_FILES[@]} valid files..."

# 创建表头
echo -n "Chr"$'\t'"Pos" > ../merged_chr12_coverage.txt
for file in "${VALID_FILES[@]}"; do
    sample=$(basename ${file} .chr12.depth)
    echo -n $'\t'"${sample}" >> ../merged_chr12_coverage.txt
done
echo "" >> ../merged_chr12_coverage.txt

# 使用改进的 awk 脚本合并数据
echo "  Processing with awk (this will take time)..."

awk '
BEGIN {
    OFS="\t"
}
FNR==1 {
    # 新文件开始，增加文件计数
    file_idx++
    file_name[file_idx] = FILENAME
}
{
    # $1=Chr, $2=Pos, $3=Depth
    pos = $2
    
    # 第一次遇到这个位置时，记录染色体信息
    if (!(pos in chr_info)) {
        chr_info[pos] = $1
        pos_list[++pos_count] = pos
    }
    
    # 存储该位置该样本的覆盖度
    depth[pos, file_idx] = $3
}
END {
    # 输出所有位置的数据
    for (i = 1; i <= pos_count; i++) {
        pos = pos_list[i]
        printf "%s\t%s", chr_info[pos], pos
        
        # 输出每个文件的覆盖度
        for (j = 1; j <= file_idx; j++) {
            if ((pos, j) in depth) {
                printf "\t%s", depth[pos, j]
            } else {
                printf "\t0"
            }
        }
        printf "\n"
    }
}' "${VALID_FILES[@]}" >> ../merged_chr12_coverage.txt

echo "  Merge completed!"

# ==================== 步骤5: 生成统计摘要 ====================
echo ""
echo "[Step 5] Generating summary statistics..."

cd ${OUTPUT_DIR}

# 创建统计摘要文件
echo -e "Sample\tTotal_Positions\tMean_Coverage\tMedian_Coverage\tMax_Coverage\tPositions_0X\tPositions_>0X\tPositions_>10X\tPositions_>30X\tCoverage_%_>0X\tCoverage_%_>10X" > coverage_summary.txt

for file in temp/*.chr12.depth; do
    # 跳过空文件
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
        
        # 计算中位数 - 排序
        n = asort(depths)
        if (n % 2 == 1) {
            median = depths[int(n/2) + 1]
        } else {
            median = (depths[n/2] + depths[n/2 + 1]) / 2
        }
        
        # 计算百分比
        pct_gt0 = (count > 0) ? (gt0 / count) * 100 : 0
        pct_gt10 = (count > 0) ? (gt10 / count) * 100 : 0
        
        printf "%s\t%d\t%.2f\t%.2f\t%d\t%d\t%d\t%d\t%d\t%.2f\t%.2f\n", 
               sample, count, mean, median, max, zero, gt0, gt10, gt30, pct_gt0, pct_gt10
    }' ${file} >> coverage_summary.txt
done

# ==================== 步骤6: 生成可视化友好的摘要 ====================
echo ""
echo "[Step 6] Creating final summary..."

cd ${OUTPUT_DIR}

# 统计合并文件的信息
MERGED_LINES=$(wc -l < merged_chr12_coverage.txt)
MERGED_SAMPLES=$(head -1 merged_chr12_coverage.txt | awk '{print NF-2}')

echo ""
echo "============================================"
echo "Analysis Complete!"
echo "============================================"
echo "Region analyzed: ${CHR}:${START}-${CHR_LENGTH}"
echo "Region size: ${REGION_SIZE} bp"
echo "Total samples: ${MERGED_SAMPLES}"
echo "Total positions: $((MERGED_LINES - 1))"
echo ""
echo "Output files:"
echo "  1. ${OUTPUT_DIR}/merged_chr12_coverage.txt"
echo "     - Merged coverage data for all samples"
echo "     - Size: $(du -h merged_chr12_coverage.txt | cut -f1)"
echo ""
echo "  2. ${OUTPUT_DIR}/coverage_summary.txt"
echo "     - Statistical summary for each sample"
echo ""
echo "  3. ${OUTPUT_DIR}/temp/*.chr12.depth"
echo "     - Individual sample coverage files"
echo ""
echo "Quick preview:"
head -3 merged_chr12_coverage.txt | cut -f1-10

echo ""
echo "============================================"
echo "Usage examples:"
echo "============================================"
echo "# View merged file header and first few lines:"
echo "head -20 ${OUTPUT_DIR}/merged_chr12_coverage.txt"
echo ""
echo "# View summary statistics:"
echo "column -t ${OUTPUT_DIR}/coverage_summary.txt | less -S"
echo ""
echo "# Compress merged file to save space:"
echo "gzip ${OUTPUT_DIR}/merged_chr12_coverage.txt"
echo ""
echo "# Extract specific sample from merged file:"
echo "cut -f1,2,3 ${OUTPUT_DIR}/merged_chr12_coverage.txt | head"
echo "============================================"
