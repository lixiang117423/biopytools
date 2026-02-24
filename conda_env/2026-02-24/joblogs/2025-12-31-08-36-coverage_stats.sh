#!/bin/bash

# 完整的Chr12覆盖度提取和合并脚本
# 使用bedtools genomecov方法

# ==================== 设置参数 ====================
BAM_DIR="/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/31.重测序数据比对到挂载的结果上/02.mapping/bam"
CHR="Chr12"
START=125000000
OUTPUT_DIR="./chr12_coverage"
THREADS=4  # 可根据服务器资源调整

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
    
    # 使用bedtools计算覆盖度（每个碱基位点）
    # -d 参数输出每个位置的覆盖度
    # -ibam 指定输入BAM文件
    bedtools genomecov -ibam ${bam} -d | \
    awk -v chr="${CHR}" -v start=${START} '$1==chr && $2>=start {print $1"\t"$2"\t"$3}' \
    > ${OUTPUT_DIR}/temp/${sample}.chr12.depth
    
    # 检查输出文件是否为空
    if [ ! -s ${OUTPUT_DIR}/temp/${sample}.chr12.depth ]; then
        echo "    Warning: No data extracted for ${sample}"
    else
        LINE_COUNT=$(wc -l < ${OUTPUT_DIR}/temp/${sample}.chr12.depth)
        echo "    Extracted ${LINE_COUNT} positions"
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

# ==================== 步骤4: 合并所有样本的覆盖度数据 ====================
echo ""
echo "[Step 4] Merging all samples..."

cd ${OUTPUT_DIR}/temp

# 创建表头
echo -n "Chr"$'\t'"Pos" > ../merged_chr12_coverage.txt
for file in *.chr12.depth; do
    sample=$(basename ${file} .chr12.depth)
    echo -n $'\t'"${sample}" >> ../merged_chr12_coverage.txt
done
echo "" >> ../merged_chr12_coverage.txt

# 使用awk合并所有文件
echo "  Merging data with awk..."

awk '
BEGIN {
    OFS="\t"
}
{
    # $1=Chr, $2=Pos, $3=Depth
    key = $1":"$2
    depth[key][FILENAME] = $3
    
    # 记录位置信息（只需要记录一次）
    if (!(key in positions)) {
        positions[key] = $1"\t"$2
        pos_order[++order_count] = key
    }
    
    # 记录所有文件名
    if (!(FILENAME in files)) {
        file_order[++file_count] = FILENAME
        files[FILENAME] = 1
    }
}
END {
    # 按位置顺序输出（已经按染色体位置排序）
    for (i = 1; i <= order_count; i++) {
        pos = pos_order[i]
        printf "%s", positions[pos]
        
        # 按文件顺序输出每个样本的覆盖度
        for (j = 1; j <= file_count; j++) {
            file = file_order[j]
            if (pos in depth && file in depth[pos]) {
                printf "\t%s", depth[pos][file]
            } else {
                printf "\t0"
            }
        }
        printf "\n"
    }
}' *.chr12.depth | sort -t$'\t' -k2,2n >> ../merged_chr12_coverage.txt

# ==================== 步骤5: 生成统计摘要 ====================
echo ""
echo "[Step 5] Generating summary statistics..."

cd ${OUTPUT_DIR}

# 计算每个样本的基本统计信息
echo "Sample"$'\t'"Total_Positions"$'\t'"Mean_Coverage"$'\t'"Median_Coverage"$'\t'"Positions_0X"$'\t'"Positions_>10X"$'\t'"Positions_>30X" > coverage_summary.txt

for file in temp/*.chr12.depth; do
    sample=$(basename ${file} .chr12.depth)
    
    awk -v sample="${sample}" '
    {
        sum += $3
        count++
        depths[count] = $3
        if ($3 == 0) zero++
        if ($3 > 10) gt10++
        if ($3 > 30) gt30++
    }
    END {
        mean = (count > 0) ? sum/count : 0
        
        # 计算中位数
        asort(depths)
        if (count % 2 == 1) {
            median = depths[int(count/2) + 1]
        } else {
            median = (depths[count/2] + depths[count/2 + 1]) / 2
        }
        
        printf "%s\t%d\t%.2f\t%.2f\t%d\t%d\t%d\n", 
               sample, count, mean, median, zero, gt10, gt30
    }' ${file} >> coverage_summary.txt
done

# ==================== 步骤6: 清理和完成 ====================
echo ""
echo "[Step 6] Finalizing..."

# 压缩合并后的文件以节省空间
echo "  Compressing merged file..."
gzip -f merged_chr12_coverage.txt

# 计算最终文件大小
FINAL_SIZE=$(du -h merged_chr12_coverage.txt.gz | cut -f1)

echo ""
echo "============================================"
echo "Analysis Complete!"
echo "============================================"
echo "Output files:"
echo "  - Merged coverage: ${OUTPUT_DIR}/merged_chr12_coverage.txt.gz"
echo "  - Summary stats:   ${OUTPUT_DIR}/coverage_summary.txt"
echo "  - Individual files: ${OUTPUT_DIR}/temp/*.chr12.depth"
echo ""
echo "Merged file size: ${FINAL_SIZE}"
echo ""
echo "To view the merged file:"
echo "  zcat ${OUTPUT_DIR}/merged_chr12_coverage.txt.gz | head"
echo ""
echo "To keep temp files, don't run the cleanup."
echo "To remove temp files: rm -rf ${OUTPUT_DIR}/temp"
echo "============================================"
