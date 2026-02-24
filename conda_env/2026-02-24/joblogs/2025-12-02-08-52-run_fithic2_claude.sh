#!/bin/bash

# ==================== 配置参数 ====================
HIC_FILE="/share/org/YZWL/yzwl_lixg/project/22.Est1_hic/03.hic/fithic2/inter_30.hic"
OUT_DIR="/share/org/YZWL/yzwl_lixg/project/22.Est1_hic/03.hic/fithic2/result_5kb"

FITHIC_BIN="/share/org/YZWL/yzwl_lixg/miniforge3/envs/fithic-v.2.0.8/bin/fithic"
UTILS_DIR="/share/org/YZWL/yzwl_lixg/software/fithic2/utils"
JUICER_JAR="/share/org/YZWL/yzwl_lixg/software/juicer/scripts/juicer_tools_1.22.01.jar"

RESOLUTION=5000
LOWER_BOUND=2000
UPPER_BOUND=3000000
QVALUE_THRESH=0.01
CHROMS=("Chr1" "Chr2" "Chr3" "Chr4" "Chr5")

# 染色体长度文件（需要创建）
CHR_LENS="${OUT_DIR}/chr_lengths.txt"

# ==================== 创建输出目录 ====================
mkdir -p ${OUT_DIR}
mkdir -p ${OUT_DIR}/temp
mkdir -p ${OUT_DIR}/contacts

echo "=== FitHiC Pipeline Started ==="
echo "Resolution: ${RESOLUTION}"
echo "Output Directory: ${OUT_DIR}"

# ==================== 步骤1: 创建染色体长度文件 ====================
echo "Step 1: Creating chromosome lengths file..."
cat > ${CHR_LENS} << EOF
Chr1	29412698
Chr2	19420461
Chr3	23119851
Chr4	18582608
Chr5	26592585
EOF

echo "Chromosome lengths file created at: ${CHR_LENS}"

# ==================== 步骤2: 生成Fragments文件 ====================
echo "Step 2: Creating fragments file..."
FRAGMENTS_FILE="${OUT_DIR}/fragments.txt"

python3 ${UTILS_DIR}/createFitHiCFragments-fixedsize.py \
    --chrLens ${CHR_LENS} \
    --resolution ${RESOLUTION} \
    --outFile ${FRAGMENTS_FILE}

# 压缩fragments文件
gzip -f ${FRAGMENTS_FILE}
echo "Fragments file created: ${FRAGMENTS_FILE}.gz"

# ==================== 步骤3: 从.hic文件提取contacts ====================
echo "Step 3: Extracting contacts from .hic file..."
ALL_CONTACTS="${OUT_DIR}/contacts/all_contacts.txt"
> ${ALL_CONTACTS}  # 清空文件

for chr1 in "${CHROMS[@]}"; do
    for chr2 in "${CHROMS[@]}"; do
        # 只处理chr1 <= chr2的情况，避免重复
        if [[ "${chr1}" < "${chr2}" ]] || [[ "${chr1}" == "${chr2}" ]]; then
            echo "  Extracting ${chr1} - ${chr2}..."
            TEMP_CONTACT="${OUT_DIR}/temp/${chr1}_${chr2}.txt"
            
            # 使用Juicer dump命令提取接触矩阵
            java -jar ${JUICER_JAR} dump observed NONE ${HIC_FILE} ${chr1} ${chr2} BP ${RESOLUTION} ${TEMP_CONTACT}
            
            if [ -f "${TEMP_CONTACT}" ]; then
                # 转换为FitHiC格式：chr1 mid1 chr2 mid2 count
                awk -v chr1="${chr1}" -v chr2="${chr2}" '{print chr1"\t"$1"\t"chr2"\t"$2"\t"$3}' ${TEMP_CONTACT} >> ${ALL_CONTACTS}
                rm ${TEMP_CONTACT}
            fi
        fi
    done
done

# 压缩contacts文件
gzip -f ${ALL_CONTACTS}
echo "Contacts file created: ${ALL_CONTACTS}.gz"

# ==================== 步骤4: 运行FitHiC ====================
echo "Step 4: Running FitHiC..."

${FITHIC_BIN} \
    -f ${FRAGMENTS_FILE}.gz \
    -i ${ALL_CONTACTS}.gz \
    -o ${OUT_DIR} \
    -r ${RESOLUTION} \
    -L ${LOWER_BOUND} \
    -U ${UPPER_BOUND} \
    -l inter_30_${RESOLUTION}bp \
    -p 2 \
    -b 100 \
    -x intraOnly \
    -v

echo "FitHiC analysis completed!"

# ==================== 步骤5: 过滤显著互作 ====================
echo "Step 5: Filtering significant interactions..."
FITHIC_OUTPUT="${OUT_DIR}/inter_30_${RESOLUTION}bp.spline_pass2.significances.txt.gz"
FILTERED_OUTPUT="${OUT_DIR}/inter_30_${RESOLUTION}bp.significant_q${QVALUE_THRESH}.txt"

if [ -f "${FITHIC_OUTPUT}" ]; then
    zcat ${FITHIC_OUTPUT} | awk -v q="${QVALUE_THRESH}" 'NR==1 || $7<=q' > ${FILTERED_OUTPUT}
    
    TOTAL_INTERACTIONS=$(zcat ${FITHIC_OUTPUT} | tail -n +2 | wc -l)
    SIGNIFICANT_INTERACTIONS=$(tail -n +2 ${FILTERED_OUTPUT} | wc -l)
    
    echo "Total interactions: ${TOTAL_INTERACTIONS}"
    echo "Significant interactions (q<=${QVALUE_THRESH}): ${SIGNIFICANT_INTERACTIONS}"
    echo "Filtered results saved to: ${FILTERED_OUTPUT}"
else
    echo "Warning: FitHiC output file not found!"
fi

# ==================== 步骤6: 生成可视化文件 (可选) ====================
echo "Step 6: Creating UCSC visualization file..."
if [ -f "${FILTERED_OUTPUT}" ]; then
    UCSC_OUTPUT="${OUT_DIR}/inter_30_${RESOLUTION}bp.ucsc_interact.txt"
    bash ${UTILS_DIR}/visualize-UCSC.sh ${FILTERED_OUTPUT} ${UCSC_OUTPUT} ${QVALUE_THRESH}
    echo "UCSC visualization file created: ${UCSC_OUTPUT}"
fi

# ==================== 清理临时文件 ====================
echo "Step 7: Cleaning up temporary files..."
rm -rf ${OUT_DIR}/temp

echo "=== FitHiC Pipeline Completed Successfully ==="
echo "Results are in: ${OUT_DIR}"