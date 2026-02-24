#!/bin/bash

# 🧬 GATK联合分型脚本 - 大豆疫霉菌 (CombineGVCFs方法)
# 使用 CombineGVCFs + GenotypeGVCFs 流程
# 日期: $(date +%Y-%m-%d)

set -e
set -u
set -o pipefail

# ========================================
# 🐍 Python环境检查
# ========================================
if ! command -v python &> /dev/null; then
    if command -v python3 &> /dev/null; then
        TEMP_BIN="${HOME}/.local/bin"
        mkdir -p "${TEMP_BIN}"
        ln -sf "$(which python3)" "${TEMP_BIN}/python"
        export PATH="${TEMP_BIN}:${PATH}"
    fi
fi

start_time=$(date +%s)

# ========================================
# 📁 路径配置
# ========================================
PROJECT_DIR="/share/org/YZWL/yzwl_lixg/project/19.大豆疫霉菌"
GVCF_DIR="${PROJECT_DIR}/02.mapping/vcf"
GENOME="${PROJECT_DIR}/01.data/genome/Phytophthora_sojae_JS2.fa"
OUTPUT_DIR="${PROJECT_DIR}/03.gatk_joint"
LOG_DIR="${OUTPUT_DIR}/logs"
CHR_DIR="${OUTPUT_DIR}/by_chromosome"
TMP_DIR="${OUTPUT_DIR}/tmp"

# ========================================
# 🔧 资源配置
# ========================================
GATK="gatk"
MAX_THREADS=88
COMBINE_MEM="20G"    # CombineGVCFs内存
GENOTYPE_MEM="40G"   # GenotypeGVCFs内存
GATHER_MEM="900G"    # GatherVcfs内存
THREADS_PER_CHR=2    # 每个染色体任务的线程数
MAX_JOBS=$((MAX_THREADS / THREADS_PER_CHR))

echo "════════════════════════════════════════════════════════════"
echo "🧬 GATK联合分型 - CombineGVCFs方法"
echo "════════════════════════════════════════════════════════════"
echo "📅 开始时间: $(date '+%Y-%m-%d %H:%M:%S')"
echo "📂 输出目录: ${OUTPUT_DIR}"
echo "⚙️  最大并发任务: ${MAX_JOBS}"
echo "════════════════════════════════════════════════════════════"
echo ""

# 创建目录
mkdir -p "${OUTPUT_DIR}" "${LOG_DIR}" "${CHR_DIR}" "${TMP_DIR}"

# ========================================
# 📝 收集所有gVCF文件
# ========================================
echo "📋 正在收集gVCF文件..."

cd "${GVCF_DIR}" || exit 1
GVCF_FILES=(*.g.vcf.gz)
cd - > /dev/null

sample_count=${#GVCF_FILES[@]}

if [ ${sample_count} -eq 0 ]; then
    echo "❌ 未找到gVCF文件!"
    exit 1
fi

echo "✅ 找到 ${sample_count} 个样品"
echo "   示例: ${GVCF_FILES[0]}, ${GVCF_FILES[1]}, ..."
echo ""

# ========================================
# 🧬 提取染色体列表
# ========================================
echo "🧬 提取染色体列表..."
CHR_LIST="${OUTPUT_DIR}/chromosome_list.txt"

[ ! -f "${GENOME}.fai" ] && samtools faidx "${GENOME}"
awk '{print $1}' "${GENOME}.fai" > "${CHR_LIST}"
chr_count=$(wc -l < "${CHR_LIST}")

echo "✅ 找到 ${chr_count} 条染色体"
head -5 "${CHR_LIST}" | while read chr; do
    echo "   • ${chr}"
done
[ ${chr_count} -gt 5 ] && echo "   ... 还有 $((chr_count - 5)) 条"
echo ""

# ========================================
# 🔗 步骤1: CombineGVCFs - 按染色体合并
# ========================================
echo "════════════════════════════════════════════════════════════"
echo "🔗 步骤1: CombineGVCFs (并发: ${MAX_JOBS})"
echo "════════════════════════════════════════════════════════════"

while read chr; do
    while [ $(jobs -r | wc -l) -ge ${MAX_JOBS} ]; do sleep 2; done
    
    echo "🔄 [$(date '+%H:%M:%S')] ${chr} - CombineGVCFs"
    
    (
        COMBINED_GVCF="${CHR_DIR}/${chr}.combined.g.vcf.gz"
        CHR_TMP="${TMP_DIR}/tmp_${chr}_combine"
        mkdir -p "${CHR_TMP}"
        
        # 构建所有样品的-V参数
        VARIANT_ARGS=""
        for gvcf in "${GVCF_DIR}"/*.g.vcf.gz; do
            VARIANT_ARGS="${VARIANT_ARGS} -V ${gvcf}"
        done
        
        ${GATK} --java-options "-Xmx${COMBINE_MEM} -Djava.io.tmpdir=${CHR_TMP}" \
            CombineGVCFs \
            -R "${GENOME}" \
            ${VARIANT_ARGS} \
            -L "${chr}" \
            -O "${COMBINED_GVCF}" \
            --tmp-dir "${CHR_TMP}" \
            > "${LOG_DIR}/combine_${chr}.log" 2>&1
        
        status=$?
        rm -rf "${CHR_TMP}"
        
        if [ ${status} -eq 0 ]; then
            echo "✅ [$(date '+%H:%M:%S')] ${chr} - CombineGVCFs完成"
        else
            echo "❌ [$(date '+%H:%M:%S')] ${chr} - CombineGVCFs失败"
            exit 1
        fi
    ) &
done < "${CHR_LIST}"

echo "⏳ 等待所有CombineGVCFs任务完成..."
wait
[ $? -ne 0 ] && echo "❌ 有任务失败!" && exit 1
echo "✅ CombineGVCFs全部完成!"
echo ""

# ========================================
# 🧬 步骤2: GenotypeGVCFs - 联合分型
# ========================================
echo "════════════════════════════════════════════════════════════"
echo "🧬 步骤2: GenotypeGVCFs (并发: ${MAX_JOBS})"
echo "════════════════════════════════════════════════════════════"

while read chr; do
    while [ $(jobs -r | wc -l) -ge ${MAX_JOBS} ]; do sleep 2; done
    
    echo "🔄 [$(date '+%H:%M:%S')] ${chr} - GenotypeGVCFs"
    
    (
        COMBINED_GVCF="${CHR_DIR}/${chr}.combined.g.vcf.gz"
        OUTPUT_VCF="${CHR_DIR}/${chr}.vcf.gz"
        CHR_TMP="${TMP_DIR}/tmp_${chr}_genotype"
        mkdir -p "${CHR_TMP}"
        
        ${GATK} --java-options "-Xmx${GENOTYPE_MEM} -Djava.io.tmpdir=${CHR_TMP}" \
            GenotypeGVCFs \
            -R "${GENOME}" \
            -V "${COMBINED_GVCF}" \
            -L "${chr}" \
            -O "${OUTPUT_VCF}" \
            --tmp-dir "${CHR_TMP}" \
            > "${LOG_DIR}/genotype_${chr}.log" 2>&1
        
        status=$?
        rm -rf "${CHR_TMP}"
        
        if [ ${status} -eq 0 ]; then
            echo "✅ [$(date '+%H:%M:%S')] ${chr} - GenotypeGVCFs完成"
            # 删除中间的combined文件以节省空间
            rm -f "${COMBINED_GVCF}" "${COMBINED_GVCF}.tbi"
        else
            echo "❌ [$(date '+%H:%M:%S')] ${chr} - GenotypeGVCFs失败"
            exit 1
        fi
    ) &
done < "${CHR_LIST}"

echo "⏳ 等待所有GenotypeGVCFs任务完成..."
wait
[ $? -ne 0 ] && echo "❌ 有任务失败!" && exit 1
echo "✅ GenotypeGVCFs全部完成!"
echo ""

# ========================================
# 🔗 步骤3: 合并所有染色体VCF
# ========================================
echo "════════════════════════════════════════════════════════════"
echo "🔗 步骤3: 合并染色体VCF"
echo "════════════════════════════════════════════════════════════"

VCF_LIST="${OUTPUT_DIR}/vcf_list.txt"
> "${VCF_LIST}"

while read chr; do
    vcf_file="${CHR_DIR}/${chr}.vcf.gz"
    if [ -f "${vcf_file}" ]; then
        echo "${vcf_file}" >> "${VCF_LIST}"
    else
        echo "⚠️  警告: ${chr}.vcf.gz 不存在"
    fi
done < "${CHR_LIST}"

FINAL_VCF="${OUTPUT_DIR}/joint_genotyped.vcf.gz"
GATHER_TMP="${TMP_DIR}/tmp_gather"
mkdir -p "${GATHER_TMP}"

echo "🔗 正在合并 $(wc -l < ${VCF_LIST}) 个VCF文件..."

${GATK} --java-options "-Xmx${GATHER_MEM} -Djava.io.tmpdir=${GATHER_TMP}" \
    GatherVcfs \
    -I "${VCF_LIST}" \
    -O "${FINAL_VCF}" \
    > "${LOG_DIR}/gather.log" 2>&1

status=$?
rm -rf "${GATHER_TMP}"

if [ ${status} -eq 0 ]; then
    echo "✅ VCF合并完成!"
else
    echo "❌ VCF合并失败，查看日志: ${LOG_DIR}/gather.log"
    exit 1
fi
echo ""

# ========================================
# 📊 统计信息
# ========================================
echo "════════════════════════════════════════════════════════════"
echo "📊 统计信息"
echo "════════════════════════════════════════════════════════════"

if [ -f "${FINAL_VCF}" ]; then
    total_vars=$(zgrep -vc "^#" "${FINAL_VCF}" || echo 0)
    snps=$(zgrep -v "^#" "${FINAL_VCF}" | awk '{if(length($4)==1 && length($5)==1) print}' | wc -l)
    indels=$(zgrep -v "^#" "${FINAL_VCF}" | awk '{if(length($4)!=length($5)) print}' | wc -l)
else
    total_vars=0
    snps=0
    indels=0
fi

end_time=$(date +%s)
runtime=$((end_time - start_time))
hours=$((runtime / 3600))
minutes=$(((runtime % 3600) / 60))
seconds=$((runtime % 60))

echo "✅ 完成时间: $(date '+%Y-%m-%d %H:%M:%S')"
echo "⏱️  运行时间: ${hours}h ${minutes}m ${seconds}s"
echo ""
echo "📊 样品数: ${sample_count}"
echo "📊 染色体数: ${chr_count}"
echo "📊 总变异: ${total_vars}"
echo "📊 SNPs: ${snps}"
echo "📊 InDels: ${indels}"
echo ""
echo "📁 最终VCF: ${FINAL_VCF}"
echo "📂 日志目录: ${LOG_DIR}"
echo "📂 染色体VCF: ${CHR_DIR}/"
echo ""
echo "💡 后续操作:"
echo "   1. 变异过滤: gatk VariantFiltration -R ${GENOME} -V ${FINAL_VCF} ..."
echo "   2. 质量统计: bcftools stats ${FINAL_VCF} > stats.txt"
echo "   3. 节省空间: rm -rf ${CHR_DIR} ${TMP_DIR}"
echo "════════════════════════════════════════════════════════════"