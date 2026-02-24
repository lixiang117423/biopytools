#!/bin/bash

# 三代转录组比对到基因组脚本
# Author: Generated Script
# Date: 2025-11-16

set -e  # 遇到错误立即退出

# ================================
# 📁 路径配置
# ================================
BASE_DIR="/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/66.三代转录组比对到组装的基因组"
PACBIO_DIR="${BASE_DIR}/pacbio"
OUTPUT_DIR="${BASE_DIR}"
REF_GENOME="${BASE_DIR}/Orychophragmus_violaceus_OV53_1_HiFi.fa"

# 输入BAM文件
JINGYEA_BAM="${PACBIO_DIR}/jingyeA.sreads.bam"
GENA_BAM="${PACBIO_DIR}/genA.sreads.bam"

# 输出目录
RESULT_DIR="${OUTPUT_DIR}/alignment_results"
LOG_DIR="${RESULT_DIR}/logs"
STATS_DIR="${RESULT_DIR}/stats"

# 创建输出目录
mkdir -p ${RESULT_DIR} ${LOG_DIR} ${STATS_DIR}

# 日志文件
LOG_FILE="${LOG_DIR}/alignment_$(date +%Y%m%d_%H%M%S).log"

# ================================
# 📝 日志函数
# ================================
log_info() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ℹ️  INFO: $1" | tee -a ${LOG_FILE}
}

log_success() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ✅ SUCCESS: $1" | tee -a ${LOG_FILE}
}

log_error() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ❌ ERROR: $1" | tee -a ${LOG_FILE}
}

log_warning() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ⚠️  WARNING: $1" | tee -a ${LOG_FILE}
}

log_step() {
    echo "" | tee -a ${LOG_FILE}
    echo "==========================================" | tee -a ${LOG_FILE}
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] 🚀 $1" | tee -a ${LOG_FILE}
    echo "==========================================" | tee -a ${LOG_FILE}
}

# ================================
# 🔍 检查输入文件
# ================================
log_step "检查输入文件和依赖"

check_file() {
    if [ ! -f "$1" ]; then
        log_error "文件不存在: $1"
        exit 1
    else
        log_success "文件存在: $1"
    fi
}

check_file ${REF_GENOME}
check_file ${JINGYEA_BAM}
check_file ${GENA_BAM}

# 检查必要的软件
check_command() {
    if ! command -v $1 &> /dev/null; then
        log_error "未找到命令: $1，请先安装"
        exit 1
    else
        log_success "命令可用: $1 ($(command -v $1))"
    fi
}

check_command minimap2
check_command samtools

# ================================
# 🧬 处理函数
# ================================
process_sample() {
    local sample_name=$1
    local input_bam=$2
    
    log_step "处理样本: ${sample_name}"
    
    # 输出文件路径
    local fastq="${RESULT_DIR}/${sample_name}.fastq"
    local aligned_sam="${RESULT_DIR}/${sample_name}_aligned.sam"
    local aligned_bam="${RESULT_DIR}/${sample_name}_aligned.bam"
    local sorted_bam="${RESULT_DIR}/${sample_name}_aligned.sorted.bam"
    local stats_file="${STATS_DIR}/${sample_name}_stats.txt"
    
    # Step 1: 将BAM转换为FASTQ
    log_info "步骤1️⃣: 将BAM转换为FASTQ..."
    samtools fastq ${input_bam} > ${fastq} 2>> ${LOG_FILE}
    log_success "FASTQ文件生成: ${fastq}"
    
    # 统计reads数量
    local read_count=$(grep -c "^@" ${fastq} || true)
    log_info "📊 总reads数: ${read_count}"
    
    # Step 2: 使用minimap2比对到参考基因组
    log_info "步骤2️⃣: 使用minimap2比对到参考基因组（使用88线程）..."
    minimap2 -ax splice -uf -k14 \
        --secondary=no \
        -t 88 \
        ${REF_GENOME} \
        ${fastq} > ${aligned_sam} 2>> ${LOG_FILE}
    log_success "比对完成: ${aligned_sam}"
    
    # Step 3: 转换为BAM格式
    log_info "步骤3️⃣: 转换为BAM格式..."
    samtools view -bS ${aligned_sam} > ${aligned_bam} 2>> ${LOG_FILE}
    log_success "BAM文件生成: ${aligned_bam}"
    
    # Step 4: 排序BAM文件
    log_info "步骤4️⃣: 排序BAM文件（使用88线程）..."
    samtools sort -@ 88 -o ${sorted_bam} ${aligned_bam} 2>> ${LOG_FILE}
    log_success "排序完成: ${sorted_bam}"
    
    # Step 5: 建立索引
    log_info "步骤5️⃣: 建立索引..."
    samtools index ${sorted_bam} 2>> ${LOG_FILE}
    log_success "索引创建完成"
    
    # Step 6: 生成比对统计
    log_info "步骤6️⃣: 生成比对统计..."
    samtools flagstat ${sorted_bam} > ${stats_file} 2>> ${LOG_FILE}
    
    # 解析统计结果
    echo "" | tee -a ${LOG_FILE}
    echo "📈 ${sample_name} 比对统计结果:" | tee -a ${LOG_FILE}
    echo "----------------------------------------" | tee -a ${LOG_FILE}
    cat ${stats_file} | tee -a ${LOG_FILE}
    echo "----------------------------------------" | tee -a ${LOG_FILE}
    
    # 计算比对率
    local total_reads=$(grep "in total" ${stats_file} | awk '{print $1}')
    local mapped_reads=$(grep "mapped (" ${stats_file} | head -1 | awk '{print $1}')
    local mapping_rate=$(grep "mapped (" ${stats_file} | head -1 | awk '{print $5}' | tr -d '()')
    
    log_success "总reads数: ${total_reads}"
    log_success "比对上的reads: ${mapped_reads}"
    log_success "比对率: ${mapping_rate}"
    
    # 清理中间文件
    log_info "步骤7️⃣: 清理中间文件..."
    rm -f ${fastq} ${aligned_sam} ${aligned_bam}
    log_success "中间文件已清理"
    
    echo "" | tee -a ${LOG_FILE}
}

# ================================
# 🏃 执行比对
# ================================
log_step "开始三代转录组比对分析"

# 处理 jingyeA 样本
process_sample "jingyeA" ${JINGYEA_BAM}

# 处理 genA 样本
process_sample "genA" ${GENA_BAM}

# ================================
# 📊 生成汇总报告
# ================================
log_step "生成汇总报告"

SUMMARY_FILE="${RESULT_DIR}/alignment_summary.txt"

cat > ${SUMMARY_FILE} << EOF
=========================================
🧬 三代转录组比对汇总报告
=========================================
分析日期: $(date '+%Y-%m-%d %H:%M:%S')
参考基因组: ${REF_GENOME}

-----------------------------------------
📂 样本信息:
-----------------------------------------
样本1: jingyeA
  - 输入: ${JINGYEA_BAM}
  - 输出: ${RESULT_DIR}/jingyeA_aligned.sorted.bam
  
样本2: genA
  - 输入: ${GENA_BAM}
  - 输出: ${RESULT_DIR}/genA_aligned.sorted.bam

-----------------------------------------
📈 比对统计:
-----------------------------------------

🔹 jingyeA样本:
$(cat ${STATS_DIR}/jingyeA_stats.txt)

🔹 genA样本:
$(cat ${STATS_DIR}/genA_stats.txt)

-----------------------------------------
📁 输出文件位置:
-----------------------------------------
结果目录: ${RESULT_DIR}
统计目录: ${STATS_DIR}
日志目录: ${LOG_DIR}

=========================================
EOF

cat ${SUMMARY_FILE} | tee -a ${LOG_FILE}

log_success "所有分析完成！"
log_info "详细日志: ${LOG_FILE}"
log_info "汇总报告: ${SUMMARY_FILE}"

echo ""
echo "🎉 分析流程全部完成！"
echo "📂 结果目录: ${RESULT_DIR}"
echo "📄 查看汇总报告: cat ${SUMMARY_FILE}"