#!/bin/bash
set -e  # 遇到错误立即停止
set -u  # 使用未定义变量时报错
set -o pipefail  # 管道命令中任一失败则失败

# ==============================================================================
# 🛠️ 用户配置区域 (请修改这里!)
# ==============================================================================

# 1. 输入文件 (请使用绝对路径)
REF_FASTA="/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/05.allhic/OV53_1.primary.fa"
HIC_R1="/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/06.puzzle-hic/OV53_1-hic_R1.fastq.gz"
HIC_R2="/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/06.puzzle-hic/OV53_1-hic_R2.fastq.gz"

# 2. 项目参数
PROJECT_NAME="OV53-1"                         # 项目/物种名称 (不要包含空格)
CHROM_NUM=12                                  # 染色体数量 (必须准确)
ENZYME="MboI"                                 # 限制性内切酶 (如 DpnII, HindIII, MboI)
THREADS=64                                    # 使用的线程数

# 3. 软件路径
PUZZLE_DIR="/share/org/YZWL/yzwl_lixg/software/puzzle-hic"
JUICER_DIR="/share/org/YZWL/yzwl_lixg/software/juicer"
PUZZLE_CONDA_ENV="/share/org/YZWL/yzwl_lixg/miniforge3/envs/puzzle-hi-c"

# 4. 可选参数
PUZZLE_BINSIZE=10000                          # Puzzle binsize (默认10k)
PUZZLE_CUTOFF=0.35                            # Puzzle cutoff (默认0.35)
PUZZLE_INIT_TRIANGLE=6                        # Puzzle init triangle (默认6)

# 5. 断点续传控制
SKIP_JUICER=false                             # 跳过 Juicer（如果 merged_nodups.txt 已存在）
SKIP_PUZZLE=false                             # 跳过 Puzzle Hi-C

# ==============================================================================
# 🔧 辅助函数
# ==============================================================================

# 日志函数
log_info() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ℹ️  $*" | tee -a ${LOG_FILE}
}

log_success() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ✅ $*" | tee -a ${LOG_FILE}
}

log_warning() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ⚠️  $*" | tee -a ${LOG_FILE}
}

log_error() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ❌ $*" | tee -a ${LOG_FILE}
}

# 错误处理函数
error_exit() {
    log_error "$1"
    log_error "流程失败，请查看日志: ${LOG_FILE}"
    exit 1
}

# 检查文件是否存在
check_file() {
    if [ ! -f "$1" ]; then
        error_exit "文件不存在: $1"
    fi
}

# 检查目录是否存在
check_dir() {
    if [ ! -d "$1" ]; then
        error_exit "目录不存在: $1"
    fi
}

# ==============================================================================
# 🏁 流程开始
# ==============================================================================

echo "╔══════════════════════════════════════════════════════════════════════╗"
echo "║         🧬 Puzzle Hi-C 基因组组装流程 v3.0                          ║"
echo "╚══════════════════════════════════════════════════════════════════════╝"
echo ""

# 初始化工作目录
# WORK_DIR=$(pwd)/${PROJECT_NAME}_analysis
WORK_DIR=$(pwd)/puzzle_hic_output
mkdir -p ${WORK_DIR}
LOG_FILE="${WORK_DIR}/pipeline_$(date +%Y%m%d_%H%M%S).log"

log_info "Pipeline started"
log_info "Project: ${PROJECT_NAME}"
log_info "Working directory: ${WORK_DIR}"
log_info "Log file: ${LOG_FILE}"

# ==============================================================================
# Step 0: 环境检查
# ==============================================================================

log_info "[Step 0] 检查环境和输入文件..."

# 检查输入文件
check_file "${REF_FASTA}"
check_file "${HIC_R1}"
check_file "${HIC_R2}"
log_success "输入文件检查通过"

# 检查软件目录
check_dir "${PUZZLE_DIR}"
check_dir "${JUICER_DIR}"
check_dir "${PUZZLE_CONDA_ENV}"
log_success "软件目录检查通过"

# 设置软件路径
JUICER_SH="${JUICER_DIR}/scripts/juicer.sh"
JUICER_TOOLS="${JUICER_DIR}/scripts/juicer_tools"

if [ ! -f "${JUICER_SH}" ]; then
    error_exit "Juicer脚本未找到: ${JUICER_SH}"
fi

# 创建目录结构
mkdir -p ${WORK_DIR}/{references,restriction_sites,fastq,logs,backup}
log_success "目录结构创建完成"

# ==============================================================================
# Step 1: 准备参考基因组和必需文件（参照您的成功脚本）
# ==============================================================================

log_info "[Step 1] 准备参考基因组和必需文件..."
cd ${WORK_DIR}

# 1.1 链接基因组文件（Juicer 需要特定文件名）
GENOME_FA="${WORK_DIR}/genome.fa"
if [ ! -f "${GENOME_FA}" ]; then
    ln -s $(readlink -f ${REF_FASTA}) ${GENOME_FA}
    log_success "基因组文件链接创建: ${GENOME_FA}"
fi

# 1.2 构建 BWA 索引
if [ ! -f "${GENOME_FA}.bwt" ]; then
    log_info "构建 BWA 索引..."
    bwa index ${GENOME_FA} 2>&1 | tee -a ${LOG_FILE}
    log_success "BWA 索引构建完成"
else
    log_warning "BWA 索引已存在，跳过"
fi

# 1.3 生成 chrom.sizes 文件
CHROM_SIZES="${WORK_DIR}/references/${PROJECT_NAME}.chrom.sizes"
if [ ! -f "${CHROM_SIZES}" ]; then
    log_info "生成 chrom.sizes 文件..."
    samtools faidx ${GENOME_FA}
    awk -v OFS='\t' '{print $1, $2}' ${GENOME_FA}.fai > ${CHROM_SIZES}
    
    NUM_CHROMS=$(wc -l < ${CHROM_SIZES})
    log_info "检测到 ${NUM_CHROMS} 个序列"
    log_success "chrom.sizes 生成完成"
else
    log_warning "chrom.sizes 已存在，跳过"
fi

# 1.4 生成酶切位点文件（参照您的成功方法）
SITE_FILE="${WORK_DIR}/restriction_sites/${PROJECT_NAME}_${ENZYME}.txt"
SITE_FILE_TEMP="${WORK_DIR}/${PROJECT_NAME}_${ENZYME}.txt"

if [ -f "${SITE_FILE}" ] && [ -s "${SITE_FILE}" ]; then
    log_warning "酶切位点文件已存在且非空，跳过"
else
    log_info "生成限制性内切酶位点文件 (${ENZYME})..."
    
    GEN_SITE_SCRIPT="${JUICER_DIR}/misc/generate_site_positions.py"
    
    if [ ! -f "${GEN_SITE_SCRIPT}" ]; then
        log_warning "generate_site_positions.py 未找到，跳过位点文件生成"
        SITE_FILE=""
    else
        # 关键：先在当前目录生成，然后移动
        python2 ${GEN_SITE_SCRIPT} ${ENZYME} ${PROJECT_NAME} ${GENOME_FA}
        
        if [ ! -s "${SITE_FILE_TEMP}" ]; then
            log_error "generate_site_positions.py 未能正确生成酶切位点文件"
            SITE_FILE=""
        else
            # 排序并移动到目标位置
            log_info "排序并移动酶切位点文件..."
            sort -k1,1V -k2,2n "${SITE_FILE_TEMP}" > "${SITE_FILE}"
            rm -f "${SITE_FILE_TEMP}"
            log_success "酶切位点文件生成完成: ${SITE_FILE}"
        fi
    fi
fi

# ==============================================================================
# Step 2: 准备 Hi-C 数据
# ==============================================================================

log_info "[Step 2] 准备 Hi-C 测序数据..."
cd ${WORK_DIR}/fastq

# 创建符合 Juicer 要求的文件链接
ln -sf $(readlink -f ${HIC_R1}) ${PROJECT_NAME}_R1.fastq.gz
ln -sf $(readlink -f ${HIC_R2}) ${PROJECT_NAME}_R2.fastq.gz

log_success "Hi-C 数据链接创建完成"

# ==============================================================================
# Step 3: 运行 Juicer（使用您成功的参数格式）
# ==============================================================================

if [ "${SKIP_JUICER}" = "false" ]; then
    log_info "[Step 3] 运行 Juicer 流程..."
    cd ${WORK_DIR}
    
    # 检查是否有旧的输出目录
    MERGED_NODUPS="${WORK_DIR}/aligned/merged_nodups.txt"
    
    if [ -f "${MERGED_NODUPS}" ] && [ -s "${MERGED_NODUPS}" ]; then
        log_warning "检测到 merged_nodups.txt 已存在"
        read -p "是否跳过 Juicer 重新运行? (y/n): " -n 1 -r
        echo
        if [[ $REPLY =~ ^[Yy]$ ]]; then
            log_info "跳过 Juicer，使用现有文件"
            SKIP_JUICER=true
        else
            log_warning "清理旧的 aligned 和 splits 目录..."
            rm -rf aligned splits
        fi
    fi
    
    if [ "${SKIP_JUICER}" = "false" ]; then
        log_info "启动 Juicer 主流程..."
        log_warning "注意: Juicer 运行时间可能很长（数小时到数天），请耐心等待..."
        
        # 使用您成功脚本的参数格式
        JUICER_CMD="bash ${JUICER_SH} \
            -g ${PROJECT_NAME} \
            -d ${WORK_DIR} \
            -s ${ENZYME} \
            -p ${CHROM_SIZES} \
            -z ${GENOME_FA} \
            -D ${JUICER_DIR} \
            -t ${THREADS}"
        
        # 如果有酶切位点文件，添加 -y 参数
        if [ ! -z "${SITE_FILE}" ] && [ -f "${SITE_FILE}" ]; then
            JUICER_CMD="${JUICER_CMD} -y ${SITE_FILE}"
            log_info "使用酶切位点文件: ${SITE_FILE}"
        fi
        
        log_info "执行命令: ${JUICER_CMD}"
        
        # 运行 Juicer
        eval ${JUICER_CMD} 2>&1 | tee -a ${LOG_FILE}
        
        JUICER_EXIT_CODE=${PIPESTATUS[0]}
        if [ ${JUICER_EXIT_CODE} -ne 0 ]; then
            error_exit "Juicer 运行失败，退出码: ${JUICER_EXIT_CODE}"
        fi
        
        # 检查输出
        if [ ! -f "${MERGED_NODUPS}" ] || [ ! -s "${MERGED_NODUPS}" ]; then
            error_exit "Juicer 未生成有效的 merged_nodups.txt"
        fi
        
        VALID_PAIRS=$(wc -l < ${MERGED_NODUPS})
        log_success "Juicer 完成，有效 reads 对数: ${VALID_PAIRS}"
    fi
else
    log_warning "[Step 3] 跳过 Juicer (SKIP_JUICER=true)"
    MERGED_NODUPS="${WORK_DIR}/aligned/merged_nodups.txt"
    check_file "${MERGED_NODUPS}"
    VALID_PAIRS=$(wc -l < ${MERGED_NODUPS})
    log_info "使用现有 merged_nodups.txt，有效 reads 对数: ${VALID_PAIRS}"
fi

# ==============================================================================
# Step 4: 运行 Puzzle Hi-C
# ==============================================================================

if [ "${SKIP_PUZZLE}" = "false" ]; then
    log_info "[Step 4] 运行 Puzzle Hi-C (染色体级别组装)..."

    # 激活 conda 环境
    log_info "激活 Puzzle Hi-C conda 环境..."
    eval "$(conda shell.bash hook)"
    conda activate ${PUZZLE_CONDA_ENV}

    if [ $? -ne 0 ]; then
        error_exit "无法激活 conda 环境: ${PUZZLE_CONDA_ENV}"
    fi
    log_success "Conda 环境激活成功"

    # 创建输出目录
    PUZZLE_OUTPUT="${WORK_DIR}/Puzzle_Output"
    mkdir -p ${PUZZLE_OUTPUT}
    cd ${PUZZLE_OUTPUT}

    # 链接必要文件
    ln -sf ${MERGED_NODUPS} ./merged_nodups.txt

    # 运行 Puzzle Hi-C
    log_info "运行 Puzzle Hi-C main.py..."
    log_info "参数: 染色体数=${CHROM_NUM}, binsize=${PUZZLE_BINSIZE}, cutoff=${PUZZLE_CUTOFF}"

    python3 ${PUZZLE_DIR}/main.py \
        -c ${CHROM_NUM} \
        -p ${PROJECT_NAME} \
        -s ${PUZZLE_BINSIZE} \
        -t ${PUZZLE_CUTOFF} \
        -i ${PUZZLE_INIT_TRIANGLE} \
        -m merged_nodups.txt \
        -f ${GENOME_FA} \
        -j ${JUICER_TOOLS} \
        -n ${THREADS} 2>&1 | tee -a ${LOG_FILE}

    if [ $? -ne 0 ]; then
        conda deactivate
        error_exit "Puzzle Hi-C 运行失败"
    fi

    # 反激活 conda 环境
    conda deactivate
    log_success "Puzzle Hi-C 运行完成"
else
    log_warning "[Step 4] 跳过 Puzzle Hi-C (SKIP_PUZZLE=true)"
    PUZZLE_OUTPUT="${WORK_DIR}/Puzzle_Output"
fi

# ==============================================================================
# Step 5: 结果验证与汇总
# ==============================================================================

log_info "[Step 5] 验证结果文件..."

# 查找输出文件
FINAL_FASTA=$(find ${PUZZLE_OUTPUT} -name "*.fasta" -o -name "*.fa" 2>/dev/null | head -n 1)
FINAL_AGP=$(find ${PUZZLE_OUTPUT} -name "*.agp" 2>/dev/null | head -n 1)

if [ -f "${FINAL_FASTA}" ]; then
    log_success "找到组装结果: ${FINAL_FASTA}"
    
    # 统计scaffold信息
    NUM_SCAFFOLDS=$(grep -c "^>" ${FINAL_FASTA})
    TOTAL_LENGTH=$(awk '/^>/ {next} {len+=length($0)} END {print len}' ${FINAL_FASTA})
    
    log_info "Scaffold 数量: ${NUM_SCAFFOLDS}"
    log_info "总长度: ${TOTAL_LENGTH} bp"
else
    log_warning "未找到 FASTA 输出文件"
    NUM_SCAFFOLDS="N/A"
    TOTAL_LENGTH="N/A"
fi

if [ -f "${FINAL_AGP}" ]; then
    log_success "找到 AGP 文件: ${FINAL_AGP}"
else
    log_warning "未找到 AGP 文件"
fi

# ==============================================================================
# 流程完成
# ==============================================================================

echo ""
echo "╔══════════════════════════════════════════════════════════════════════╗"
echo "║                    🎉 流程成功完成！                                 ║"
echo "╚══════════════════════════════════════════════════════════════════════╝"
echo ""

log_success "Pipeline completed successfully"
log_info "结果目录: ${WORK_DIR}"
log_info "  - Juicer 输出: ${WORK_DIR}/aligned/"
log_info "  - Puzzle 输出: ${PUZZLE_OUTPUT}/"
log_info "  - 日志文件: ${LOG_FILE}"

# 生成结果摘要文件
SUMMARY_FILE="${WORK_DIR}/RESULTS_SUMMARY.txt"
cat > ${SUMMARY_FILE} << EOF
================================================================================
Puzzle Hi-C 流程运行摘要
================================================================================
运行时间: $(date)
项目名称: ${PROJECT_NAME}
工作目录: ${WORK_DIR}

输入文件:
  - 参考基因组: ${REF_FASTA}
  - Hi-C R1: ${HIC_R1}
  - Hi-C R2: ${HIC_R2}

参数设置:
  - 染色体数量: ${CHROM_NUM}
  - 限制性内切酶: ${ENZYME}
  - 线程数: ${THREADS}
  - Puzzle binsize: ${PUZZLE_BINSIZE}
  - Puzzle cutoff: ${PUZZLE_CUTOFF}

结果文件:
  - 组装结果: ${FINAL_FASTA}
  - AGP 文件: ${FINAL_AGP}
  - 日志文件: ${LOG_FILE}

统计信息:
  - Scaffold 数量: ${NUM_SCAFFOLDS}
  - 总长度: ${TOTAL_LENGTH} bp
  - 有效 reads 对数: ${VALID_PAIRS}

================================================================================
EOF

log_info "结果摘要已保存到: ${SUMMARY_FILE}"
cat ${SUMMARY_FILE}

echo ""
echo "✨ 所有任务完成！请查看上述路径中的结果文件。"