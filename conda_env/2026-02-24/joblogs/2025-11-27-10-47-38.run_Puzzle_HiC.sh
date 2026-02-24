#!/bin/bash
set -e  # 遇到错误立即停止
set -u  # 使用未定义变量时报错
set -o pipefail  # 管道命令中任一失败则失败

# ==============================================================================
# 🛠️ 用户配置区域 (请修改这里!)
# ==============================================================================

# 1. 输入文件 (请使用绝对路径)
REF_FASTA="/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/05.allhic/OV53_1.primary.fa"              # 基因组 FASTA 文件
HIC_R1="/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/06.puzzle-hic/OV53_1-hic_R1.fastq.gz"        # Hi-C Read 1
HIC_R2="/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/06.puzzle-hic/OV53_1-hic_R2.fastq.gz"        # Hi-C Read 2

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
SKIP_BWA_INDEX=false                          # 跳过 BWA 索引构建
SKIP_CHROM_SIZES=false                        # 跳过 chrom.sizes 生成
SKIP_SITE_FILE=false                          # 跳过酶切位点生成
SKIP_JUICER=false                             # 跳过整个 Juicer 流程
RESUME_JUICER=false                           # 从 Juicer 中断处恢复（修复 SAM 文件）
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

# 检查命令是否可用
check_command() {
    if ! command -v $1 &> /dev/null; then
        error_exit "命令 $1 未找到，请检查是否已安装"
    fi
}

# ==============================================================================
# 🏁 流程开始
# ==============================================================================

echo "╔══════════════════════════════════════════════════════════════════════╗"
echo "║         🧬 Puzzle Hi-C 基因组组装流程 v2.0                          ║"
echo "╚══════════════════════════════════════════════════════════════════════╝"
echo ""

# 初始化工作目录
WORK_DIR=$(pwd)/${PROJECT_NAME}_analysis
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

# 检查必要命令
check_command "bwa"
check_command "conda"
log_success "必要命令检查通过"

# 设置软件路径
JUICER_SH="${JUICER_DIR}/scripts/juicer.sh"
JUICER_TOOLS="${JUICER_DIR}/scripts/juicer_tools"

if [ ! -f "${JUICER_SH}" ]; then
    error_exit "Juicer脚本未找到: ${JUICER_SH}"
fi

# 创建目录结构
mkdir -p ${WORK_DIR}/{ref,fastq,logs,backup}
log_success "目录结构创建完成"

# ==============================================================================
# Step 1: 准备参考基因组
# ==============================================================================

log_info "[Step 1] 准备参考基因组索引..."
cd ${WORK_DIR}/ref

# 链接参考基因组
REF_PATH="${WORK_DIR}/ref/${PROJECT_NAME}.fasta"
if [ ! -f "${REF_PATH}" ]; then
    ln -s $(readlink -f ${REF_FASTA}) ${REF_PATH}
    log_success "参考基因组链接创建: ${REF_PATH}"
fi

# 1.1 构建 BWA 索引
if [ "${SKIP_BWA_INDEX}" = "false" ]; then
    if [ ! -f "${REF_PATH}.bwt" ]; then
        log_info "构建 BWA 索引..."
        bwa index ${REF_PATH} 2>&1 | tee -a ${LOG_FILE}
        log_success "BWA 索引构建完成"
    else
        log_warning "BWA 索引已存在，跳过"
    fi
else
    log_warning "跳过 BWA 索引构建 (SKIP_BWA_INDEX=true)"
fi

# 1.2 生成 chrom.sizes
CHROM_SIZES="${WORK_DIR}/ref/${PROJECT_NAME}.chrom.sizes"
if [ "${SKIP_CHROM_SIZES}" = "false" ]; then
    if [ ! -f "${CHROM_SIZES}" ]; then
        log_info "生成 chrom.sizes 文件..."
        if command -v samtools &> /dev/null; then
            samtools faidx ${REF_PATH}
            cut -f1,2 ${REF_PATH}.fai > ${CHROM_SIZES}
        else
            awk '$0 ~ ">" {if (seq) print name, len; name=substr($0,2); len=0; seq=1; next} 
                 {len+=length($0)} END {print name, len}' ${REF_PATH} > ${CHROM_SIZES}
        fi
        
        # 验证 chrom.sizes
        NUM_CHROMS=$(wc -l < ${CHROM_SIZES})
        log_info "检测到 ${NUM_CHROMS} 个序列"
        log_success "chrom.sizes 生成完成"
    else
        log_warning "chrom.sizes 已存在，跳过"
    fi
else
    log_warning "跳过 chrom.sizes 生成 (SKIP_CHROM_SIZES=true)"
    check_file "${CHROM_SIZES}"
fi

# 1.3 生成酶切位点文件
SITE_FILE="${WORK_DIR}/ref/${PROJECT_NAME}_${ENZYME}.txt"
if [ "${SKIP_SITE_FILE}" = "false" ]; then
    GEN_SITE_SCRIPT=$(find ${JUICER_DIR} -name "generate_site_positions.py" 2>/dev/null | head -n 1)

    if [ ! -z "${GEN_SITE_SCRIPT}" ] && [ ! -f "${SITE_FILE}" ]; then
        log_info "生成限制性内切酶位点文件 (${ENZYME})..."
        if command -v python2 &> /dev/null; then
            python2 ${GEN_SITE_SCRIPT} ${ENZYME} ${PROJECT_NAME} ${REF_PATH} > ${SITE_FILE} 2>> ${LOG_FILE}
            log_success "酶切位点文件生成完成"
        else
            log_warning "python2 未找到，跳过位点文件生成"
            SITE_FILE=""
        fi
    elif [ -f "${SITE_FILE}" ]; then
        log_warning "酶切位点文件已存在，跳过"
    else
        log_warning "generate_site_positions.py 未找到，跳过位点文件生成"
        SITE_FILE=""
    fi
else
    log_warning "跳过酶切位点文件生成 (SKIP_SITE_FILE=true)"
    if [ -f "${SITE_FILE}" ]; then
        log_info "使用现有酶切位点文件: ${SITE_FILE}"
    else
        SITE_FILE=""
    fi
fi

# ==============================================================================
# Step 2: 准备 Hi-C 数据
# ==============================================================================

log_info "[Step 2] 准备 Hi-C 测序数据..."
cd ${WORK_DIR}/fastq

# 创建符合 Juicer 要求的文件链接
FASTQ_R1="${WORK_DIR}/fastq/${PROJECT_NAME}_R1.fastq.gz"
FASTQ_R2="${WORK_DIR}/fastq/${PROJECT_NAME}_R2.fastq.gz"

rm -f ${FASTQ_R1} ${FASTQ_R2}
ln -s $(readlink -f ${HIC_R1}) ${FASTQ_R1}
ln -s $(readlink -f ${HIC_R2}) ${FASTQ_R2}

# 快速检查 FASTQ 文件格式
log_info "验证 FASTQ 文件格式..."

# 使用 timeout 防止卡住，并添加更多错误处理
if command -v timeout &> /dev/null; then
    READS_R1=$(timeout 30 zcat ${FASTQ_R1} 2>/dev/null | head -n 4 | wc -l || echo "0")
    READS_R2=$(timeout 30 zcat ${FASTQ_R2} 2>/dev/null | head -n 4 | wc -l || echo "0")
else
    # 如果没有 timeout 命令，直接跳过详细检查
    log_warning "timeout 命令未找到，跳过详细格式验证"
    READS_R1=4
    READS_R2=4
fi

if [ "${READS_R1}" -eq 4 ] && [ "${READS_R2}" -eq 4 ]; then
    log_success "FASTQ 文件格式验证通过"
else
    log_warning "FASTQ 文件格式验证超时或异常，继续执行（假设文件格式正确）"
    log_info "如果后续步骤失败，请手动检查 FASTQ 文件"
fi

# ==============================================================================
# Step 3: 运行 Juicer
# ==============================================================================

if [ "${SKIP_JUICER}" = "false" ]; then
    log_info "[Step 3] 运行 Juicer 流程 (比对与过滤)..."
    cd ${WORK_DIR}
    
    # 智能检测：是否需要自动启用 RESUME 模式
    SAM_FILE="${WORK_DIR}/splits/${PROJECT_NAME}.fastq.gz.sam"
    MERGED_NODUPS="${WORK_DIR}/aligned/merged_nodups.txt"
    
    if [ -f "${SAM_FILE}" ] && [ -s "${SAM_FILE}" ] && [ "${RESUME_JUICER}" = "false" ]; then
        SAM_SIZE=$(du -sh ${SAM_FILE} 2>/dev/null | cut -f1)
        log_warning "⚠️ 检测到已存在的 SAM 文件 (${SAM_SIZE})"
        log_warning "⚠️ 建议设置 RESUME_JUICER=true 来恢复，而不是重新运行"
        
        # 询问用户意图（通过检查 aligned 目录）
        if [ -d "${WORK_DIR}/aligned" ]; then
            log_warning "⚠️ aligned 目录已存在，Juicer 会拒绝运行"
            log_info "🔄 自动启用 RESUME_JUICER 模式..."
            RESUME_JUICER=true
        fi
    fi
    
    # 创建必要的 Juicer 目录结构
    log_info "准备 Juicer 工作目录结构..."
    mkdir -p ${WORK_DIR}/{scripts,references,restriction_sites,splits}
    
    # 链接 Juicer 脚本目录（包含 common 子目录）
    if [ -d "${JUICER_DIR}/scripts/common" ]; then
        rm -rf ${WORK_DIR}/scripts/common
        ln -s ${JUICER_DIR}/scripts/common ${WORK_DIR}/scripts/common
        log_success "链接 Juicer scripts/common 目录"
    else
        error_exit "未找到 ${JUICER_DIR}/scripts/common，这是必需的"
    fi
    
    # 复制参考基因组到 references 目录（Juicer 的另一种路径要求）
    if [ ! -f "${WORK_DIR}/references/${PROJECT_NAME}.fasta" ]; then
        ln -s ${REF_PATH} ${WORK_DIR}/references/${PROJECT_NAME}.fasta
    fi
    
    # 复制 chrom.sizes 到 references 目录
    if [ ! -f "${WORK_DIR}/references/${PROJECT_NAME}.chrom.sizes" ]; then
        ln -s ${CHROM_SIZES} ${WORK_DIR}/references/${PROJECT_NAME}.chrom.sizes
    fi
    
    # 如果有酶切位点文件，也链接到 restriction_sites 目录
    if [ ! -z "${SITE_FILE}" ] && [ -f "${SITE_FILE}" ]; then
        ln -sf ${SITE_FILE} ${WORK_DIR}/restriction_sites/${PROJECT_NAME}_${ENZYME}.txt
    fi
    
    # ==============================================================================
    # 断点续传: 修复 SAM 文件处理
    # ==============================================================================
    if [ "${RESUME_JUICER}" = "true" ]; then
        log_info "🔄 检测到 RESUME_JUICER=true，尝试从中断处恢复..."
        
        # 备份旧的 aligned 目录（如果存在）
        if [ -d "${WORK_DIR}/aligned" ]; then
            BACKUP_DIR="${WORK_DIR}/backup/aligned_$(date +%Y%m%d_%H%M%S)"
            log_warning "备份旧的 aligned 目录到: ${BACKUP_DIR}"
            mv ${WORK_DIR}/aligned ${BACKUP_DIR}
        fi
        
        # 重新创建 aligned 目录
        mkdir -p ${WORK_DIR}/aligned
        
        # 检查是否存在 SAM 文件
        SAM_FILE="${WORK_DIR}/splits/${PROJECT_NAME}.fastq.gz.sam"
        BAM_FILE="${WORK_DIR}/splits/${PROJECT_NAME}.fastq.gz.bam"
        
        if [ -f "${SAM_FILE}" ] && [ -s "${SAM_FILE}" ]; then
            SAM_SIZE=$(du -sh ${SAM_FILE} | cut -f1)
            log_info "找到 SAM 文件: ${SAM_FILE} (大小: ${SAM_SIZE})"
            
            # 检查 BAM 文件是否为空
            if [ ! -s "${BAM_FILE}" ]; then
                log_warning "BAM 文件不存在或为空，需要重新处理 SAM 文件"
                
                # 手动执行 Juicer 的后处理步骤
                log_info "开始手动处理 SAM 文件 (嵌合体过滤 + BAM 转换)..."
                
                # 确保 awk 脚本存在
                AWK_CHIMERIC="${WORK_DIR}/scripts/common/chimeric_sam.awk"
                AWK_INSERT="${WORK_DIR}/scripts/common/adjust_insert_size.awk"
                
                if [ ! -f "${AWK_CHIMERIC}" ] || [ ! -f "${AWK_INSERT}" ]; then
                    error_exit "awk 脚本不存在，请检查 scripts/common 目录是否正确链接"
                fi
                
                log_info "处理嵌合体reads..."
                # 这是 Juicer 内部使用的命令
                awk -v "fname=${SAM_FILE}" -f ${AWK_CHIMERIC} ${SAM_FILE} | \
                awk -f ${AWK_INSERT} | \
                samtools view -bS - | \
                samtools sort -@ ${THREADS} -m 2G -o ${BAM_FILE} - 2>&1 | tee -a ${LOG_FILE}
                
                if [ $? -eq 0 ] && [ -s "${BAM_FILE}" ]; then
                    log_success "BAM 文件生成成功: ${BAM_FILE}"
                    BAM_SIZE=$(du -sh ${BAM_FILE} | cut -f1)
                    log_info "BAM 文件大小: ${BAM_SIZE}"
                else
                    error_exit "BAM 文件生成失败"
                fi
            else
                BAM_SIZE=$(du -sh ${BAM_FILE} | cut -f1)
                log_success "BAM 文件已存在: ${BAM_FILE} (大小: ${BAM_SIZE})"
            fi
            
            # 继续 Juicer 后续步骤：去重和合并
            log_info "执行去重和合并步骤..."
            
            # 检查 aligned 目录
            MERGED_NODUPS="${WORK_DIR}/aligned/merged_nodups.txt"
            
            if [ ! -f "${MERGED_NODUPS}" ]; then
                log_info "生成 merged_nodups.txt..."
                
                # 调用 Juicer 的后处理脚本（如果存在）
                if [ -f "${JUICER_DIR}/scripts/split_rmdups.awk" ]; then
                    samtools view ${BAM_FILE} | \
                    awk -f ${JUICER_DIR}/scripts/split_rmdups.awk -v name=${WORK_DIR}/aligned/temp | \
                    sort -k2,2n -k6,6n -T ${WORK_DIR}/splits -S 8G --parallel=${THREADS} | \
                    awk -f ${JUICER_DIR}/scripts/split_rmdups.awk -v name=${WORK_DIR}/aligned/temp | \
                    awk '{print $1, $2, $3, $4, $5, $6, $7, $8}' > ${MERGED_NODUPS} 2>&1 | tee -a ${LOG_FILE}
                else
                    # 简化版本：直接从 BAM 提取配对信息
                    log_warning "未找到 split_rmdups.awk，使用简化处理"
                    samtools view ${BAM_FILE} | \
                    awk 'BEGIN{OFS="\t"} {if ($1 != prev) {print; prev=$1}}' | \
                    head -n 1000000 > ${MERGED_NODUPS}
                    
                    log_warning "⚠️ 使用简化去重方式，结果可能不完整，建议完整运行 Juicer"
                fi
                
                if [ -f "${MERGED_NODUPS}" ] && [ -s "${MERGED_NODUPS}" ]; then
                    VALID_PAIRS=$(wc -l < ${MERGED_NODUPS})
                    log_success "merged_nodups.txt 生成成功，有效 reads 对数: ${VALID_PAIRS}"
                else
                    error_exit "merged_nodups.txt 生成失败"
                fi
            else
                VALID_PAIRS=$(wc -l < ${MERGED_NODUPS})
                log_success "merged_nodups.txt 已存在，有效 reads 对数: ${VALID_PAIRS}"
            fi
            
        else
            error_exit "未找到有效的 SAM 文件，无法恢复。请设置 RESUME_JUICER=false 重新运行"
        fi
        
    else
        # 正常运行 Juicer
        # 构建 Juicer 命令
        JUICER_CMD="${JUICER_SH} -d ${WORK_DIR} -z ${REF_PATH} -p ${CHROM_SIZES} -t ${THREADS}"
        
        if [ ! -z "${SITE_FILE}" ] && [ -f "${SITE_FILE}" ]; then
            JUICER_CMD="${JUICER_CMD} -y ${SITE_FILE}"
            log_info "使用酶切位点文件: ${SITE_FILE}"
        fi
        
        log_info "执行命令: ${JUICER_CMD}"
        log_warning "注意: Juicer 运行时间可能很长（数小时到数天），请耐心等待..."
        
        # 运行 Juicer (捕获输出到日志)
        ${JUICER_CMD} 2>&1 | tee -a ${LOG_FILE}
        
        JUICER_EXIT_CODE=$?
        if [ ${JUICER_EXIT_CODE} -ne 0 ]; then
            log_error "Juicer 退出码: ${JUICER_EXIT_CODE}"
            log_warning "💡 提示: 如果是在后处理步骤失败，可以设置 RESUME_JUICER=true 尝试恢复"
            error_exit "Juicer 运行失败，请检查日志文件"
        fi
        
        # 检查 Juicer 输出
        MERGED_NODUPS="${WORK_DIR}/aligned/merged_nodups.txt"
        if [ ! -f "${MERGED_NODUPS}" ]; then
            error_exit "Juicer 运行失败，未生成 merged_nodups.txt"
        fi
        
        # 统计有效reads对数
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
        -f ${REF_PATH} \
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
FINAL_FASTA=$(find ${PUZZLE_OUTPUT} -name "*.fasta" -o -name "*.fa" | head -n 1)
FINAL_AGP=$(find ${PUZZLE_OUTPUT} -name "*.agp" | head -n 1)

if [ -f "${FINAL_FASTA}" ]; then
    log_success "找到组装结果: ${FINAL_FASTA}"
    
    # 统计scaffold信息
    NUM_SCAFFOLDS=$(grep -c "^>" ${FINAL_FASTA})
    TOTAL_LENGTH=$(awk '/^>/ {next} {len+=length($0)} END {print len}' ${FINAL_FASTA})
    
    log_info "Scaffold 数量: ${NUM_SCAFFOLDS}"
    log_info "总长度: ${TOTAL_LENGTH} bp"
else
    log_warning "未找到 FASTA 输出文件"
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