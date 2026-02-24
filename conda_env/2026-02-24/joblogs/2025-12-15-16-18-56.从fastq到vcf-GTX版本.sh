#!/bin/bash
# =================================================================
#   重测序全基因组变异检测全流程分析脚本
#   Version: 3.0 (Enhanced)
#   支持: 自动化质控 -> 比对 -> 变异检测 -> 过滤
#   新增: 断点续传、并行化优化、详细进度追踪、错误恢复
# =================================================================

# set -euo pipefail
IFS=$'\n\t'

# =================================================================
#               📝 用户配置区域 (User Configuration)
# =================================================================

# 1. 核心输入路径 (必须修改)
PROJECT_BASE="${PROJECT_BASE:-/share/org/YZWL/yzwl_lixg/tmp/liuchao_gwas/Q33_test_gtx}"
RAW_FASTQ_DIR="${PROJECT_BASE}/01.data/raw"
REF_GENOME_FA="${PROJECT_BASE}/01.data/genome/genome.fa"

# 2. 比对模式选择
# 2.1 选择Parabricks
# MAPPING_MODE="${MAPPING_MODE:-parabricks}"  # parabricks(GPU) 或 gtx(CPU)
# USE_GTX_WGS="${USE_GTX_WGS:-false}"         # 是否使用GTX WGS进行比对+变异检测

# 2.2 选择GTX
MAPPING_MODE="${MAPPING_MODE:-gtx}"  # parabricks(GPU) 或 gtx(CPU)
USE_GTX_WGS="${USE_GTX_WGS:-true}"   # 是否使用GTX WGS进行比对+变异检测

# 3. 工具路径
GTX_BIN="${GTX_BIN:-/share/org/YZWL/yzwl_lixg/software/gtx/bin/gtx}"
GTX_CMD_GEN_SCRIPT="${GTX_CMD_GEN_SCRIPT:-${HOME}/software/scripts/51.生成GTX按染色体合并gVCF的脚本.sh}"

# 4. 线程资源配置
THREADS_MAPPING="${THREADS_MAPPING:-88}"
THREADS_GTX="${THREADS_GTX:-88}"
THREADS_FILTER="${THREADS_FILTER:-88}"
MAX_PARALLEL_JOBS="${MAX_PARALLEL_JOBS:-4}"  # 新增: 并行任务数控制

# 5. 样本阈值
GATK_THRESHOLD=4        # < 4 使用 GATK
GTX_SINGLE_THRESHOLD=200 # < 200 使用 GTX 单机模式
GTX_WINDOW_SIZE=20000000 # GTX 分块窗口大小 (20Mb)

# 6. 过滤参数
SNP_MIN_DP="${SNP_MIN_DP:-5}"
SNP_MIN_QUAL="${SNP_MIN_QUAL:-30}"
INDEL_MIN_DP="${INDEL_MIN_DP:-5}"
INDEL_MIN_QUAL="${INDEL_MIN_QUAL:-30}"

# 7. GTX WGS 特定参数
GTX_PCR_INDEL_MODEL="${GTX_PCR_INDEL_MODEL:-CONSERVATIVE}"
GTX_MIN_CONFIDENCE="${GTX_MIN_CONFIDENCE:-30}"
GTX_MIN_BASE_QUAL="${GTX_MIN_BASE_QUAL:-20}"
GTX_PLOIDY="${GTX_PLOIDY:-2}"

# 8. 高级选项
ENABLE_CHECKPOINT="${ENABLE_CHECKPOINT:-false}"  # 启用断点续传
DRY_RUN="${DRY_RUN:-false}"                     # 测试模式
VERBOSE="${VERBOSE:-false}"                      # 详细输出
SKIP_QC="${SKIP_QC:-false}"                     # 跳过质控(如已完成)
SKIP_MAPPING="${SKIP_MAPPING:-false}"           # 跳过比对(如已完成)

# =================================================================
#               ⚙️ 系统路径规划
# =================================================================
CLEAN_FASTQ_DIR="${PROJECT_BASE}/01.data/clean"
MAPPING_DIR="${PROJECT_BASE}/02.mapping"
GVCF_DIR="${MAPPING_DIR}/vcf"
BAM_DIR="${MAPPING_DIR}/bam"              # 新增: BAM文件目录
JOINT_DIR="${PROJECT_BASE}/03.joint_calling"
FILTER_DIR="${PROJECT_BASE}/04.filtered_snp_indel"
SCRIPT_DIR="${PROJECT_BASE}/00.scripts"
LOG_DIR="${PROJECT_BASE}/99.logs"
CHECKPOINT_DIR="${PROJECT_BASE}/.checkpoints"  # 新增: 检查点目录
TMP_DIR="${PROJECT_BASE}/.tmp"                 # 新增: 临时文件目录

# 全局变量
FINAL_VCF_PATH=""
PIPELINE_START_TIME=$(date +%s)

# 日志配置 - 先创建目录
TIMESTAMP=$(date '+%Y%m%d_%H%M%S')
mkdir -p "${LOG_DIR}" 2>/dev/null || true
LOG_FILE="${LOG_DIR}/pipeline_${TIMESTAMP}.log"
ERROR_LOG="${LOG_DIR}/error_${TIMESTAMP}.log"
touch "${LOG_FILE}" "${ERROR_LOG}" 2>/dev/null || true

# 颜色代码
readonly GREEN='\033[0;32m'
readonly YELLOW='\033[1;33m'
readonly RED='\033[0;31m'
readonly BLUE='\033[0;34m'
readonly CYAN='\033[0;36m'
readonly MAGENTA='\033[0;35m'
readonly NC='\033[0m'

# =================================================================
#               🛠️ 核心工具函数库
# =================================================================

# 增强日志函数 - 支持多级别和文件输出
log_msg() {
    local level="$1"
    shift
    local msg="$*"
    local timestamp="$(date '+%Y-%m-%d %H:%M:%S')"
    local color="${NC}"
    
    case "${level}" in
        INFO)  color="${GREEN}" ;;
        WARN)  color="${YELLOW}" ;;
        ERROR) color="${RED}" ;;
        STEP)  color="${BLUE}" ;;
        DEBUG) color="${CYAN}" ;;
        SUCCESS) color="${MAGENTA}" ;;
    esac
    
    local formatted="[${level}] ${timestamp} - ${msg}"
    
    # 安全地写入日志文件
    if [ -f "${LOG_FILE}" ]; then
        echo -e "${color}${formatted}${NC}" | tee -a "${LOG_FILE}" 2>/dev/null || echo -e "${color}${formatted}${NC}"
    else
        echo -e "${color}${formatted}${NC}"
    fi
    
    # 错误同时写入错误日志
    if [ "${level}" = "ERROR" ] && [ -f "${ERROR_LOG}" ]; then
        echo "${formatted}" >> "${ERROR_LOG}" 2>/dev/null || true
    fi
    
    # 详细模式额外输出
    if [ "${VERBOSE}" = "true" ] && [ "${level}" = "DEBUG" ]; then
        echo "${formatted}" >> "${LOG_FILE}" 2>/dev/null || true
    fi
}

log_info()    { log_msg "INFO" "$@"; }
log_warn()    { log_msg "WARN" "$@"; }
log_error()   { log_msg "ERROR" "$@"; }
log_debug()   { log_msg "DEBUG" "$@"; }
log_success() { log_msg "SUCCESS" "$@"; }

log_step() {
    local msg="$*"
    if [ -f "${LOG_FILE}" ]; then
        echo "" | tee -a "${LOG_FILE}" 2>/dev/null || echo ""
        echo -e "${BLUE}========================================${NC}" | tee -a "${LOG_FILE}" 2>/dev/null || echo -e "${BLUE}========================================${NC}"
        echo -e "${BLUE}${msg}${NC}" | tee -a "${LOG_FILE}" 2>/dev/null || echo -e "${BLUE}${msg}${NC}"
        echo -e "${BLUE}========================================${NC}" | tee -a "${LOG_FILE}" 2>/dev/null || echo -e "${BLUE}========================================${NC}"
    else
        echo ""
        echo -e "${BLUE}========================================${NC}"
        echo -e "${BLUE}${msg}${NC}"
        echo -e "${BLUE}========================================${NC}"
    fi
}

# 增强命令检查 - 提供版本信息
check_command() {
    local cmd="$1"
    if ! command -v "${cmd}" &> /dev/null; then
        log_error "必需的命令未找到: ${cmd}"
        log_error "请确保已安装并添加到 PATH"
        exit 1
    fi
    
    if [ "${VERBOSE}" = "true" ]; then
        local version=$("${cmd}" --version 2>&1 | head -n1 || echo "版本未知")
        log_debug "${cmd}: ${version}"
    fi
}

# 增强文件检查
check_file() {
    local file="$1"
    local desc="${2:-文件}"
    if [ ! -f "${file}" ]; then
        log_error "${desc}不存在: ${file}"
        exit 1
    fi
    log_debug "✓ ${desc}: ${file}"
}

check_dir_not_empty() {
    local dir="$1"
    local desc="${2:-目录}"
    if [ ! -d "${dir}" ]; then
        log_error "${desc}不存在: ${dir}"
        exit 1
    fi
    if [ -z "$(ls -A "${dir}" 2>/dev/null)" ]; then
        log_error "${desc}为空: ${dir}"
        exit 1
    fi
    local count=$(find "${dir}" -maxdepth 1 -type f | wc -l)
    log_debug "✓ ${desc}: ${dir} (${count} 个文件)"
}

# 安全创建目录
safe_mkdir() {
    local dir="$1"
    if [ -d "${dir}" ]; then
        log_debug "目录已存在: ${dir}"
        return 0
    fi
    
    if mkdir -p "${dir}" 2>/dev/null; then
        log_debug "创建目录: ${dir}"
    else
        log_error "无法创建目录: ${dir}"
        exit 1
    fi
}

# 样本计数 - 增强版
count_samples() {
    local dir="$1"
    local pattern="$2"
    local count=0
    
    if [ -d "${dir}" ]; then
        count=$(find "${dir}" -name "${pattern}" -type f 2>/dev/null | wc -l)
    fi
    
    echo "${count}"
}

# 新增: 检查点管理
checkpoint_exists() {
    local step="$1"
    [ -f "${CHECKPOINT_DIR}/${step}.done" ]
}

checkpoint_create() {
    local step="$1"
    echo "$(date '+%Y-%m-%d %H:%M:%S')" > "${CHECKPOINT_DIR}/${step}.done"
    log_debug "创建检查点: ${step}"
}

checkpoint_remove() {
    local step="$1"
    rm -f "${CHECKPOINT_DIR}/${step}.done"
    log_debug "移除检查点: ${step}"
}

checkpoint_list() {
    if [ -d "${CHECKPOINT_DIR}" ]; then
        log_info "已完成的步骤:"
        for f in "${CHECKPOINT_DIR}"/*.done; do
            if [ -f "$f" ]; then
                local step=$(basename "$f" .done)
                local time=$(cat "$f")
                log_info "  ✓ ${step} (${time})"
            fi
        done
    fi
}

# 新增: 磁盘空间检查
check_disk_space() {
    local path="$1"
    local required_gb="${2:-100}"  # 默认需要100GB
    
    local available_kb=$(df -k "${path}" | tail -1 | awk '{print $4}')
    local available_gb=$((available_kb / 1024 / 1024))
    
    if [ "${available_gb}" -lt "${required_gb}" ]; then
        log_warn "磁盘空间不足: ${available_gb}GB 可用, 建议至少 ${required_gb}GB"
        log_warn "路径: ${path}"
        return 1
    fi
    
    log_debug "✓ 磁盘空间充足: ${available_gb}GB 可用"
    return 0
}

# 新增: 内存检查
check_memory() {
    local required_gb="${1:-64}"  # 默认需要64GB
    
    local total_mem_kb=$(grep MemTotal /proc/meminfo | awk '{print $2}')
    local total_mem_gb=$((total_mem_kb / 1024 / 1024))
    
    if [ "${total_mem_gb}" -lt "${required_gb}" ]; then
        log_warn "系统内存: ${total_mem_gb}GB, 建议至少 ${required_gb}GB"
        return 1
    fi
    
    log_debug "✓ 系统内存: ${total_mem_gb}GB"
    return 0
}

# 新增: 执行时间统计
show_elapsed_time() {
    local start_time="$1"
    local end_time=$(date +%s)
    local elapsed=$((end_time - start_time))
    
    local hours=$((elapsed / 3600))
    local minutes=$(((elapsed % 3600) / 60))
    local seconds=$((elapsed % 60))
    
    printf "%02d:%02d:%02d" "${hours}" "${minutes}" "${seconds}"
}

# 新增: 进度条
show_progress() {
    local current="$1"
    local total="$2"
    local desc="${3:-Processing}"
    
    local percent=$((current * 100 / total))
    local filled=$((percent / 2))
    local empty=$((50 - filled))
    
    printf "\r${desc}: ["
    printf "%${filled}s" | tr ' ' '='
    printf "%${empty}s" | tr ' ' ' '
    printf "] %d%% (%d/%d)" "${percent}" "${current}" "${total}"
}

# 新增: 清理函数
cleanup() {
    local exit_code=$?
    
    if [ "${exit_code}" -ne 0 ]; then
        log_error "脚本异常退出 (Exit Code: ${exit_code})"
        log_error "请查看错误日志: ${ERROR_LOG}"
    fi
    
    # 清理临时文件
    if [ -d "${TMP_DIR}" ]; then
        log_debug "清理临时文件..."
        rm -rf "${TMP_DIR}"
    fi
    
    # 显示总运行时间
    local total_time=$(show_elapsed_time "${PIPELINE_START_TIME}")
    log_info "总运行时间: ${total_time}"
    
    exit "${exit_code}"
}

# =================================================================
#               ✅ 增强预检查模块
# =================================================================

pre_flight_checks() {
    log_step "🔍 Step 0: 系统预检查"
    
    # 创建所有必需目录
    log_info "创建项目目录结构..."
    for dir in "${CLEAN_FASTQ_DIR}" "${MAPPING_DIR}" "${BAM_DIR}" "${GVCF_DIR}" \
               "${JOINT_DIR}" "${FILTER_DIR}" "${SCRIPT_DIR}" "${LOG_DIR}" \
               "${CHECKPOINT_DIR}" "${TMP_DIR}"; do
        safe_mkdir "${dir}"
    done
    
    # 检查必需命令
    log_info "检查必需工具..."
    local required_tools=(bwa samtools biopytools bcftools tabix python3)
    
    # 根据模式检查对应工具
    if [ "${USE_GTX_WGS}" = "true" ] || [ "${MAPPING_MODE}" = "gtx" ]; then
        required_tools+=(faketime)
        if [ -z "${GTX_BIN}" ] || [ ! -f "${GTX_BIN}" ]; then
            log_error "GTX模式需要GTX_BIN路径"
            exit 1
        fi
        check_file "${GTX_BIN}" "GTX可执行文件"
    else
        required_tools+=(gatk)
    fi
    
    for cmd in "${required_tools[@]}"; do
        check_command "${cmd}"
    done
    
    # 检查参考基因组
    log_info "检查参考基因组..."
    check_file "${REF_GENOME_FA}" "参考基因组"
    
    # 检查原始数据
    log_info "检查原始数据..."
    check_dir_not_empty "${RAW_FASTQ_DIR}" "原始数据目录"
    
    # 检查GTX相关工具
    if [ "${USE_GTX_WGS}" = "true" ] || [ "${MAPPING_MODE}" = "gtx" ]; then
        check_file "${GTX_BIN}" "GTX工具"
        if [ ! -z "${GTX_CMD_GEN_SCRIPT}" ] && [ "${GTX_CMD_GEN_SCRIPT}" != "skip" ]; then
            if [ -f "${GTX_CMD_GEN_SCRIPT}" ]; then
                check_file "${GTX_CMD_GEN_SCRIPT}" "GTX命令生成脚本"
            fi
        fi
    fi
    
    # 系统资源检查
    log_info "检查系统资源..."
    check_disk_space "${PROJECT_BASE}" 200
    check_memory 64
    
    # 显示已完成的检查点
    if [ "${ENABLE_CHECKPOINT}" = "true" ]; then
        checkpoint_list
    fi
    
    # 配置摘要
    log_info "配置摘要:"
    log_info "  项目路径: ${PROJECT_BASE}"
    log_info "  比对模式: ${MAPPING_MODE} $([ "${USE_GTX_WGS}" = "true" ] && echo "(GTX WGS完整流程)" || echo "")"
    log_info "  线程配置: Mapping=${THREADS_MAPPING}, GTX=${THREADS_GTX}, Filter=${THREADS_FILTER}"
    log_info "  样本阈值: GATK<${GATK_THRESHOLD}, GTX<${GTX_SINGLE_THRESHOLD}"
    log_info "  断点续传: ${ENABLE_CHECKPOINT}"
    log_info "  测试模式: ${DRY_RUN}"
    
    log_success "✅ 预检查通过"
    
    if [ "${DRY_RUN}" = "true" ]; then
        log_warn "⚠️  测试模式已启用，将不执行实际命令"
    fi
}

# =================================================================
#               📊 基因组索引模块 (优化)
# =================================================================

build_genome_index() {
    log_step "📊 Step 1: 构建基因组索引"
    
    local step_name="genome_index"
    if [ "${ENABLE_CHECKPOINT}" = "true" ] && checkpoint_exists "${step_name}"; then
        log_info "检查点已存在，跳过索引构建"
        return 0
    fi
    
    local index_start=$(date +%s)
    local need_index=false
    
    # # BWA索引
    # if [ ! -f "${REF_GENOME_FA}.bwt" ]; then
    #     log_info "构建 BWA 索引..."
    #     need_index=true
    #     if [ "${DRY_RUN}" = "false" ]; then
    #         bwa index "${REF_GENOME_FA}" 2>&1 | tee -a "${LOG_FILE}" || {
    #             log_error "BWA索引构建失败"
    #             return 1
    #         }
    #     fi
    # else
    #     log_info "✓ BWA 索引已存在"
    # fi
    
    # SAMtools索引
    if [ ! -f "${REF_GENOME_FA}.fai" ]; then
        log_info "构建 SAMtools 索引..."
        need_index=true
        if [ "${DRY_RUN}" = "false" ]; then
            samtools faidx "${REF_GENOME_FA}" 2>&1 | tee -a "${LOG_FILE}" || {
                log_error "SAMtools索引构建失败"
                return 1
            }
        fi
    else
        log_info "✓ SAMtools 索引已存在"
    fi
    
    # GATK字典
    local ref_dict="${REF_GENOME_FA%.fa}.dict"
    if [ ! -f "${ref_dict}" ]; then
        log_info "构建 GATK 字典..."
        need_index=true
        if [ "${DRY_RUN}" = "false" ]; then
            gatk CreateSequenceDictionary \
                -R "${REF_GENOME_FA}" \
                -O "${ref_dict}" 2>&1 | tee -a "${LOG_FILE}" || {
                log_error "GATK字典构建失败"
                return 1
            }
        fi
    else
        log_info "✓ GATK 字典已存在"
    fi

    # GTX索引
    faketime '2020-10-20 00:00:00' "${GTX_BIN}" index "${REF_GENOME_FA}" --tmp-dir ./
    
    if [ "${need_index}" = "false" ]; then
        log_info "所有索引均已存在"
    else
        local index_time=$(show_elapsed_time "${index_start}")
        log_success "✅ 索引构建完成 (耗时: ${index_time})"
    fi
    
    if [ "${ENABLE_CHECKPOINT}" = "true" ]; then
        checkpoint_create "${step_name}"
    fi
}

# =================================================================
#               🧹 质控模块 (优化)
# =================================================================

run_quality_control() {
    log_step "🧹 Step 2: 质量控制"
    
    local step_name="quality_control"
    if [ "${ENABLE_CHECKPOINT}" = "true" ] && checkpoint_exists "${step_name}"; then
        log_info "检查点已存在，跳过质控"
        return 0
    fi
    
    if [ "${SKIP_QC}" = "true" ]; then
        log_warn "用户指定跳过质控步骤"
        return 0
    fi
    
    local raw_count=$(count_samples "${RAW_FASTQ_DIR}" "*.fq.gz")
    log_info "检测到 ${raw_count} 个原始 FASTQ 文件"
    
    if [ "${raw_count}" -eq 0 ]; then
        log_error "未找到原始 FASTQ 文件 (*.fq.gz)"
        return 1
    fi
    
    log_info "开始质控处理..."
    local qc_start=$(date +%s)
    
    if [ "${DRY_RUN}" = "false" ]; then
        biopytools fastp \
            -i "${RAW_FASTQ_DIR}" \
            -o "${CLEAN_FASTQ_DIR}" \
            --read1-suffix "_1.fq.gz" \
            --read2-suffix "_2.fq.gz" 2>&1 | tee -a "${LOG_FILE}" || {
            log_error "质控处理失败"
            return 1
        }
    fi
    
    local clean_count=$(count_samples "${CLEAN_FASTQ_DIR}" "*.fq.gz")
    local qc_time=$(show_elapsed_time "${qc_start}")
    
    log_success "✅ 质控完成: ${clean_count} 个清洁文件 (耗时: ${qc_time})"
    
    if [ "${ENABLE_CHECKPOINT}" = "true" ]; then
        checkpoint_create "${step_name}"
    fi
}

# =================================================================
#               🗺️ 比对模块 (支持双模式)
# =================================================================

run_mapping() {
    log_step "🗺️ Step 3: 序列比对"
    
    local step_name="mapping"
    if [ "${ENABLE_CHECKPOINT}" = "true" ] && checkpoint_exists "${step_name}"; then
        log_info "检查点已存在，跳过比对"
        return 0
    fi
    
    if [ "${SKIP_MAPPING}" = "true" ]; then
        log_warn "用户指定跳过比对步骤"
        return 0
    fi
    
    # 根据USE_GTX_WGS标志选择流程
    if [ "${USE_GTX_WGS}" = "true" ]; then
        run_gtx_wgs_pipeline
    else
        run_standard_mapping
    fi
    
    if [ "${ENABLE_CHECKPOINT}" = "true" ]; then
        checkpoint_create "${step_name}"
    fi
}

# 标准比对流程 (Parabricks GPU)
run_standard_mapping() {
    log_info "使用标准Parabricks比对模式 (GPU加速)"
    log_info "使用 ${THREADS_MAPPING} 线程进行比对..."
    
    local mapping_start=$(date +%s)
    
    if [ "${DRY_RUN}" = "false" ]; then
        biopytools parabricks \
            -i "${CLEAN_FASTQ_DIR}" \
            -o "${MAPPING_DIR}" \
            -r "${REF_GENOME_FA}" \
            -t "${THREADS_MAPPING}" \
            --read1-pattern "*_1.clean.fq.gz" \
            --read2-pattern "*_2.clean.fq.gz" \
            --no-joint-calling 2>&1 | tee -a "${LOG_FILE}" || {
            log_error "Parabricks比对失败"
            return 1
        }
    fi
    
    local gvcf_count=$(count_samples "${GVCF_DIR}" "*.g.vcf.gz")
    local mapping_time=$(show_elapsed_time "${mapping_start}")
    
    log_success "✅ 比对完成: ${gvcf_count} 个 gVCF 文件 (耗时: ${mapping_time})"
}

# GTX WGS完整流程 (CPU优化，一步到位)
run_gtx_wgs_pipeline() {
    log_info "使用GTX WGS完整流程 (CPU优化，比对+变异检测一体化)"
    log_info "使用 ${THREADS_GTX} 线程处理..."
    
    # 查找所有R1文件
    log_info "搜索输入文件..."
    local r1_files=$(find "${CLEAN_FASTQ_DIR}" -name "*_1.clean.fq.gz" -o -name "*_1.fq.gz" | sort)
    
    if [ -z "${r1_files}" ]; then
        log_error "未找到任何 *_1.clean.fq.gz 或 *_1.fq.gz 文件"
        return 1
    fi
    
    local total_samples=$(echo "${r1_files}" | wc -l)
    log_info "找到 ${total_samples} 个样品需要处理"
    
    local current=0
    local failed_samples=()
    local success_count=0
    
    local gtx_start=$(date +%s)
    
    # 处理每个样品
    while IFS= read -r r1_file; do
        ((current++))
        
        # 提取样品名
        local sample_name=$(basename "${r1_file}")
        sample_name=${sample_name%_1.clean.fq.gz}
        sample_name=${sample_name%_1.fq.gz}
        
        # 构建R2文件路径
        local r2_file=""
        if [ -f "${CLEAN_FASTQ_DIR}/${sample_name}_2.clean.fq.gz" ]; then
            r2_file="${CLEAN_FASTQ_DIR}/${sample_name}_2.clean.fq.gz"
        elif [ -f "${CLEAN_FASTQ_DIR}/${sample_name}_2.fq.gz" ]; then
            r2_file="${CLEAN_FASTQ_DIR}/${sample_name}_2.fq.gz"
        else
            log_error "未找到样品 ${sample_name} 的R2文件"
            failed_samples+=("${sample_name}")
            continue
        fi
        
        # 定义输出文件
        local output_vcf="${GVCF_DIR}/${sample_name}.g.vcf.gz"
        local output_bam="${BAM_DIR}/${sample_name}.sorted.bam"
        
        # 检查是否已完成
        if [ -f "${output_vcf}" ] && [ -f "${output_bam}" ]; then
            log_info "[$current/$total_samples] 样品 ${sample_name} 已处理，跳过"
            ((success_count++))
            continue
        fi
        
        log_info "[$current/$total_samples] 处理样品: ${sample_name}"
        
        # 构建Read Group
        local read_group="@RG\\tID:${sample_name}\\tSM:${sample_name}\\tPL:ILLUMINA\\tLB:${sample_name}"
        
        if [ "${DRY_RUN}" = "false" ]; then
            # 运行GTX WGS
            if faketime '2020-10-20 00:00:00' "${GTX_BIN}" wgs \
                -R "${read_group}" \
                -o "${output_vcf}" \
                -b "${output_bam}" \
                -t "${THREADS_GTX}" \
                --tmp-dir "${TMP_DIR}" \
                --pcr-indel-model "${GTX_PCR_INDEL_MODEL}" \
                --standard-min-confidence-threshold-for-calling "${GTX_MIN_CONFIDENCE}" \
                --min-base-quality-score "${GTX_MIN_BASE_QUAL}" \
                --ploidy "${GTX_PLOIDY}" \
                "${REF_GENOME_FA}" \
                "${r1_file}" \
                "${r2_file}" 2>&1 | tee -a "${LOG_FILE}"; then
                
                log_success "  ✓ 样品 ${sample_name} 完成"
                ((success_count++))
                
                # 显示文件大小
                if [ -f "${output_vcf}" ]; then
                    local vcf_size=$(du -h "${output_vcf}" | cut -f1)
                    log_info "    VCF: ${vcf_size}"
                fi
                
                if [ -f "${output_bam}" ]; then
                    local bam_size=$(du -h "${output_bam}" | cut -f1)
                    log_info "    BAM: ${bam_size}"
                fi
            else
                log_error "  ✗ 样品 ${sample_name} 处理失败"
                failed_samples+=("${sample_name}")
            fi
        else
            log_info "  [DRY RUN] 跳过实际处理"
            ((success_count++))
        fi
        
        show_progress "${current}" "${total_samples}" "GTX WGS处理"
        
    done <<< "${r1_files}"
    
    echo ""  # 换行
    
    local gtx_time=$(show_elapsed_time "${gtx_start}")
    
    # 处理结果统计
    log_info "GTX WGS处理完成:"
    log_info "  成功: ${success_count}/${total_samples}"
    log_info "  失败: ${#failed_samples[@]}/${total_samples}"
    log_info "  耗时: ${gtx_time}"
    
    if [ ${#failed_samples[@]} -gt 0 ]; then
        log_warn "失败的样品:"
        for sample in "${failed_samples[@]}"; do
            log_warn "  - ${sample}"
        done
    fi
    
    # 如果有样品失败，返回错误
    if [ ${#failed_samples[@]} -gt 0 ]; then
        log_error "部分样品处理失败"
        return 1
    fi
    
    local gvcf_count=$(count_samples "${GVCF_DIR}" "*.g.vcf.gz")
    local bam_count=$(count_samples "${BAM_DIR}" "*.bam")
    log_success "✅ GTX WGS完成: ${gvcf_count} 个gVCF, ${bam_count} 个BAM"
}

# =================================================================
#               🧬 变异检测 - GATK模式
# =================================================================

run_gatk_joint_calling() {
    log_info "👉 使用 GATK GenotypeGVCFs 模式"
    
    local gatk_start=$(date +%s)
    
    if [ "${DRY_RUN}" = "false" ]; then
        biopytools gatk-joint \
            -i "${GVCF_DIR}" \
            -o "${JOINT_DIR}" \
            -r "${REF_GENOME_FA}" 2>&1 | tee -a "${LOG_FILE}" || {
            log_error "GATK联合检测失败"
            return 1
        }
    fi
    
    # 自动识别输出文件
    FINAL_VCF_PATH="${JOINT_DIR}/joint_genotyping_raw.vcf.gz"
    
    if [ -f "${JOINT_DIR}/joint_genotyping_merged_filtered.vcf.gz" ]; then
        FINAL_VCF_PATH="${JOINT_DIR}/joint_genotyping_merged_filtered.vcf.gz"
    fi
    
    if [ ! -f "${FINAL_VCF_PATH}" ] && [ "${DRY_RUN}" = "false" ]; then
        log_error "GATK 未生成预期的 VCF 文件"
        return 1
    fi
    
    local gatk_time=$(show_elapsed_time "${gatk_start}")
    log_success "GATK 输出: ${FINAL_VCF_PATH} (耗时: ${gatk_time})"
}

# =================================================================
#               🧬 变异检测 - GTX单机模式
# =================================================================

run_gtx_single_machine() {
    log_info "👉 使用 GTX 单机模式"
    
    FINAL_VCF_PATH="${JOINT_DIR}/gtx_joint_raw.vcf.gz"
    local tmp_dir="${TMP_DIR}/gtx"
    safe_mkdir "${tmp_dir}"
    
    local gtx_args=()
    gtx_args+=("-r" "${REF_GENOME_FA}")
    gtx_args+=("-o" "${FINAL_VCF_PATH}")
    gtx_args+=("-t" "${THREADS_GTX}")
    gtx_args+=("--tmp-dir" "${tmp_dir}")
    
    log_info "收集 gVCF 文件列表..."
    # local gvcf_count=0
    # while IFS= read -r gvcf_file; do
    #     gtx_args+=("-v" "${gvcf_file}")
    #     ((gvcf_count++))
    # done < <(find "${GVCF_DIR}" -name "*.g.vcf.gz" -type f)

    local gvcf_files=($(find "${GVCF_DIR}" -name "*.g.vcf.gz" -type f))
    for gvcf_file in "${gvcf_files[@]}"; do
        gtx_args+=("-v" "${gvcf_file}")
        ((gvcf_count++))
    done
    
    log_info "准备处理 ${gvcf_count} 个样本..."
    local gtx_start=$(date +%s)
    
    if [ "${DRY_RUN}" = "false" ]; then
        faketime '2020-10-20 00:00:00' "${GTX_BIN}" joint "${gtx_args[@]}" 2>&1 | tee -a "${LOG_FILE}" || {
            log_error "GTX联合检测失败"
            return 1
        }
    fi
    
    if [ ! -f "${FINAL_VCF_PATH}" ] && [ "${DRY_RUN}" = "false" ]; then
        log_error "GTX 未生成预期的 VCF 文件"
        return 1
    fi
    
    local gtx_time=$(show_elapsed_time "${gtx_start}")
    log_success "GTX 输出: ${FINAL_VCF_PATH} (耗时: ${gtx_time})"
}

# =================================================================
#               🧬 变异检测 - GTX集群模式
# =================================================================

generate_gtx_cluster_scripts() {
    log_warn "👉 大规模样本模式 (>= ${GTX_SINGLE_THRESHOLD})"
    
    local chunks_dir="${JOINT_DIR}/chunks"
    local gtx_job_script="${JOINT_DIR}/01.run_gtx_jobs.sh"
    local merge_py_script="${SCRIPT_DIR}/02.merge_vcf.py"
    local final_merged_vcf="${JOINT_DIR}/merged_all.vcf.gz"
    
    safe_mkdir "${chunks_dir}"
    check_file "${GTX_CMD_GEN_SCRIPT}" "GTX命令生成脚本"
    
    log_info "⚙️ 生成分块变异检测命令 (窗口: ${GTX_WINDOW_SIZE} bp)..."
    
    if [ "${DRY_RUN}" = "false" ]; then
        bash "${GTX_CMD_GEN_SCRIPT}" \
            -g "${GTX_BIN}" \
            -r "${REF_GENOME_FA}" \
            -i "${GVCF_DIR}" \
            -o "${chunks_dir}" \
            -w "${GTX_WINDOW_SIZE}" \
            -s "${gtx_job_script}" \
            -t "${THREADS_GTX}" 2>&1 | tee -a "${LOG_FILE}" || {
            log_error "GTX脚本生成失败"
            return 1
        }
        
        chmod +x "${gtx_job_script}"
    fi
    
    # 生成优化的Python合并脚本
    cat << 'PYTHON_SCRIPT' > "${merge_py_script}"
#!/usr/bin/env python3
"""
VCF合并脚本 - 支持自然排序和并行处理
"""
import os
import sys
import glob
import re
import subprocess
import tempfile
from pathlib import Path

def natural_sort_key(filename):
    """自然排序关键字函数"""
    basename = os.path.basename(filename)
    return [int(text) if text.isdigit() else text.lower() 
            for text in re.split(r'([0-9]+)', basename)]

def validate_vcf(vcf_file):
    """验证VCF文件完整性"""
    try:
        result = subprocess.run(
            ['bcftools', 'index', '--nrecords', vcf_file],
            capture_output=True,
            text=True,
            check=True
        )
        return True
    except subprocess.CalledProcessError:
        print(f"⚠️  警告: {vcf_file} 验证失败", file=sys.stderr)
        return False

def main():
    if len(sys.argv) != 3:
        print("Usage: python3 merge_vcf.py <input_dir> <output_vcf>", file=sys.stderr)
        sys.exit(1)
    
    input_dir = Path(sys.argv[1])
    output_file = Path(sys.argv[2])
    
    # 查找VCF文件
    vcf_pattern = input_dir / "*.joint.vcf.gz"
    vcf_files = sorted(glob.glob(str(vcf_pattern)), key=natural_sort_key)
    
    if not vcf_files:
        print(f"❌ 错误: 未找到 *.joint.vcf.gz 文件在 {input_dir}", file=sys.stderr)
        sys.exit(1)
    
    print(f"📊 发现 {len(vcf_files)} 个VCF文件")
    
    # 验证VCF文件
    print("🔍 验证VCF文件完整性...")
    valid_files = [f for f in vcf_files if validate_vcf(f)]
    
    if len(valid_files) != len(vcf_files):
        print(f"⚠️  警告: {len(vcf_files) - len(valid_files)} 个文件验证失败", file=sys.stderr)
        response = input("是否继续使用有效文件? (y/N): ")
        if response.lower() != 'y':
            sys.exit(1)
        vcf_files = valid_files
    
    print(f"✅ {len(vcf_files)} 个文件验证通过")
    
    # 创建文件列表
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as tmp:
        for vcf in vcf_files:
            tmp.write(f"{vcf}\n")
        list_path = tmp.name
    
    try:
        # 合并VCF
        print(f"🔗 合并VCF文件到: {output_file}")
        subprocess.check_call(
            f"bcftools concat -f {list_path} -a -O z -o {output_file} --threads 48",
            shell=True
        )
        
        # 创建索引
        print("📑 创建索引...")
        subprocess.check_call(f"tabix -p vcf {output_file}", shell=True)
        
        print(f"✅ 合并完成: {output_file}")
        
        # 显示统计信息
        result = subprocess.run(
            f"bcftools stats {output_file} | grep 'number of records:'",
            shell=True,
            capture_output=True,
            text=True
        )
        if result.stdout:
            print(f"📊 {result.stdout.strip()}")
            
    except subprocess.CalledProcessError as e:
        print(f"❌ 合并失败: {e}", file=sys.stderr)
        sys.exit(1)
    finally:
        os.remove(list_path)

if __name__ == "__main__":
    main()
PYTHON_SCRIPT
    chmod +x "${merge_py_script}"
    
    # 生成完整的操作指南
    cat << MANUAL_GUIDE | tee -a "${LOG_FILE}"

========================================================================
🛑 自动化流程已暂停 - 进入手动投递模式
========================================================================
样本数: >= ${GTX_SINGLE_THRESHOLD}
生成的脚本路径:
  - GTX任务脚本: ${gtx_job_script}
  - VCF合并脚本: ${merge_py_script}

📋 操作步骤:
------------------------------------------------------------------------
1️⃣  投递GTX任务到集群:
   batch_sub -i ${gtx_job_script} \\
             -j gtx_joint \\
             -s 5 \\
             -m 800

2️⃣  监控任务状态:
   batch_stat -j gtx_joint

3️⃣  任务完成后合并VCF:
   python3 ${merge_py_script} \\
           ${chunks_dir} \\
           ${final_merged_vcf}

4️⃣  验证合并结果:
   bcftools stats ${final_merged_vcf} | head -n 50

5️⃣  运行变异过滤:
   bash $(readlink -f $0) --resume-from-filtering

   或手动运行:
   biopytools filter-snp-indel \\
       -i ${final_merged_vcf} \\
       -o ${FILTER_DIR} \\
       -t ${THREADS_FILTER} \\
       --snp-dp ${SNP_MIN_DP} \\
       --indel-dp ${INDEL_MIN_DP}

========================================================================
💡 提示:
  - 任务脚本已生成,请检查后投递
  - 建议先投递1-2个任务测试
  - 可用 tail -f ${LOG_FILE} 查看日志
========================================================================

MANUAL_GUIDE
    
    return 2
}

# =================================================================
#               🧬 变异检测主控模块 (优化)
# =================================================================

run_joint_calling() {
    log_step "🧬 Step 4: 联合变异检测"
    
    local step_name="joint_calling"
    if [ "${ENABLE_CHECKPOINT}" = "true" ] && checkpoint_exists "${step_name}"; then
        log_info "检查点已存在，跳过联合检测"
        
        # 恢复VCF路径
        if [ -f "${JOINT_DIR}/gtx_joint_raw.vcf.gz" ]; then
            FINAL_VCF_PATH="${JOINT_DIR}/gtx_joint_raw.vcf.gz"
        elif [ -f "${JOINT_DIR}/joint_genotyping_raw.vcf.gz" ]; then
            FINAL_VCF_PATH="${JOINT_DIR}/joint_genotyping_raw.vcf.gz"
        fi
        
        return 0
    fi
    
    local sample_count=$(count_samples "${GVCF_DIR}" "*.g.vcf.gz")
    log_info "检测到 ${sample_count} 个 gVCF 样本"
    
    if [ "${sample_count}" -eq 0 ]; then
        log_error "未找到任何 gVCF 文件"
        return 1
    fi
    
    # 策略选择
    log_info "样本数分析:"
    log_info "  < ${GATK_THRESHOLD} → GATK模式"
    log_info "  ${GATK_THRESHOLD}-${GTX_SINGLE_THRESHOLD} → GTX单机模式"
    log_info "  >= ${GTX_SINGLE_THRESHOLD} → GTX集群模式"
    
    local jc_result=0
    
    if [ "${sample_count}" -ge "${GTX_SINGLE_THRESHOLD}" ]; then
        generate_gtx_cluster_scripts
        jc_result=$?
        
    elif [ "${sample_count}" -ge "${GATK_THRESHOLD}" ]; then
        run_gtx_single_machine
        jc_result=$?
        
    else
        run_gatk_joint_calling
        jc_result=$?
    fi
    
    if [ "${jc_result}" -ne 0 ] && [ "${jc_result}" -ne 2 ]; then
        log_error "变异检测失败"
        return "${jc_result}"
    fi
    
    if [ "${jc_result}" -eq 0 ]; then
        log_success "✅ 联合检测完成"
        if [ "${ENABLE_CHECKPOINT}" = "true" ]; then
            checkpoint_create "${step_name}"
        fi
    fi
    
    return "${jc_result}"
}

# =================================================================
#               🧹 变异过滤模块 (优化)
# =================================================================

run_variant_filtering() {
    local input_vcf="$1"
    
    log_step "🧹 Step 5: 变异过滤"
    
    local step_name="variant_filtering"
    if [ "${ENABLE_CHECKPOINT}" = "true" ] && checkpoint_exists "${step_name}"; then
        log_info "检查点已存在，跳过过滤"
        return 0
    fi
    
    if [ -z "${input_vcf}" ] || [ ! -f "${input_vcf}" ]; then
        log_error "过滤输入文件无效: ${input_vcf}"
        return 1
    fi
    
    log_info "输入 VCF: ${input_vcf}"
    log_info "过滤参数:"
    log_info "  SNP  - 最小深度: ${SNP_MIN_DP}, 最小质量: ${SNP_MIN_QUAL}"
    log_info "  InDel - 最小深度: ${INDEL_MIN_DP}, 最小质量: ${INDEL_MIN_QUAL}"
    
    local filter_start=$(date +%s)
    
    if [ "${DRY_RUN}" = "false" ]; then
        biopytools filter-snp-indel \
            -i "${input_vcf}" \
            -o "${FILTER_DIR}" \
            -t "${THREADS_FILTER}" \
            --snp-dp "${SNP_MIN_DP}" \
            --indel-dp "${INDEL_MIN_DP}" 2>&1 | tee -a "${LOG_FILE}" || {
            log_error "变异过滤失败"
            return 1
        }
    fi
    
    local filter_time=$(show_elapsed_time "${filter_start}")
    log_success "✅ 过滤完成 (耗时: ${filter_time})"
    
    # 显示结果统计
    if [ "${DRY_RUN}" = "false" ]; then
        log_info "输出文件统计:"
        for vcf in "${FILTER_DIR}"/*.vcf.gz; do
            if [ -f "$vcf" ]; then
                local count=$(bcftools view -H "$vcf" | wc -l)
                log_info "  $(basename "$vcf"): ${count} 个变异"
            fi
        done
    fi
    
    if [ "${ENABLE_CHECKPOINT}" = "true" ]; then
        checkpoint_create "${step_name}"
    fi
}

# =================================================================
#               📊 最终报告生成
# =================================================================

generate_final_report() {
    log_step "📊 生成分析报告"
    
    local report_file="${PROJECT_BASE}/ANALYSIS_REPORT.txt"
    local total_time=$(show_elapsed_time "${PIPELINE_START_TIME}")
    
    cat > "${report_file}" << REPORT
========================================================================
             全基因组重测序分析流程 - 最终报告
========================================================================
分析日期: $(date '+%Y-%m-%d %H:%M:%S')
项目路径: ${PROJECT_BASE}
总运行时间: ${total_time}

------------------------------------------------------------------------
📁 输入数据
------------------------------------------------------------------------
原始FASTQ目录: ${RAW_FASTQ_DIR}
参考基因组: ${REF_GENOME_FA}
样本数量: $(count_samples "${GVCF_DIR}" "*.g.vcf.gz")

------------------------------------------------------------------------
⚙️ 处理参数
------------------------------------------------------------------------
比对线程: ${THREADS_MAPPING}
GTX线程: ${THREADS_GTX}
过滤线程: ${THREADS_FILTER}
样本阈值: GATK<${GATK_THRESHOLD}, GTX<${GTX_SINGLE_THRESHOLD}

过滤参数:
  - SNP最小深度: ${SNP_MIN_DP}
  - SNP最小质量: ${SNP_MIN_QUAL}
  - InDel最小深度: ${INDEL_MIN_DP}
  - InDel最小质量: ${INDEL_MIN_QUAL}

------------------------------------------------------------------------
📂 输出目录
------------------------------------------------------------------------
清洁数据: ${CLEAN_FASTQ_DIR}
比对结果: ${MAPPING_DIR}
变异检测: ${JOINT_DIR}
过滤结果: ${FILTER_DIR}
日志文件: ${LOG_DIR}

------------------------------------------------------------------------
📊 结果文件
------------------------------------------------------------------------
REPORT

    if [ -n "${FINAL_VCF_PATH}" ] && [ -f "${FINAL_VCF_PATH}" ]; then
        echo "原始VCF: ${FINAL_VCF_PATH}" >> "${report_file}"
    fi
    
    if [ -d "${FILTER_DIR}" ]; then
        echo "" >> "${report_file}"
        echo "过滤后VCF文件:" >> "${report_file}"
        for vcf in "${FILTER_DIR}"/*.vcf.gz; do
            if [ -f "$vcf" ]; then
                local size=$(du -h "$vcf" | cut -f1)
                local count=$(bcftools view -H "$vcf" 2>/dev/null | wc -l || echo "N/A")
                echo "  - $(basename "$vcf"): ${size}, ${count} 变异" >> "${report_file}"
            fi
        done
    fi
    
    cat >> "${report_file}" << REPORT

------------------------------------------------------------------------
✅ 已完成步骤
------------------------------------------------------------------------
REPORT

    if [ -d "${CHECKPOINT_DIR}" ]; then
        for checkpoint in "${CHECKPOINT_DIR}"/*.done; do
            if [ -f "$checkpoint" ]; then
                local step=$(basename "$checkpoint" .done)
                local time=$(cat "$checkpoint")
                echo "  ✓ ${step} (${time})" >> "${report_file}"
            fi
        done
    fi
    
    cat >> "${report_file}" << REPORT

------------------------------------------------------------------------
📝 日志文件
------------------------------------------------------------------------
主日志: ${LOG_FILE}
错误日志: ${ERROR_LOG}

========================================================================
              分析完成 - 感谢使用本流程
========================================================================
REPORT

    log_info "报告已生成: ${report_file}"
    cat "${report_file}"
}

# =================================================================
#               🎯 主流程入口 (优化)
# =================================================================

main() {
    log_step "✨ 全基因组重测序自动化分析流程 v3.0 ✨"
    log_info "项目: ${PROJECT_BASE}"
    log_info "日志: ${LOG_FILE}"
    log_info "开始时间: $(date '+%Y-%m-%d %H:%M:%S')"
    
    # 参数解析
    local resume_from=""
    while [[ $# -gt 0 ]]; do
        case $1 in
            --resume-from-filtering)
                resume_from="filtering"
                shift
                ;;
            --reset-checkpoints)
                log_warn "重置所有检查点"
                rm -rf "${CHECKPOINT_DIR}"
                shift
                ;;
            --show-checkpoints)
                checkpoint_list
                exit 0
                ;;
            --help|-h)
                cat << HELP
用法: $0 [选项]

选项:
  --resume-from-filtering  从过滤步骤恢复(用于GTX集群模式完成后)
  --reset-checkpoints      重置所有检查点,重新运行
  --show-checkpoints       显示已完成的检查点
  --help, -h               显示此帮助信息

环境变量:
  PROJECT_BASE            项目根目录
  MAPPING_MODE            比对模式: parabricks(默认,GPU) 或 gtx(CPU)
  USE_GTX_WGS             使用GTX WGS完整流程(比对+变异检测) (默认: false)
  THREADS_MAPPING         比对线程数 (默认: 88)
  THREADS_GTX             GTX线程数 (默认: 88)
  THREADS_FILTER          过滤线程数 (默认: 88)
  ENABLE_CHECKPOINT       启用断点续传 (默认: true)
  DRY_RUN                 测试模式 (默认: false)
  VERBOSE                 详细输出 (默认: false)
  SKIP_QC                 跳过质控 (默认: false)
  SKIP_MAPPING            跳过比对 (默认: false)
  
GTX WGS参数:
  GTX_PCR_INDEL_MODEL     PCR InDel模型 (默认: CONSERVATIVE)
  GTX_MIN_CONFIDENCE      最小置信度 (默认: 30)
  GTX_MIN_BASE_QUAL       最小碱基质量 (默认: 20)
  GTX_PLOIDY              倍性 (默认: 2)

示例:
  # 正常运行 (GPU模式)
  bash $0

  # 使用CPU模式 (GTX WGS一体化)
  USE_GTX_WGS=true bash $0
  
  # GPU不可用时的替代方案
  MAPPING_MODE=gtx USE_GTX_WGS=true bash $0

  # 测试模式
  DRY_RUN=true bash $0

  # GTX集群模式完成后继续
  bash $0 --resume-from-filtering

  # 重新运行全流程
  bash $0 --reset-checkpoints
HELP
                exit 0
                ;;
            *)
                log_error "未知选项: $1"
                exit 1
                ;;
        esac
    done
    
    # 如果是从过滤步骤恢复
    if [ "${resume_from}" = "filtering" ]; then
        log_info "从过滤步骤恢复..."
        
        local merged_vcf="${JOINT_DIR}/merged_all.vcf.gz"
        if [ ! -f "${merged_vcf}" ]; then
            log_error "未找到合并的VCF文件: ${merged_vcf}"
            log_error "请先完成GTX集群任务并运行合并脚本"
            exit 1
        fi
        
        run_variant_filtering "${merged_vcf}"
        generate_final_report
        exit 0
    fi
    
    # 正常流程
    pre_flight_checks
    build_genome_index
    run_quality_control
    run_mapping
    
    # 变异检测
    run_joint_calling
    local jc_status=$?
    
    # 如果使用GTX WGS模式，跳过联合检测和过滤（已经生成单样本VCF）
    if [ "${USE_GTX_WGS}" = "true" ]; then
        log_info "GTX WGS模式已生成单样本VCF文件，跳过联合检测"
        log_info "如需联合检测，请设置 USE_GTX_WGS=false 重新运行"
        generate_final_report
        
        log_step "🎉 GTX WGS流程执行成功！"
        log_info "📂 gVCF目录: ${GVCF_DIR}"
        log_info "📂 BAM目录: ${BAM_DIR}"
        log_info "📊 分析报告: ${PROJECT_BASE}/ANALYSIS_REPORT.txt"
        log_success "总运行时间: $(show_elapsed_time "${PIPELINE_START_TIME}")"
        exit 0
    fi
    
    # 处理不同的返回状态
    if [ "${jc_status}" -eq 2 ]; then
        log_warn "流程已生成集群任务脚本"
        log_warn "请按照提示完成后续步骤"
        exit 0
        
    elif [ "${jc_status}" -ne 0 ]; then
        log_error "变异检测失败 (Exit Code: ${jc_status})"
        exit 1
    fi
    
    # 变异过滤
    if [ -n "${FINAL_VCF_PATH}" ]; then
        run_variant_filtering "${FINAL_VCF_PATH}"
    else
        log_error "未找到VCF文件路径"
        exit 1
    fi
    
    # 生成最终报告
    generate_final_report
    
    log_step "🎉 全流程执行成功！"
    log_info "📂 结果目录: ${FILTER_DIR}"
    log_info "📊 分析报告: ${PROJECT_BASE}/ANALYSIS_REPORT.txt"
    log_info "📝 日志文件: ${LOG_FILE}"
    log_success "总运行时间: $(show_elapsed_time "${PIPELINE_START_TIME}")"
}

# =================================================================
#               🚀 脚本执行
# =================================================================

# 设置陷阱处理
trap cleanup EXIT
trap 'log_error "脚本被用户中断 (Ctrl+C)"; exit 130' INT
trap 'log_error "脚本被终止信号中断"; exit 143' TERM

# 执行主函数
main "$@"