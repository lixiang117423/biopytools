#!/bin/bash

# ==============================================================================
# 脚本名称: run_glnexus_by_chrom.sh
# 功能: 按染色体分步运行 GLnexus 联合变异检测，防止内存溢出，最后合并结果
# 版本: 2.0
# ==============================================================================

set -e
set -o pipefail

# --- [颜色定义] ---
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# --- [日志函数] ---
log_info() { echo -e "${BLUE}[INFO]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_success() { echo -e "${GREEN}[SUCCESS]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_warn() { echo -e "${YELLOW}[WARNING]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_error() { echo -e "${RED}[ERROR]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1"; }

# --- [用户配置区域] ---
REF_GENOME="${REF_GENOME:-/share/org/YZWL/yzwl_lixg/project/21.野生大豆群体/01.data/genome/genome.fa}"
GVCF_DIR="${GVCF_DIR:-/share/org/YZWL/yzwl_lixg/project/21.野生大豆群体/03.glnexus_joint/vcf}"
OUTPUT_DIR="${OUTPUT_DIR:-/share/org/YZWL/yzwl_lixg/project/21.野生大豆群体/03.glnexus_joint}"
GLNEXUS_CONFIG="${GLNEXUS_CONFIG:-gatk}"

# 线程设置
INDEX_THREADS="${INDEX_THREADS:-16}"
GLNEXUS_THREADS="${GLNEXUS_THREADS:-32}"
CONVERT_THREADS="${CONVERT_THREADS:-64}"

# 高级选项
KEEP_INTERMEDIATE="${KEEP_INTERMEDIATE:-false}"  # 是否保留中间文件
RESUME_MODE="${RESUME_MODE:-false}"              # 断点续传模式
PARALLEL_CHROMS="${PARALLEL_CHROMS:-1}"          # 并行处理的染色体数（1=串行）
MIN_CHROM_LENGTH="${MIN_CHROM_LENGTH:-0}"        # 忽略小于此长度的染色体
MEMORY_LIMIT_GB="${MEMORY_LIMIT_GB:-0}"          # 内存限制(GB)，0表示不限制
GENERATE_UNCOMPRESSED="${GENERATE_UNCOMPRESSED:-false}"  # 是否生成未压缩的VCF

# ==============================================================================

# --- [帮助信息] ---
show_help() {
    cat << EOF
用法: $0 [选项]

选项:
  -r, --ref-genome PATH       参考基因组路径 (默认: 配置中的路径)
  -g, --gvcf-dir PATH         gVCF 文件目录 (默认: 配置中的路径)
  -o, --output-dir PATH       输出目录 (默认: 配置中的路径)
  -c, --config NAME           GLnexus 配置 (默认: gatk)
  -t, --threads NUM           GLnexus 线程数 (默认: 32)
  -p, --parallel NUM          并行处理染色体数 (默认: 1)
  -k, --keep-intermediate     保留中间文件
  -R, --resume                断点续传模式
  -m, --min-length NUM        最小染色体长度 (默认: 0)
  -u, --uncompressed          同时生成未压缩的 VCF 文件
  -h, --help                  显示此帮助信息

环境变量:
  所有配置项都可以通过环境变量设置，例如: REF_GENOME, GVCF_DIR 等

示例:
  $0 -t 64 -p 2 -k
  REF_GENOME=/path/to/ref.fa $0 --resume

EOF
    exit 0
}

# --- [参数解析] ---
while [[ $# -gt 0 ]]; do
    case $1 in
        -r|--ref-genome) REF_GENOME="$2"; shift 2 ;;
        -g|--gvcf-dir) GVCF_DIR="$2"; shift 2 ;;
        -o|--output-dir) OUTPUT_DIR="$2"; shift 2 ;;
        -c|--config) GLNEXUS_CONFIG="$2"; shift 2 ;;
        -t|--threads) GLNEXUS_THREADS="$2"; shift 2 ;;
        -p|--parallel) PARALLEL_CHROMS="$2"; shift 2 ;;
        -k|--keep-intermediate) KEEP_INTERMEDIATE=true; shift ;;
        -R|--resume) RESUME_MODE=true; shift ;;
        -m|--min-length) MIN_CHROM_LENGTH="$2"; shift 2 ;;
        -u|--uncompressed) GENERATE_UNCOMPRESSED=true; shift ;;
        -h|--help) show_help ;;
        *) log_error "未知参数: $1"; show_help ;;
    esac
done

# ==============================================================================
# --- [环境检查] ---
# ==============================================================================

log_info "开始环境检查..."

check_command() {
    if ! command -v "$1" &> /dev/null; then
        log_error "$1 未找到，请先安装"
        exit 1
    fi
    log_success "$1 已找到: $(command -v $1)"
}

check_command bcftools
check_command glnexus_cli
check_command samtools

# 检查文件和目录
if [ ! -f "${REF_GENOME}" ]; then
    log_error "参考基因组文件不存在: ${REF_GENOME}"
    exit 1
fi

if [ ! -d "${GVCF_DIR}" ]; then
    log_error "gVCF 目录不存在: ${GVCF_DIR}"
    exit 1
fi

# 统计 gVCF 文件数量
GVCF_COUNT=$(find "${GVCF_DIR}" -name "*.g.vcf.gz" | wc -l)
if [ "${GVCF_COUNT}" -eq 0 ]; then
    log_error "在 ${GVCF_DIR} 中未找到任何 .g.vcf.gz 文件"
    exit 1
fi
log_success "找到 ${GVCF_COUNT} 个 gVCF 文件"

# ==============================================================================
# --- [准备工作] ---
# ==============================================================================

log_info "准备工作目录..."
mkdir -p "${OUTPUT_DIR}"
TMP_OUT_DIR="${OUTPUT_DIR}/per_chrom_output"
mkdir -p "${TMP_OUT_DIR}"

# 日志文件
LOG_DIR="${OUTPUT_DIR}/logs"
mkdir -p "${LOG_DIR}"
MAIN_LOG="${LOG_DIR}/glnexus_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${MAIN_LOG}") 2>&1

log_info "日志文件: ${MAIN_LOG}"

VCF_OUT="${OUTPUT_DIR}/merged_output.vcf.gz"
BCF_LIST_FILE="${OUTPUT_DIR}/bcf_file_list.txt"
PROGRESS_FILE="${OUTPUT_DIR}/.progress"

# 断点续传: 读取已完成的染色体
if [ "${RESUME_MODE}" = true ] && [ -f "${PROGRESS_FILE}" ]; then
    log_info "启用断点续传模式，读取进度文件..."
    mapfile -t COMPLETED_CHROMS < "${PROGRESS_FILE}"
    log_info "已完成 ${#COMPLETED_CHROMS[@]} 条染色体"
else
    COMPLETED_CHROMS=()
    > "${PROGRESS_FILE}"
fi

# 初始化 BCF 列表
if [ "${RESUME_MODE}" = false ]; then
    > "${BCF_LIST_FILE}"
fi

# ==============================================================================
# --- [参考基因组索引] ---
# ==============================================================================

if [ ! -f "${REF_GENOME}.fai" ]; then
    log_info "为参考基因组创建 .fai 索引..."
    samtools faidx "${REF_GENOME}"
    log_success "索引创建完成"
else
    log_success "参考基因组索引已存在"
fi

# ==============================================================================
# --- [系统资源监控] ---
# ==============================================================================

check_system_resources() {
    local chrom=$1
    
    # 检查可用内存
    if [ "${MEMORY_LIMIT_GB}" -gt 0 ]; then
        local available_mem=$(free -g | awk '/^Mem:/{print $7}')
        if [ "${available_mem}" -lt "${MEMORY_LIMIT_GB}" ]; then
            log_warn "可用内存不足 (${available_mem}GB < ${MEMORY_LIMIT_GB}GB)，等待60秒..."
            sleep 60
            return 1
        fi
    fi
    
    # 检查磁盘空间 (至少需要 50GB)
    local available_space=$(df -BG "${OUTPUT_DIR}" | awk 'NR==2 {print $4}' | sed 's/G//')
    if [ "${available_space}" -lt 50 ]; then
        log_error "磁盘空间不足 (< 50GB)，请清理空间后重试"
        exit 1
    fi
    
    return 0
}

# ==============================================================================
# --- [处理单条染色体] ---
# ==============================================================================

process_chromosome() {
    local CHROM=$1
    local LENGTH=$2
    
    # 检查是否已完成
    if [[ " ${COMPLETED_CHROMS[@]} " =~ " ${CHROM} " ]]; then
        log_info "染色体 ${CHROM} 已处理过，跳过"
        return 0
    fi
    
    # 检查染色体长度
    if [ "${LENGTH}" -lt "${MIN_CHROM_LENGTH}" ]; then
        log_warn "染色体 ${CHROM} 长度 ${LENGTH} 小于最小值 ${MIN_CHROM_LENGTH}，跳过"
        return 0
    fi
    
    log_info "========== 处理染色体: ${CHROM} (长度: ${LENGTH} bp) =========="
    
    local START_TIME=$(date +%s)
    local CURRENT_BED="${TMP_OUT_DIR}/${CHROM}.bed"
    local CURRENT_DB="${TMP_OUT_DIR}/GLnexus.DB.${CHROM}"
    local CURRENT_BCF="${TMP_OUT_DIR}/${CHROM}.bcf"
    local CHROM_LOG="${LOG_DIR}/${CHROM}.log"
    
    # 等待系统资源
    while ! check_system_resources "${CHROM}"; do
        sleep 10
    done
    
    # 生成 BED 文件
    echo -e "${CHROM}\t0\t${LENGTH}" > "${CURRENT_BED}"
    log_success "BED 文件已创建: ${CURRENT_BED}"
    
    # 清理旧 DB
    if [ -d "${CURRENT_DB}" ]; then
        log_warn "删除旧数据库: ${CURRENT_DB}"
        rm -rf "${CURRENT_DB}"
    fi
    
    # 运行 GLnexus
    log_info "启动 GLnexus (线程: ${GLNEXUS_THREADS})..."
    
    if glnexus_cli \
        --config "${GLNEXUS_CONFIG}" \
        --bed "${CURRENT_BED}" \
        --dir "${CURRENT_DB}" \
        --threads "${GLNEXUS_THREADS}" \
        "${GVCF_DIR}"/*.g.vcf.gz \
        > "${CURRENT_BCF}" 2> "${CHROM_LOG}"; then
        
        # 验证输出
        if [ -s "${CURRENT_BCF}" ]; then
            local BCF_SIZE=$(du -h "${CURRENT_BCF}" | cut -f1)
            local END_TIME=$(date +%s)
            local ELAPSED=$((END_TIME - START_TIME))
            
            log_success "染色体 ${CHROM} 处理成功 (用时: ${ELAPSED}s, 大小: ${BCF_SIZE})"
            
            # 记录到文件列表
            echo "${CURRENT_BCF}" >> "${BCF_LIST_FILE}"
            
            # 记录进度
            echo "${CHROM}" >> "${PROGRESS_FILE}"
            
            # 清理临时文件
            if [ "${KEEP_INTERMEDIATE}" = false ]; then
                log_info "清理临时文件..."
                rm -rf "${CURRENT_DB}"
                rm -f "${CURRENT_BED}"
            fi
            
            return 0
        else
            log_error "染色体 ${CHROM} 输出文件为空"
            return 1
        fi
    else
        log_error "染色体 ${CHROM} GLnexus 运行失败，详见日志: ${CHROM_LOG}"
        return 1
    fi
}

# ==============================================================================
# --- [主处理循环] ---
# ==============================================================================

log_info "=========================================="
log_info " 开始按染色体处理 (并行数: ${PARALLEL_CHROMS})"
log_info "=========================================="

# 读取染色体列表
CHROM_LIST=()
while read -r CHROM LENGTH REST; do
    CHROM_LIST+=("${CHROM}:${LENGTH}")
done < "${REF_GENOME}.fai"

TOTAL_CHROMS=${#CHROM_LIST[@]}
log_info "共有 ${TOTAL_CHROMS} 条染色体需要处理"

# 处理染色体
PROCESSED=0
FAILED=0

if [ "${PARALLEL_CHROMS}" -eq 1 ]; then
    # 串行处理
    for CHROM_INFO in "${CHROM_LIST[@]}"; do
        IFS=':' read -r CHROM LENGTH <<< "${CHROM_INFO}"
        if process_chromosome "${CHROM}" "${LENGTH}"; then
            ((PROCESSED++))
        else
            ((FAILED++))
        fi
    done
else
    # 并行处理
    log_info "使用并行模式 (最多同时处理 ${PARALLEL_CHROMS} 条染色体)"
    
    for CHROM_INFO in "${CHROM_LIST[@]}"; do
        IFS=':' read -r CHROM LENGTH <<< "${CHROM_INFO}"
        
        # 等待直到有空闲槽位
        while [ $(jobs -r | wc -l) -ge "${PARALLEL_CHROMS}" ]; do
            sleep 5
        done
        
        # 后台处理
        (
            if process_chromosome "${CHROM}" "${LENGTH}"; then
                echo "SUCCESS:${CHROM}"
            else
                echo "FAILED:${CHROM}"
            fi
        ) &
    done
    
    # 等待所有任务完成
    wait
    
    # 统计结果
    PROCESSED=$(wc -l < "${PROGRESS_FILE}")
    FAILED=$((TOTAL_CHROMS - PROCESSED))
fi

log_info "=========================================="
log_info " 染色体处理完成: 成功 ${PROCESSED}, 失败 ${FAILED}"
log_info "=========================================="

if [ "${FAILED}" -gt 0 ]; then
    log_error "有 ${FAILED} 条染色体处理失败"
    exit 1
fi

# ==============================================================================
# --- [合并结果] ---
# ==============================================================================

if [ ! -s "${BCF_LIST_FILE}" ]; then
    log_error "没有生成任何 BCF 文件，无法合并"
    exit 1
fi

BCF_COUNT=$(wc -l < "${BCF_LIST_FILE}")
log_info "=========================================="
log_info " 开始合并 ${BCF_COUNT} 个 BCF 文件"
log_info "=========================================="

log_info "使用 bcftools concat 合并 (线程: ${CONVERT_THREADS})..."

if bcftools concat \
    --threads "${CONVERT_THREADS}" \
    -f "${BCF_LIST_FILE}" \
    -O z \
    -o "${VCF_OUT}"; then
    
    log_success "合并完成: ${VCF_OUT}"
    
    # 建立索引
    log_info "建立索引..."
    if bcftools index -t --threads "${INDEX_THREADS}" "${VCF_OUT}"; then
        log_success "索引创建完成"
    else
        log_error "索引创建失败"
        exit 1
    fi
else
    log_error "文件合并失败"
    exit 1
fi

# ==============================================================================
# --- [格式验证和转换] ---
# ==============================================================================

log_info "=========================================="
log_info " 验证输出格式"
log_info "=========================================="

# 验证 VCF 格式
log_info "验证 VCF 文件格式..."
if bcftools view -h "${VCF_OUT}" | head -1 | grep -q "^##fileformat=VCF"; then
    log_success "VCF 格式验证通过"
else
    log_error "输出文件不是有效的 VCF 格式"
    exit 1
fi

# 可选：生成未压缩的 VCF 文件
if [ "${GENERATE_UNCOMPRESSED:-false}" = true ]; then
    log_info "生成未压缩的 VCF 文件..."
    UNCOMPRESSED_VCF="${OUTPUT_DIR}/merged_output.vcf"
    bcftools view "${VCF_OUT}" > "${UNCOMPRESSED_VCF}"
    log_success "未压缩 VCF: ${UNCOMPRESSED_VCF}"
fi

# ==============================================================================
# --- [输出统计信息] ---
# ==============================================================================

log_info "=========================================="
log_info " 生成统计信息"
log_info "=========================================="

VARIANT_COUNT=$(bcftools view -H "${VCF_OUT}" | wc -l)
SAMPLE_COUNT=$(bcftools query -l "${VCF_OUT}" | wc -l)
FILE_SIZE=$(du -h "${VCF_OUT}" | cut -f1)

cat << EOF

#################################################
#              任务完成！                        
#################################################

输出文件: 
  - 压缩 VCF: ${VCF_OUT}
  - 索引文件: ${VCF_OUT}.tbi
$([ "${GENERATE_UNCOMPRESSED}" = true ] && echo "  - 未压缩 VCF: ${OUTPUT_DIR}/merged_output.vcf")

文件大小: ${FILE_SIZE}
样本数量: ${SAMPLE_COUNT}
变异数量: ${VARIANT_COUNT}

中间文件: ${TMP_OUT_DIR}
日志目录: ${LOG_DIR}

$([ "${KEEP_INTERMEDIATE}" = false ] && echo "提示: 确认结果无误后，可手动删除中间文件")

说明:
  - VCF.GZ 是压缩的标准 VCF 格式，推荐使用
  - 可用 bcftools view 或 zcat 查看内容
  - 如需未压缩版本，使用 -u 参数重新运行

#################################################
EOF

log_success "全部任务完成！"