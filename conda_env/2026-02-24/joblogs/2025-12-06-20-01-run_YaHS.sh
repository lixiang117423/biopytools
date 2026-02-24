#!/bin/bash
# =============================================================================
#  🧬 YaHS 高速染色体挂载流程 - 优化增强版 (v4.0)
#  日期: 2025-12-06
#  优化: 模块化设计、智能资源管理、增强错误处理、断点续传
# =============================================================================

# =============================================================================
# --- 🔧 基础配置 ---
# =============================================================================
set -euo pipefail

# 颜色输出定义
readonly RED='\033[0;31m'
readonly GREEN='\033[0;32m'
readonly YELLOW='\033[1;33m'
readonly BLUE='\033[0;34m'
readonly MAGENTA='\033[0;35m'
readonly CYAN='\033[0;36m'
readonly NC='\033[0m'

# 全局变量
PIPELINE_START_TIME=$(date +%s)
LOG_FILE=""

# =============================================================================
# --- 📝 日志系统 ---
# =============================================================================
init_logging() {
    local log_dir="$1/logs"
    mkdir -p "${log_dir}"
    LOG_FILE="${log_dir}/pipeline_$(date +%Y%m%d_%H%M%S).log"
    exec > >(tee -a "${LOG_FILE}")
    exec 2>&1
}

log_msg() {
    local level="$1"
    local color="$2"
    local icon="$3"
    shift 3
    local msg="$*"
    echo -e "${color}${icon} [$(date '+%Y-%m-%d %H:%M:%S')] ${level}:${NC} ${msg}" | tee -a "${LOG_FILE}"
}

log_info() { log_msg "INFO" "${BLUE}" "ℹ️ " "$@"; }
log_success() { log_msg "SUCCESS" "${GREEN}" "✅" "$@"; }
log_warning() { log_msg "WARNING" "${YELLOW}" "⚠️ " "$@"; }
log_error() { log_msg "ERROR" "${RED}" "❌" "$@"; }
log_step() { log_msg "STEP" "${CYAN}" "🔹" "$@"; }

# 错误处理
error_handler() {
    local line_no=$1
    local exit_code=$2
    log_error "脚本在第 ${line_no} 行失败 (退出码: ${exit_code})"
    log_error "请检查日志文件: ${LOG_FILE}"
    exit "${exit_code}"
}

trap 'error_handler ${LINENO} $?' ERR

# =============================================================================
# --- 💻 配置管理 ---
# =============================================================================
declare -A CONFIG

load_config() {
    # 1. 📂 路径设置
    CONFIG[WORK_DIR]="/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/04.YaHS"
    CONFIG[REF_FA]="${CONFIG[WORK_DIR]}/OV53_1.primary.fa"
    CONFIG[R1_FQ]="${CONFIG[WORK_DIR]}/fastq/OV53_1-hic_R1.fastq.gz"
    CONFIG[R2_FQ]="${CONFIG[WORK_DIR]}/fastq/OV53_1-hic_R2.fastq.gz"
    
    # 2. ⚙️ 软件路径
    CONFIG[JAVA_CMD]="java"
    CONFIG[JUICER_JAR]="/share/org/YZWL/yzwl_lixg/software/juicer/scripts/common/juicer_tools.jar"
    CONFIG[YAHS_JUICER]="/share/org/YZWL/yzwl_lixg/miniforge3/envs/yahs_v.1.2.2/bin/juicer"
    
    # 3. 🔧 参数配置
    CONFIG[ENZYME_SEQ]="GATC"
    CONFIG[MIN_LEN]=10000
    CONFIG[MIN_MAPQ]=30
    
    # 4. 💾 资源配置（自动检测系统资源）
    local total_mem_gb=$(free -g | awk '/^Mem:/{print $2}')
    local cpu_cores=$(nproc)
    
    CONFIG[THREADS]=${THREADS:-$((cpu_cores - 2))}  # 预留2核
    CONFIG[JAVA_RAM]=${JAVA_RAM:-"$((total_mem_gb * 6 / 10))G"}  # 60%内存给Java
    CONFIG[SORT_RAM]=${SORT_RAM:-"$((total_mem_gb * 2 / 10))G"}  # 20%内存给Sort
    CONFIG[SAM_MEM]="4G"
    
    # 5. 环境变量
    export PATH="/share/org/YZWL/yzwl_lixg/miniforge3/envs/yahs_v.1.2.2/bin:${PATH}"
}

print_config() {
    cat << EOF

╔══════════════════════════════════════════════════════════════╗
║                  📋 流程配置信息                              ║
╠══════════════════════════════════════════════════════════════╣
║ 工作目录: ${CONFIG[WORK_DIR]}
║ 参考基因组: ${CONFIG[REF_FA]}
║ Hi-C R1: ${CONFIG[R1_FQ]}
║ Hi-C R2: ${CONFIG[R2_FQ]}
║ 
║ 线程数: ${CONFIG[THREADS]}
║ Java内存: ${CONFIG[JAVA_RAM]}
║ Sort内存: ${CONFIG[SORT_RAM]}
║ 限制性酶: ${CONFIG[ENZYME_SEQ]}
║ 最小长度: ${CONFIG[MIN_LEN]} bp
║ 最小MAPQ: ${CONFIG[MIN_MAPQ]}
╚══════════════════════════════════════════════════════════════╝

EOF
}

# =============================================================================
# --- 🔍 环境检查 ---
# =============================================================================
check_environment() {
    log_step "开始环境检查"
    
    # 检查工作目录
    if [[ ! -d "${CONFIG[WORK_DIR]}" ]]; then
        log_error "工作目录不存在: ${CONFIG[WORK_DIR]}"
        return 1
    fi
    
    cd "${CONFIG[WORK_DIR]}" || return 1
    mkdir -p logs tmp_files results
    
    # 检查 Java
    log_info "检查 Java 版本:"
    if ! ${CONFIG[JAVA_CMD]} -version 2>&1 | head -n 3 | tee logs/java_version.log; then
        log_error "Java 未找到或无法执行"
        return 1
    fi
    
    # 检查必需工具
    local required_tools=("bwa" "samtools" "yahs" "awk" "sort" "gzip")
    for tool in "${required_tools[@]}"; do
        if ! command -v "${tool}" &> /dev/null; then
            log_error "未找到必需工具: ${tool}"
            return 1
        fi
    done
    log_info "✓ 所有必需工具已安装"
    
    # 检查 YaHS Juicer
    if [[ ! -x "${CONFIG[YAHS_JUICER]}" ]]; then
        log_error "YaHS juicer工具不可执行: ${CONFIG[YAHS_JUICER]}"
        return 1
    fi
    
    # 检查 Juicer JAR
    if [[ ! -f "${CONFIG[JUICER_JAR]}" ]]; then
        log_error "Juicer JAR文件未找到: ${CONFIG[JUICER_JAR]}"
        return 1
    fi
    
    # 检查输入文件
    local files=("${CONFIG[REF_FA]}" "${CONFIG[R1_FQ]}" "${CONFIG[R2_FQ]}")
    for file in "${files[@]}"; do
        if [[ ! -f "${file}" ]] || [[ ! -s "${file}" ]]; then
            log_error "输入文件不存在或为空: ${file}"
            return 1
        fi
    done
    log_info "✓ 所有输入文件验证通过"
    
    # 检查磁盘空间（至少需要输入文件大小的5倍）
    local required_space=$(( $(du -sb "${CONFIG[R1_FQ]}" "${CONFIG[R2_FQ]}" | awk '{sum+=$1} END {print sum}') * 5 / 1024 / 1024 / 1024 ))
    local available_space=$(df -BG "${CONFIG[WORK_DIR]}" | tail -1 | awk '{print $4}' | sed 's/G//')
    
    if (( available_space < required_space )); then
        log_warning "磁盘空间可能不足 (需要约 ${required_space}G，可用 ${available_space}G)"
    else
        log_info "✓ 磁盘空间充足 (可用 ${available_space}G)"
    fi
    
    log_success "环境检查完成"
    return 0
}

# =============================================================================
# --- 🔧 工具函数 ---
# =============================================================================

# 检查文件是否存在且有效
check_file() {
    local file="$1"
    local min_size="${2:-100}"  # 默认最小100字节
    
    if [[ -f "${file}" ]] && [[ $(stat -c%s "${file}") -gt ${min_size} ]]; then
        return 0
    fi
    return 1
}

# 安全删除文件
safe_remove() {
    local file="$1"
    if [[ -f "${file}" ]]; then
        rm -f "${file}"
        log_info "已删除临时文件: $(basename ${file})"
    fi
}

# 计算运行时间
calculate_runtime() {
    local start_time=$1
    local end_time=$(date +%s)
    local duration=$((end_time - start_time))
    
    local hours=$((duration / 3600))
    local minutes=$(((duration % 3600) / 60))
    local seconds=$((duration % 60))
    
    printf "%02d:%02d:%02d" ${hours} ${minutes} ${seconds}
}

# =============================================================================
# --- 步骤 1: 基因组索引 ---
# =============================================================================
step_indexing() {
    echo ""
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    log_step "步骤 1/6: 构建基因组索引"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    
    local step_start=$(date +%s)
    local ref_fa="${CONFIG[REF_FA]}"
    
    # BWA 索引
    if [[ ! -f "${ref_fa}.bwt" ]]; then
        log_info "构建 BWA 索引..."
        bwa index "${ref_fa}" 2>&1 | tee logs/bwa_index.log
        log_success "BWA 索引构建完成"
    else
        log_info "✓ BWA 索引已存在，跳过"
    fi
    
    # SAMtools 索引
    if [[ ! -f "${ref_fa}.fai" ]]; then
        log_info "构建 SAMtools 索引..."
        samtools faidx "${ref_fa}" 2>&1 | tee logs/samtools_faidx.log
        log_success "SAMtools 索引构建完成"
    else
        log_info "✓ SAMtools 索引已存在，跳过"
    fi
    
    log_success "步骤 1 完成 [用时: $(calculate_runtime ${step_start})]"
}

# =============================================================================
# --- 步骤 2: Hi-C 比对 ---
# =============================================================================
step_mapping() {
    echo ""
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    log_step "步骤 2/6: Hi-C 数据比对与处理"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    
    local step_start=$(date +%s)
    local final_bam="results/aligned_sorted_dedup.bam"
    
    # 检查是否已完成
    if check_file "${final_bam}" 1000000 && check_file "${final_bam}.bai"; then
        log_info "✓ BAM 文件已存在，跳过比对步骤"
        return 0
    fi
    
    # 清理并创建临时目录
    rm -rf tmp_files/nsort tmp_files/sort
    mkdir -p tmp_files/{nsort,sort}
    
    # BWA 比对
    log_info "[1/6] BWA MEM 比对 (参数: -5SP)..."
    bwa mem -5SP -t "${CONFIG[THREADS]}" "${CONFIG[REF_FA]}" \
        "${CONFIG[R1_FQ]}" "${CONFIG[R2_FQ]}" 2> logs/bwa_mem.log | \
    samtools view -@ "${CONFIG[THREADS]}" -bS - > tmp_files/aligned.bam
    
    if ! check_file "tmp_files/aligned.bam" 1000000; then
        log_error "BWA 比对失败"
        return 1
    fi
    log_success "比对完成"
    
    # 按名称排序
    log_info "[2/6] 按 read name 排序..."
    samtools sort -n -@ "${CONFIG[THREADS]}" -m "${CONFIG[SAM_MEM]}" \
        -T tmp_files/nsort/split \
        -o tmp_files/aligned_nsorted.bam tmp_files/aligned.bam \
        2> logs/samtools_nsort.log
    safe_remove "tmp_files/aligned.bam"
    
    # 修复 mate pair
    log_info "[3/6] 修复 mate pair 信息..."
    samtools fixmate -m -@ "${CONFIG[THREADS]}" \
        tmp_files/aligned_nsorted.bam tmp_files/aligned_fixmate.bam \
        2> logs/samtools_fixmate.log
    safe_remove "tmp_files/aligned_nsorted.bam"
    
    # 按坐标排序
    log_info "[4/6] 按坐标排序..."
    samtools sort -@ "${CONFIG[THREADS]}" -m "${CONFIG[SAM_MEM]}" \
        -T tmp_files/sort/split \
        -o tmp_files/aligned_sorted.bam tmp_files/aligned_fixmate.bam \
        2> logs/samtools_sort.log
    safe_remove "tmp_files/aligned_fixmate.bam"
    
    # 标记重复
    log_info "[5/6] 标记并移除 PCR 重复..."
    samtools markdup -r -@ "${CONFIG[THREADS]}" \
        tmp_files/aligned_sorted.bam "${final_bam}" \
        2> logs/samtools_markdup.log
    safe_remove "tmp_files/aligned_sorted.bam"
    
    # 构建索引
    log_info "[6/6] 构建 BAM 索引..."
    samtools index -@ "${CONFIG[THREADS]}" "${final_bam}"
    
    # 统计比对结果
    log_info "比对统计信息:"
    samtools flagstat -@ "${CONFIG[THREADS]}" "${final_bam}" | tee logs/flagstat.txt
    
    log_success "步骤 2 完成 [用时: $(calculate_runtime ${step_start})]"
}

# =============================================================================
# --- 步骤 3: YaHS 组装 ---
# =============================================================================
step_scaffolding() {
    echo ""
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    log_step "步骤 3/6: YaHS 染色体挂载"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    
    local step_start=$(date +%s)
    local out_prefix="results/yahs_out"
    local final_fa="${out_prefix}_scaffolds_final.fa"
    
    if check_file "${final_fa}" 1000; then
        log_info "✓ Scaffold 文件已存在，跳过 YaHS"
        return 0
    fi
    
    log_info "运行 YaHS 组装..."
    log_info "参数: -e ${CONFIG[ENZYME_SEQ]} -q ${CONFIG[MIN_MAPQ]} -l ${CONFIG[MIN_LEN]}"
    
    yahs -e "${CONFIG[ENZYME_SEQ]}" \
         -q "${CONFIG[MIN_MAPQ]}" \
         -l "${CONFIG[MIN_LEN]}" \
         -o "${out_prefix}" \
         --no-contig-ec \
         "${CONFIG[REF_FA]}" \
         "results/aligned_sorted_dedup.bam" \
         2>&1 | tee logs/yahs.log
    
    if ! check_file "${out_prefix}_scaffolds_final.agp"; then
        log_error "YaHS 运行失败，未生成 AGP 文件"
        return 1
    fi
    
    log_success "步骤 3 完成 [用时: $(calculate_runtime ${step_start})]"
}

# =============================================================================
# --- 步骤 4: 生成标准 .hic 文件 ---
# =============================================================================
step_hic_standard() {
    echo ""
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    log_step "步骤 4/6: 生成标准 Hi-C 热图"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    
    local step_start=$(date +%s)
    local out_prefix="results/yahs_out"
    local hic_file="results/yahs_out_final.hic"
    
    if check_file "${hic_file}" 100000; then
        log_info "✓ 标准 .hic 文件已存在，跳过"
        return 0
    fi
    
    # 检查依赖文件
    if ! check_file "${out_prefix}.bin"; then
        log_error "YaHS .bin 文件未找到"
        return 1
    fi
    
    # 生成并排序比对链接
    log_info "[1/3] 生成并排序比对链接..."
    "${CONFIG[YAHS_JUICER]}" pre "${out_prefix}.bin" \
        "${out_prefix}_scaffolds_final.agp" \
        "${CONFIG[REF_FA]}.fai" 2> logs/juicer_pre.log | \
    sort -k2,2d -k6,6d -k3,3n -k7,7n \
        -T tmp_files/ \
        --parallel="${CONFIG[THREADS]}" \
        -S"${CONFIG[SORT_RAM]}" | \
    awk 'NF' > tmp_files/alignments_sorted.txt
    
    if ! check_file "tmp_files/alignments_sorted.txt"; then
        log_error "比对链接生成失败"
        cat logs/juicer_pre.log
        return 1
    fi
    
    # 准备染色体大小
    log_info "[2/3] 准备染色体大小文件..."
    samtools faidx "${out_prefix}_scaffolds_final.fa"
    cut -f1,2 "${out_prefix}_scaffolds_final.fa.fai" > results/chrom.sizes.final
    
    # 生成 .hic
    log_info "[3/3] 运行 Juicer Tools..."
    ${CONFIG[JAVA_CMD]} -Xmx"${CONFIG[JAVA_RAM]}" -Xms8G \
        -jar "${CONFIG[JUICER_JAR]}" pre \
        tmp_files/alignments_sorted.txt \
        results/out.hic.part \
        results/chrom.sizes.final \
        2>&1 | tee logs/juicer_tools.log
    
    if check_file "results/out.hic.part"; then
        mv results/out.hic.part "${hic_file}"
        log_success "标准 .hic 文件生成成功"
        safe_remove "tmp_files/alignments_sorted.txt"
    else
        log_error ".hic 生成失败，请检查日志"
        return 1
    fi
    
    log_success "步骤 4 完成 [用时: $(calculate_runtime ${step_start})]"
}

# =============================================================================
# --- 步骤 5: 生成 JBAT 文件 ---
# =============================================================================
step_jbat() {
    echo ""
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    log_step "步骤 5/6: 生成 JBAT 文件（手动纠错用）"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    
    local step_start=$(date +%s)
    local out_prefix="results/yahs_out"
    local jbat_hic="results/out_JBAT.hic"
    
    if check_file "${jbat_hic}" 100000; then
        log_info "✓ JBAT 文件已存在，跳过"
        return 0
    fi
    
    # 生成 JBAT 文本
    log_info "[1/4] 生成 JBAT 文本格式..."
    "${CONFIG[YAHS_JUICER]}" pre -a -o results/out_JBAT \
        "${out_prefix}.bin" \
        "${out_prefix}_scaffolds_final.agp" \
        "${CONFIG[REF_FA]}.fai" \
        > logs/out_JBAT.log 2>&1
    
    if ! check_file "results/out_JBAT.txt"; then
        log_error "JBAT 文本生成失败"
        return 1
    fi
    
    # 提取 assembly 大小
    log_info "[2/4] 提取 assembly 大小信息..."
    if grep -q "PRE_C_SIZE" logs/out_JBAT.log; then
        grep "PRE_C_SIZE" logs/out_JBAT.log | \
            awk '{print $2" "$3}' > results/jbat_chrom_sizes.txt
    else
        log_warning "未找到 PRE_C_SIZE，使用备用方法"
        local total_bp=$(awk '{s+=$2} END {print s}' results/chrom.sizes.final)
        echo "assembly ${total_bp}" > results/jbat_chrom_sizes.txt
    fi
    
    # 关键：排序 JBAT 文件
    log_info "[3/4] 对 JBAT 文件进行排序（防止 OOM）..."
    sort -k2,2d -k6,6d -k3,3n -k7,7n \
        --parallel="${CONFIG[THREADS]}" \
        -S"${CONFIG[SORT_RAM]}" \
        -T tmp_files/ \
        results/out_JBAT.txt > results/out_JBAT_sorted.txt
    
    # 生成 .hic
    log_info "[4/4] 生成 JBAT .hic 文件..."
    ${CONFIG[JAVA_CMD]} -Xmx"${CONFIG[JAVA_RAM]}" -Xms8G \
        -jar "${CONFIG[JUICER_JAR]}" pre \
        results/out_JBAT_sorted.txt \
        results/out_JBAT.hic.part \
        results/jbat_chrom_sizes.txt \
        2>&1 | tee logs/juicer_jbat.log
    
    if check_file "results/out_JBAT.hic.part"; then
        mv results/out_JBAT.hic.part "${jbat_hic}"
        log_success "JBAT 文件生成成功"
        
        # 清理大文件（可选）
        if [[ "${KEEP_TEMP:-false}" != "true" ]]; then
            safe_remove "results/out_JBAT.txt"
            safe_remove "results/out_JBAT_sorted.txt"
        fi
    else
        log_error "JBAT .hic 生成失败"
        return 1
    fi
    
    log_success "步骤 5 完成 [用时: $(calculate_runtime ${step_start})]"
}

# =============================================================================
# --- 步骤 6: 质量评估 ---
# =============================================================================
step_assessment() {
    echo ""
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    log_step "步骤 6/6: 组装质量评估"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    
    local step_start=$(date +%s)
    local scaffold_fa="results/yahs_out_scaffolds_final.fa"
    
    if ! check_file "${scaffold_fa}"; then
        log_error "Scaffold 文件未找到"
        return 1
    fi
    
    # 提取长度
    log_info "计算统计指标..."
    awk '/^>/ {if (seq) print length(seq); seq=""; next} 
         {seq=seq$0} 
         END {if (seq) print length(seq)}' "${scaffold_fa}" | \
    sort -rn > tmp_files/scaffold_lengths.txt
    
    # 基础统计
    local total_length=$(awk '{sum += $1} END {print sum}' tmp_files/scaffold_lengths.txt)
    local n_scaffolds=$(wc -l < tmp_files/scaffold_lengths.txt)
    local max_scaffold=$(head -1 tmp_files/scaffold_lengths.txt)
    
    # N50/N90
    awk -v total="${total_length}" '
    BEGIN { n50=0; n90=0; cum=0 }
    {
        cum += $1
        if (cum >= total*0.5 && n50==0) n50=$1
        if (cum >= total*0.9 && n90==0) {n90=$1; exit}
    }
    END { 
        print "N50\t" n50
        print "N90\t" n90
    }' tmp_files/scaffold_lengths.txt > results/assembly_metrics.txt
    
    # 输出报告
    echo ""
    echo "╔════════════════════════════════════════════════════════════╗"
    echo "║           📊 组装质量统计报告                              ║"
    echo "╠════════════════════════════════════════════════════════════╣"
    printf "║ %-30s %28s ║\n" "Scaffold 总数:" "${n_scaffolds}"
    printf "║ %-30s %28s ║\n" "总长度:" "$(printf "%'d" ${total_length}) bp"
    printf "║ %-30s %28s ║\n" "最长 Scaffold:" "$(printf "%'d" ${max_scaffold}) bp"
    
    while IFS=$'\t' read -r metric value; do
        printf "║ %-30s %28s ║\n" "${metric}:" "$(printf "%'d" ${value}) bp"
    done < results/assembly_metrics.txt
    
    echo "╚════════════════════════════════════════════════════════════╝"
    
    log_success "步骤 6 完成 [用时: $(calculate_runtime ${step_start})]"
}

# =============================================================================
# --- 主流程 ---
# =============================================================================
main() {
    echo ""
    echo "╔══════════════════════════════════════════════════════════════╗"
    echo "║      🧬 YaHS 高速染色体挂载流程 - 优化版 v4.0              ║"
    echo "╚══════════════════════════════════════════════════════════════╝"
    echo ""
    
    # 初始化
    load_config
    init_logging "${CONFIG[WORK_DIR]}"
    
    log_info "流程启动于: $(date '+%Y-%m-%d %H:%M:%S')"
    log_info "运行节点: $(hostname)"
    log_info "当前用户: $(whoami)"
    log_info "日志文件: ${LOG_FILE}"
    
    print_config
    
    # 环境检查
    if ! check_environment; then
        log_error "环境检查失败，退出"
        exit 1
    fi
    
    # 执行各步骤
    step_indexing || exit 1
    step_mapping || exit 1
    step_scaffolding || exit 1
    step_hic_standard || exit 1
    step_jbat || exit 1
    step_assessment || exit 1
    
    # 最终报告
    local total_time=$(calculate_runtime ${PIPELINE_START_TIME})
    
    echo ""
    echo "╔════════════════════════════════════════════════════════════╗"
    echo "║                  🎉 流程执行完成                             ║"
    echo "╠════════════════════════════════════════════════════════════╣"
    printf "║ %-30s %28s ║\n" "总运行时间:" "${total_time}"
    printf "║ %-30s %28s ║\n" "工作目录:" "$(basename ${CONFIG[WORK_DIR]})"
    printf "║ %-30s %28s ║\n" "日志文件:" "$(basename ${LOG_FILE})"
    echo "║                                                            ║"
    echo "║ 主要输出文件:                                                ║"
    echo "║   • results/yahs_out_scaffolds_final.fa                    ║"
    echo "║   • results/yahs_out_final.hic                             ║"
    echo "║   • results/out_JBAT.hic                                   ║"
    echo "║   • results/out_JBAT.assembly                              ║"
    echo "║   • results/assembly_metrics.txt                           ║"
    echo "╚════════════════════════════════════════════════════════════╝"
    echo ""
    
    log_success "全流程执行完毕！"
}

# =============================================================================
# --- 脚本入口 ---
# =============================================================================
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    main "$@"
fi