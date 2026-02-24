#!/bin/bash
set -eo pipefail  # 管道错误退出，但允许未定义变量（兼容性更好）

# ================= 配置区域 =================
REF_GENOME="/share/org/YZWL/yzwl_lixg/project/21.野生大豆群体/01.data/genome/genome.fa"
GVCF_DIR="/share/org/YZWL/yzwl_lixg/project/21.野生大豆群体/02.mapping/vcf"
OUTPUT_DIR="/share/org/YZWL/yzwl_lixg/project/21.野生大豆群体/03.glnexus_joint"
QUARANTINE_DIR="${OUTPUT_DIR}/quarantine_files"
LOG_FILE="${OUTPUT_DIR}/file_check.log"
GLNEXUS_CONFIG="gatk"
THREADS=88
MIN_FILE_SIZE=1000  # 最小文件大小（字节），小于此值视为异常
# ===========================================

# 颜色输出
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# ================= 工具函数 =================

log_info() {
    echo -e "${GREEN}[INFO]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1" | tee -a ${LOG_FILE}
}

log_warn() {
    echo -e "${YELLOW}[WARN]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1" | tee -a ${LOG_FILE}
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1" | tee -a ${LOG_FILE}
}

# 检查必需的软件工具
check_dependencies() {
    log_info "检查依赖软件..."
    local missing_tools=()
    
    for tool in bcftools samtools bgzip glnexus_cli; do
        if ! command -v $tool &> /dev/null; then
            missing_tools+=($tool)
        fi
    done
    
    if [ ${#missing_tools[@]} -gt 0 ]; then
        log_error "缺少必需工具: ${missing_tools[*]}"
        log_error "请先安装这些工具后再运行脚本"
        exit 1
    fi
    log_info "所有依赖工具已就绪 ✓"
}

# 检查参考基因组文件
check_reference() {
    log_info "检查参考基因组..."
    if [ ! -f "${REF_GENOME}" ]; then
        log_error "参考基因组不存在: ${REF_GENOME}"
        exit 1
    fi
    
    if [ ! -f "${REF_GENOME}.fai" ]; then
        log_info "正在为参考基因组建立索引..."
        samtools faidx ${REF_GENOME}
    fi
    log_info "参考基因组检查完成 ✓"
}

# 初始化目录结构
initialize_directories() {
    log_info "初始化输出目录..."
    mkdir -p ${OUTPUT_DIR}
    mkdir -p ${QUARANTINE_DIR}
    
    # 创建时间戳备份目录
    BACKUP_DIR="${OUTPUT_DIR}/backup_$(date '+%Y%m%d_%H%M%S')"
    mkdir -p ${BACKUP_DIR}
    
    log_info "输出目录: ${OUTPUT_DIR}"
    log_info "隔离目录: ${QUARANTINE_DIR}"
    log_info "备份目录: ${BACKUP_DIR}"
}

# ================= 核心功能函数 =================

# 检查并修复单个 gVCF 文件
check_and_fix_vcf() {
    local vcf_file=$1
    local vcf_name=$(basename "$vcf_file")
    local tbi_file="${vcf_file}.tbi"
    
    # 使用 trap 捕获所有错误，防止静默退出
    (
        # 1. 检查文件存在性
        if [ ! -f "$vcf_file" ]; then
            echo "[ERROR] 文件不存在: $vcf_name" >> ${LOG_FILE}
            return 1
        fi

        # 2. 检查文件大小（排除空文件或损坏文件）
        local file_size
        file_size=$(stat -c%s "$vcf_file" 2>/dev/null || stat -f%z "$vcf_file" 2>/dev/null || echo "0")
        if [ "$file_size" -lt "$MIN_FILE_SIZE" ]; then
            echo "[WARN] 文件过小 (${file_size} bytes): $vcf_name -> 移至隔离区" >> ${LOG_FILE}
            mv "$vcf_file" "${QUARANTINE_DIR}/" 2>/dev/null || true
            [ -f "$tbi_file" ] && mv "$tbi_file" "${QUARANTINE_DIR}/" 2>/dev/null || true
            return 1
        fi

        # 3. 检查是否为有效的 gzip 文件
        if ! gzip -t "$vcf_file" 2>/dev/null; then
            echo "[ERROR] GZIP格式损坏: $vcf_name -> 移至隔离区" >> ${LOG_FILE}
            mv "$vcf_file" "${QUARANTINE_DIR}/" 2>/dev/null || true
            [ -f "$tbi_file" ] && mv "$tbi_file" "${QUARANTINE_DIR}/" 2>/dev/null || true
            return 1
        fi

        # 4. 检查并重建索引（如果需要）
        if [ ! -f "$tbi_file" ] || [ "$vcf_file" -nt "$tbi_file" ]; then
            echo "[WARN] 索引过期或缺失，正在重建: $vcf_name" >> ${LOG_FILE}
            if ! bcftools index -t "$vcf_file" 2>>"${LOG_FILE}"; then
                echo "[ERROR] 索引重建失败: $vcf_name -> 移至隔离区" >> ${LOG_FILE}
                mv "$vcf_file" "${QUARANTINE_DIR}/" 2>/dev/null || true
                [ -f "$tbi_file" ] && mv "$tbi_file" "${QUARANTINE_DIR}/" 2>/dev/null || true
                return 1
            fi
        fi

        # 5. 快速完整性检查
        if ! bcftools index -n "$vcf_file" >/dev/null 2>&1; then
            echo "[ERROR] 文件结构损坏（索引无法读取）: $vcf_name -> 移至隔离区" >> ${LOG_FILE}
            mv "$vcf_file" "$tbi_file" "${QUARANTINE_DIR}/" 2>/dev/null || true
            return 1
        fi
        
        # 6. 检查VCF头部完整性（快速检查前100行）
        if ! bcftools view -h "$vcf_file" 2>/dev/null | head -100 >/dev/null 2>&1; then
            echo "[ERROR] VCF头部损坏: $vcf_name -> 移至隔离区" >> ${LOG_FILE}
            mv "$vcf_file" "$tbi_file" "${QUARANTINE_DIR}/" 2>/dev/null || true
            return 1
        fi
        
        return 0
    ) || return 1  # 子shell退出码传递
}

# 批量检查所有 gVCF 文件
batch_check_gvcfs() {
    log_info "=========================================="
    log_info "步骤 1: 批量检查 gVCF 文件"
    log_info "=========================================="
    
    # 统计文件总数
    local total
    total=$(find ${GVCF_DIR} -name "*.g.vcf.gz" 2>/dev/null | wc -l)
    if [ "$total" -eq 0 ]; then
        log_error "未找到任何 .g.vcf.gz 文件在目录: ${GVCF_DIR}"
        exit 1
    fi
    
    log_info "发现 $total 个 gVCF 文件"
    log_info "开始检查和修复..."
    
    local count=0
    local success=0
    local failed=0
    local start_time=$(date +%s)
    
    # 创建临时文件列表
    local temp_file_list="${OUTPUT_DIR}/.temp_vcf_list.txt"
    find ${GVCF_DIR} -name "*.g.vcf.gz" > "$temp_file_list" 2>/dev/null
    
    while IFS= read -r vcf; do
        ((count++)) || true
        
        # 显示进度（更频繁的更新，每10个文件）
        if (( count % 10 == 0 )) || (( count == total )); then
            local percent=$((count * 100 / total))
            local elapsed=$(($(date +%s) - start_time))
            local eta=0
            if [ $count -gt 0 ]; then
                eta=$(awk "BEGIN {printf \"%.0f\", $elapsed / $count * ($total - $count)}")
            fi
            echo -ne "\r进度: ${count}/${total} (${percent}%) | 成功: ${success} | 失败: ${failed} | 预计剩余: ${eta}s    "
        fi
        
        # 检查文件，即使失败也继续
        if check_and_fix_vcf "$vcf"; then
            ((success++)) || true
        else
            ((failed++)) || true
            # 记录失败的文件名
            echo "$vcf" >> "${OUTPUT_DIR}/failed_files.txt"
        fi
    done < "$temp_file_list"
    
    # 清理临时文件
    rm -f "$temp_file_list"
    
    echo ""  # 换行
    log_info "检查完成！"
    log_info "总文件数: $total | 成功: $success | 失败: $failed"
    
    if [ $failed -gt 0 ]; then
        log_warn "有 $failed 个文件存在问题，已移至: ${QUARANTINE_DIR}"
        log_warn "失败文件列表: ${OUTPUT_DIR}/failed_files.txt"
    fi
    
    if [ $success -eq 0 ]; then
        log_error "没有有效的 gVCF 文件可供合并！"
        exit 1
    fi
}

# 准备 BED 文件
prepare_bed_file() {
    log_info "=========================================="
    log_info "步骤 2: 准备 BED 文件"
    log_info "=========================================="
    
    BED_FILE="${OUTPUT_DIR}/genome.bed"
    
    awk -v OFS='\t' '{print $1, 0, $2}' ${REF_GENOME}.fai > ${BED_FILE}
    
    local chr_count=$(wc -l < ${BED_FILE})
    log_info "BED文件已生成: ${BED_FILE}"
    log_info "染色体/scaffold数量: $chr_count"
}

# 清理旧的GLnexus缓存
cleanup_old_cache() {
    log_info "=========================================="
    log_info "步骤 3: 清理旧的 GLnexus 缓存"
    log_info "=========================================="
    
    if [ -d "GLnexus.DB" ]; then
        log_warn "发现旧的 GLnexus.DB，正在移动到备份目录..."
        mv GLnexus.DB ${BACKUP_DIR}/GLnexus.DB_old 2>/dev/null || true
    fi
    
    BCF_OUT="${OUTPUT_DIR}/glnexus_output.bcf"
    if [ -f "${BCF_OUT}" ]; then
        log_warn "发现旧的输出文件，正在备份..."
        mv ${BCF_OUT} ${BACKUP_DIR}/glnexus_output.bcf.old 2>/dev/null || true
    fi
    
    log_info "缓存清理完成 ✓"
}

# 运行 GLnexus
run_glnexus() {
    log_info "=========================================="
    log_info "步骤 4: 运行 GLnexus 联合检测"
    log_info "=========================================="
    
    local valid_count=$(find ${GVCF_DIR} -name "*.g.vcf.gz" | wc -l)
    log_info "参与合并的文件数: $valid_count"
    log_info "配置模式: ${GLNEXUS_CONFIG}"
    log_info "线程数: ${THREADS}"
    log_info "开始时间: $(date '+%Y-%m-%d %H:%M:%S')"
    
    BCF_OUT="${OUTPUT_DIR}/glnexus_output.bcf"
    
    # 运行 GLnexus 并捕获错误
    if glnexus_cli \
        --config ${GLNEXUS_CONFIG} \
        --bed ${BED_FILE} \
        ${GVCF_DIR}/*.g.vcf.gz \
        --threads ${THREADS} \
        > ${BCF_OUT} 2>&1 | tee -a ${LOG_FILE}; then
        
        log_info "GLnexus 运行成功 ✓"
        log_info "结束时间: $(date '+%Y-%m-%d %H:%M:%S')"
    else
        log_error "GLnexus 运行失败！"
        log_error "错误信息已记录到: ${LOG_FILE}"
        
        # 检查是否是负数PL错误
        if grep -q "negative PL entry" ${LOG_FILE}; then
            log_error "检测到 'negative PL entry' 错误"
            log_error "请检查日志找出问题文件，手动移至隔离区后重新运行"
        fi
        
        exit 1
    fi
}

# 后处理：转换为VCF并压缩
post_process() {
    log_info "=========================================="
    log_info "步骤 5: 后处理 - 转换与压缩"
    log_info "=========================================="
    
    VCF_OUT="${OUTPUT_DIR}/glnexus_output.vcf.gz"
    
    log_info "正在转换 BCF -> VCF.GZ..."
    bcftools view ${BCF_OUT} | bgzip -@ 4 -c > ${VCF_OUT}
    
    log_info "正在建立索引..."
    bcftools index -t ${VCF_OUT}
    
    # 生成统计信息
    log_info "正在生成统计信息..."
    bcftools stats ${VCF_OUT} > ${OUTPUT_DIR}/vcf_stats.txt
    
    local var_count=$(bcftools view -H ${VCF_OUT} | wc -l)
    local file_size=$(du -h ${VCF_OUT} | cut -f1)
    
    log_info "变异位点数: $var_count"
    log_info "文件大小: $file_size"
    log_info "输出文件: ${VCF_OUT}"
    log_info "统计信息: ${OUTPUT_DIR}/vcf_stats.txt"
}

# 生成最终报告
generate_report() {
    log_info "=========================================="
    log_info "流程完成汇总"
    log_info "=========================================="
    
    REPORT_FILE="${OUTPUT_DIR}/pipeline_report.txt"
    
    cat > ${REPORT_FILE} <<EOF
GLnexus 联合变异检测流程报告
=========================================
运行时间: $(date '+%Y-%m-%d %H:%M:%S')
参考基因组: ${REF_GENOME}
输入目录: ${GVCF_DIR}
输出目录: ${OUTPUT_DIR}

文件统计:
- 输入 gVCF 文件数: $(find ${GVCF_DIR} -name "*.g.vcf.gz" | wc -l)
- 隔离文件数: $(find ${QUARANTINE_DIR} -name "*.g.vcf.gz" 2>/dev/null | wc -l)
- 输出 VCF: ${VCF_OUT}
- 变异位点数: $(bcftools view -H ${VCF_OUT} | wc -l)

配置参数:
- GLnexus 配置: ${GLNEXUS_CONFIG}
- 线程数: ${THREADS}

详细日志: ${LOG_FILE}
=========================================
EOF
    
    cat ${REPORT_FILE}
    log_info "报告已保存至: ${REPORT_FILE}"
}

# ================= 主流程 =================

main() {
    log_info "=========================================="
    log_info "GLnexus 联合变异检测流程启动"
    log_info "=========================================="
    
    # 捕获所有错误并记录
    trap 'handle_error $? $LINENO' ERR
    
    # 前置检查
    check_dependencies
    check_reference
    initialize_directories
    
    # 核心流程
    batch_check_gvcfs
    prepare_bed_file
    cleanup_old_cache
    run_glnexus
    post_process
    
    # 生成报告
    generate_report
    
    log_info "=========================================="
    log_info "所有步骤完成！ ✓"
    log_info "=========================================="
}

# 错误处理函数
handle_error() {
    local exit_code=$1
    local line_no=$2
    log_error "脚本在第 ${line_no} 行异常退出，退出码: ${exit_code}"
    log_error "请检查日志文件: ${LOG_FILE}"
    exit $exit_code
}

# 执行主流程
main "$@"