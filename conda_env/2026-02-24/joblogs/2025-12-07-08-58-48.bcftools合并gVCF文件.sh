#!/bin/bash

# set -euo pipefail  # 严格错误处理

# ================= 配置区域 =================
# 1. 输入包含 g.vcf.gz 的文件夹路径
IN_DIR="/share/org/YZWL/yzwl_lixg/project/21.野生大豆群体/03.glnexus_joint/vcf"

# 2. 输出文件夹路径
OUT_DIR="/share/org/YZWL/yzwl_lixg/project/21.野生大豆群体/03.第一批720个样品/03.bcftools_joint"

# 3. 参考基因组路径 (***必须修改此处***)
# 请填写你的参考基因组 .fa 或 .fasta 文件的绝对路径
# 例如: /share/database/soybean/Wm82.a2.v1.genome.fa
REF_GENOME="/share/org/YZWL/yzwl_lixg/project/21.野生大豆群体/01.data/genome/genome.fa"

# 4. 输出文件名
OUT_FILE="${OUT_DIR}/all_samples_joint_called.vcf.gz"

# 5. 线程分配
THREADS_MERGE=40
THREADS_CALL=40
THREADS_INDEX=20

# 6. 日志和临时文件配置
LOG_FILE="${OUT_DIR}/joint_calling_$(date +%Y%m%d_%H%M%S).log"
TMP_DIR="${OUT_DIR}/tmp"
LIST_FILE="${OUT_DIR}/vcf_file_list.txt"

# 7. 高级选项
MIN_OPEN_FILES=2000          # 最小文件句柄限制
BACKUP_EXISTING=true         # 是否备份已存在的输出文件
CHECK_VCF_INTEGRITY=true     # 是否检查输入VCF完整性
CALL_OUTPUT_ALL_SITES=false  # 是否输出所有位点(包括0/0)
PLOIDY=2                     # 倍性设置 (大豆是二倍体,设为2)
USE_NAIVE_MERGE=false        # 是否使用简单合并模式(不进行joint calling)
AUTO_REMOVE_QS=true          # 自动移除有问题的QS标签
STRIP_PROBLEMATIC_TAGS=true  # 移除其他可能有问题的标签

# ===========================================

# 颜色输出
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# 日志函数
log() {
    local level=$1
    shift
    local color=$NC
    case $level in
        ERROR) color=$RED ;;
        WARN)  color=$YELLOW ;;
        INFO)  color=$GREEN ;;
        DEBUG) color=$BLUE ;;
    esac
    echo -e "${color}[$(date '+%Y-%m-%d %H:%M:%S')] [$level]${NC} $*" | tee -a "$LOG_FILE"
}

error_exit() {
    log "ERROR" "$1"
    log "ERROR" "任务失败,退出码: ${2:-1}"
    exit "${2:-1}"
}

# 打印分隔线
print_separator() {
    echo "=======================================================================" | tee -a "$LOG_FILE"
}

# =================== 预检查阶段 ===================

log "INFO" "Joint Calling 任务启动"
print_separator

# 1. 创建必要目录
log "INFO" "创建工作目录..."
mkdir -p "$OUT_DIR" || error_exit "无法创建输出目录: $OUT_DIR"
mkdir -p "$TMP_DIR" || error_exit "无法创建临时目录: $TMP_DIR"

# 2. 检查依赖工具
log "INFO" "检查依赖工具..."
for tool in bcftools find; do
    if ! command -v $tool &> /dev/null; then
        error_exit "$tool 未安装或不在 PATH 中"
    fi
done

# 检查 samtools (用于索引参考基因组)
if ! command -v samtools &> /dev/null; then
    log "WARN" "未找到 samtools,如果缺少 .fai 索引文件将无法自动创建"
fi

# 获取 bcftools 版本
BCFTOOLS_VERSION=$(bcftools --version | head -n1)
log "INFO" "使用工具版本: $BCFTOOLS_VERSION"

# 检查 bcftools 版本是否支持所需功能
BCFTOOLS_MAJOR_VERSION=$(echo "$BCFTOOLS_VERSION" | grep -oP '\d+\.\d+' | head -n1 | cut -d. -f1)
if [ "$BCFTOOLS_MAJOR_VERSION" -lt 1 ]; then
    log "WARN" "bcftools 版本较旧,建议使用 1.10 或更高版本"
fi

# 3. 检查参考基因组
log "INFO" "检查参考基因组..."
if [ ! -f "$REF_GENOME" ] || [[ "$REF_GENOME" == "/path/to/your/reference_genome.fa" ]]; then
    print_separator
    error_exit "未找到参考基因组文件!\n请编辑脚本中的 REF_GENOME 变量。\nJoint Calling 必须依赖参考基因组进行 Indel 矫正和分型。"
fi

# 验证参考基因组文件格式
log "INFO" "验证参考基因组文件格式..."
if ! head -n 1 "$REF_GENOME" | grep -q "^>"; then
    error_exit "参考基因组文件格式错误,不是有效的 FASTA 格式: $REF_GENOME"
fi

# 检查参考基因组索引
if [ ! -f "${REF_GENOME}.fai" ]; then
    log "WARN" "参考基因组缺少索引文件 (.fai),正在创建..."
    
    # 检查是否有 samtools
    if command -v samtools &> /dev/null; then
        samtools faidx "$REF_GENOME" || error_exit "无法为参考基因组创建索引"
        log "INFO" "参考基因组索引创建完成"
    else
        error_exit "未找到 samtools 工具,无法创建参考基因组索引\n请先运行: samtools faidx $REF_GENOME"
    fi
else
    log "INFO" "参考基因组索引文件存在"
fi

log "INFO" "参考基因组: $REF_GENOME"
log "INFO" "参考基因组大小: $(du -h "$REF_GENOME" | awk '{print $1}')"

# 4. 检查输入目录
log "INFO" "检查输入目录..."
if [ ! -d "$IN_DIR" ]; then
    error_exit "输入目录不存在: $IN_DIR"
fi

# 5. 生成文件列表
log "INFO" "搜索 g.vcf.gz 文件..."
find "$IN_DIR" -type f -name "*.g.vcf.gz" | sort > "$LIST_FILE"

FILE_COUNT=$(wc -l < "$LIST_FILE")
log "INFO" "找到 $FILE_COUNT 个 g.vcf.gz 文件"

if [ "$FILE_COUNT" -eq 0 ]; then
    error_exit "在输入目录未找到 .g.vcf.gz 文件: $IN_DIR"
fi

# 6. 检查索引文件
log "INFO" "检查输入文件索引..."
MISSING_INDEX=0
MISSING_INDEX_FILES=()

while IFS= read -r vcf_file; do
    if [ ! -f "${vcf_file}.tbi" ] && [ ! -f "${vcf_file}.csi" ]; then
        ((MISSING_INDEX++))
        MISSING_INDEX_FILES+=("$vcf_file")
    fi
done < "$LIST_FILE"

if [ $MISSING_INDEX -gt 0 ]; then
    log "WARN" "发现 $MISSING_INDEX 个文件缺少索引"
    log "INFO" "正在创建缺失的索引文件..."
    
    INDEX_COUNT=0
    for vcf_file in "${MISSING_INDEX_FILES[@]}"; do
        ((INDEX_COUNT++))
        log "DEBUG" "[$INDEX_COUNT/$MISSING_INDEX] 为 $(basename "$vcf_file") 创建索引..."
        bcftools index -t "$vcf_file" --threads 4 || log "ERROR" "索引创建失败: $vcf_file"
    done
    log "INFO" "索引创建完成"
else
    log "INFO" "所有输入文件索引完整"
fi

# 7. VCF 完整性检查 (可选)
if [ "$CHECK_VCF_INTEGRITY" = true ]; then
    log "INFO" "进行 VCF 完整性抽样检查 (检查前10个文件)..."
    CHECK_COUNT=0
    CORRUPT_FILES=0
    HAS_QS_HEADER=0
    
    while IFS= read -r vcf_file && [ $CHECK_COUNT -lt 10 ]; do
        ((CHECK_COUNT++))
        if ! bcftools view -h "$vcf_file" > /dev/null 2>&1; then
            log "ERROR" "文件损坏或格式错误: $vcf_file"
            ((CORRUPT_FILES++))
        fi
        
        # 检查 Header 中是否定义了 QS 字段
        if bcftools view -h "$vcf_file" | grep -q "##FORMAT=.*<ID=QS"; then
            ((HAS_QS_HEADER++))
        fi
    done < "$LIST_FILE"
    
    if [ $CORRUPT_FILES -gt 0 ]; then
        error_exit "发现 $CORRUPT_FILES 个损坏的 VCF 文件,请检查输入数据"
    fi
    
    # 判断是否需要处理 QS 问题
    if [ $HAS_QS_HEADER -gt 0 ]; then
        log "WARN" "检测到 $HAS_QS_HEADER 个文件在 Header 中定义了 QS 字段"
        log "WARN" "这可能导致 bcftools merge 失败 (QS annotation not present)"
        
        if [ "$AUTO_REMOVE_QS" = true ]; then
            log "INFO" "将预处理所有文件,移除 QS 及其他问题标签"
            STRIP_PROBLEMATIC_TAGS=true
        else
            log "ERROR" "请设置 AUTO_REMOVE_QS=true 以自动移除 QS 标签"
            error_exit "VCF 文件存在 QS 标签问题"
        fi
    fi
    
    log "INFO" "VCF 完整性检查完成"
fi

# 8. 检查文件句柄限制
log "INFO" "检查系统资源限制..."
LIMIT_OPEN_FILES=$(ulimit -n)
log "INFO" "当前最大打开文件数: $LIMIT_OPEN_FILES"

if [ "$LIMIT_OPEN_FILES" -lt "$MIN_OPEN_FILES" ]; then
    log "WARN" "文件句柄限制过低,尝试提升到 65535..."
    if ulimit -n 65535 2>/dev/null; then
        log "INFO" "文件句柄限制已提升到: $(ulimit -n)"
    else
        log "WARN" "无法自动提升限制,建议手动运行: ulimit -n 65535"
    fi
fi

# 9. 磁盘空间检查
log "INFO" "检查磁盘空间..."
AVAILABLE_SPACE_KB=$(df -k "$OUT_DIR" | awk 'NR==2 {print $4}')
AVAILABLE_SPACE_GB=$((AVAILABLE_SPACE_KB / 1024 / 1024))

# 估算输入总大小
TOTAL_INPUT_SIZE_KB=$(du -sk $(cat "$LIST_FILE") 2>/dev/null | awk '{s+=$1} END {print s}')
TOTAL_INPUT_SIZE_GB=$((TOTAL_INPUT_SIZE_KB / 1024 / 1024))

# Joint Calling 后大小通常会减小,预留1.5倍空间
REQUIRED_SPACE_GB=$((TOTAL_INPUT_SIZE_GB * 3 / 2))

log "INFO" "输入总大小: ${TOTAL_INPUT_SIZE_GB}GB"
log "INFO" "可用磁盘空间: ${AVAILABLE_SPACE_GB}GB"
log "INFO" "预估所需空间: ${REQUIRED_SPACE_GB}GB"

if [ "$AVAILABLE_SPACE_GB" -lt "$REQUIRED_SPACE_GB" ]; then
    log "WARN" "磁盘空间可能不足,但继续执行"
else
    log "INFO" "磁盘空间充足"
fi

# 10. 备份已存在的输出文件
if [ -f "$OUT_FILE" ]; then
    if [ "$BACKUP_EXISTING" = true ]; then
        BACKUP_FILE="${OUT_FILE}.backup_$(date +%Y%m%d_%H%M%S)"
        log "WARN" "输出文件已存在,备份到: $BACKUP_FILE"
        mv "$OUT_FILE" "$BACKUP_FILE"
        [ -f "${OUT_FILE}.tbi" ] && mv "${OUT_FILE}.tbi" "${BACKUP_FILE}.tbi"
    else
        log "WARN" "输出文件已存在,将被覆盖: $OUT_FILE"
        rm -f "$OUT_FILE" "${OUT_FILE}.tbi"
    fi
fi

# =================== 预处理阶段 (移除问题标签) ===================

PROCESSED_LIST="$LIST_FILE"

if [ "$STRIP_PROBLEMATIC_TAGS" = true ]; then
    print_separator
    log "INFO" "预处理阶段: 移除问题标签 (QS, RGQ, MIN_DP)"
    log "INFO" "这是解决 'QS annotation not present' 错误的关键步骤"
    
    # 创建预处理输出目录
    PROCESSED_DIR="${TMP_DIR}/processed_vcf"
    mkdir -p "$PROCESSED_DIR"
    
    # 生成处理后的文件列表
    PROCESSED_LIST="${TMP_DIR}/processed_vcf_list.txt"
    > "$PROCESSED_LIST"  # 清空文件
    
    log "INFO" "正在处理 $FILE_COUNT 个文件..."
    
    PROCESS_FAILED=0
    
    # 并行处理多个文件
    PARALLEL_JOBS=$((THREADS_MERGE / 4))  # 使用部分线程进行预处理
    [ $PARALLEL_JOBS -lt 4 ] && PARALLEL_JOBS=4
    [ $PARALLEL_JOBS -gt 20 ] && PARALLEL_JOBS=20  # 限制最大并行数
    
    log "INFO" "使用 $PARALLEL_JOBS 个并行任务进行预处理"
    
    # 使用 GNU parallel 或者循环处理
    if command -v parallel &> /dev/null; then
        log "INFO" "使用 GNU parallel 加速处理"
        
        # 创建处理函数
        export PROCESSED_DIR
        export -f log 2>/dev/null || true
        
        # 使用 parallel 处理,每个任务的输出文件路径只写入一次
        cat "$LIST_FILE" | parallel -j "$PARALLEL_JOBS" --bar --joblog "${TMP_DIR}/parallel.log" \
            'vcf_basename=$(basename {}); \
             output_file="'"$PROCESSED_DIR"'/${vcf_basename}"; \
             bcftools annotate -x FORMAT/QS,FORMAT/RGQ,FORMAT/MIN_DP -O z -o "$output_file" {} 2>/dev/null && \
             bcftools index -t "$output_file" 2>/dev/null && \
             echo "$output_file"' > "$PROCESSED_LIST" 2>> "$LOG_FILE"
        
        PROCESS_FAILED=$?
        
        if [ $PROCESS_FAILED -ne 0 ]; then
            log "ERROR" "GNU parallel 处理过程中出现错误"
            log "INFO" "查看详细日志: ${TMP_DIR}/parallel.log"
        fi
        
    else
        log "INFO" "使用串行处理 (建议安装 GNU parallel 以加速)"
        
        PROCESS_COUNT=0
        
        while IFS= read -r vcf_file; do
            ((PROCESS_COUNT++))
            
            vcf_basename=$(basename "$vcf_file")
            output_file="${PROCESSED_DIR}/${vcf_basename}"
            
            # 显示进度
            if [ $((PROCESS_COUNT % 100)) -eq 0 ] || [ $PROCESS_COUNT -eq 1 ]; then
                log "INFO" "处理进度: $PROCESS_COUNT/$FILE_COUNT ($(awk "BEGIN {printf \"%.1f\", $PROCESS_COUNT/$FILE_COUNT*100}")%)"
            fi
            
            # 移除问题标签并重新压缩
            if bcftools annotate \
                -x FORMAT/QS,FORMAT/RGQ,FORMAT/MIN_DP \
                -O z \
                -o "$output_file" \
                "$vcf_file" 2>> "$LOG_FILE"; then
                
                # 创建索引
                if bcftools index -t "$output_file" 2>> "$LOG_FILE"; then
                    # 添加到处理后的列表
                    echo "$output_file" >> "$PROCESSED_LIST"
                else
                    log "ERROR" "索引创建失败: $output_file"
                    ((PROCESS_FAILED++))
                fi
            else
                log "ERROR" "处理失败: $vcf_file"
                ((PROCESS_FAILED++))
            fi
            
        done < "$LIST_FILE"
    fi
    
    # 验证处理后的文件数量
    PROCESSED_COUNT=$(wc -l < "$PROCESSED_LIST")
    
    log "INFO" "预处理统计: 输入 $FILE_COUNT 个, 输出 $PROCESSED_COUNT 个"
    
    if [ $PROCESS_FAILED -gt 0 ]; then
        log "ERROR" "有 $PROCESS_FAILED 个文件预处理失败"
    fi
    
    if [ "$PROCESSED_COUNT" -ne "$FILE_COUNT" ]; then
        log "ERROR" "预处理后文件数量不匹配: 期望 $FILE_COUNT, 实际 $PROCESSED_COUNT"
        
        # 显示一些调试信息
        log "DEBUG" "原始文件列表前5行:"
        head -n 5 "$LIST_FILE" | tee -a "$LOG_FILE"
        
        log "DEBUG" "处理后文件列表前5行:"
        head -n 5 "$PROCESSED_LIST" | tee -a "$LOG_FILE"
        
        log "DEBUG" "处理后文件列表后5行:"
        tail -n 5 "$PROCESSED_LIST" | tee -a "$LOG_FILE"
        
        # 检查是否有重复
        DUPLICATE_COUNT=$(sort "$PROCESSED_LIST" | uniq -d | wc -l)
        if [ $DUPLICATE_COUNT -gt 0 ]; then
            log "ERROR" "发现 $DUPLICATE_COUNT 个重复条目"
            log "DEBUG" "重复文件示例:"
            sort "$PROCESSED_LIST" | uniq -d | head -n 5 | tee -a "$LOG_FILE"
        fi
        
        error_exit "预处理失败,文件数量不匹配"
    fi
    
    log "INFO" "预处理完成: 成功处理 $PROCESSED_COUNT 个文件"
    log "INFO" "处理后的文件位于: $PROCESSED_DIR"
    log "INFO" "处理后的文件列表: $PROCESSED_LIST"
    print_separator
fi

# =================== Joint Calling 执行阶段 ===================

print_separator

# 根据文件类型选择处理方式
if [ "$USE_NAIVE_MERGE" = true ]; then
    log "INFO" "开始合并 VCF 文件 (简单合并模式)"
    log "WARN" "注意: 将直接合并变异位点,不重新计算基因型"
else
    log "INFO" "开始 Joint Calling (Merge + Call)"
fi

log "INFO" "输入文件数: $FILE_COUNT"
log "INFO" "使用文件列表: $PROCESSED_LIST"
log "INFO" "倍性: ${PLOIDY}倍体"
log "INFO" "Merge 线程数: $THREADS_MERGE"

if [ "$USE_NAIVE_MERGE" = false ]; then
    log "INFO" "Call 线程数: $THREADS_CALL"
fi

log "INFO" "输出文件: $OUT_FILE"
print_separator

START_TIME=$(date +%s)

# 使用临时文件避免中断时损坏输出
TEMP_OUT="${TMP_DIR}/$(basename "$OUT_FILE").tmp"

# 根据模式执行不同的命令
if [ "$USE_NAIVE_MERGE" = true ]; then
    # ========== 简单合并模式 ==========
    log "INFO" "执行简单合并..."
    
    bcftools merge \
        -l "$PROCESSED_LIST" \
        --threads "$THREADS_MERGE" \
        -O z \
        -o "$TEMP_OUT" 2>&1 | tee -a "$LOG_FILE"
    
    EXIT_CODE=${PIPESTATUS[0]}
    if [ $EXIT_CODE -ne 0 ]; then
        rm -f "$TEMP_OUT"
        error_exit "bcftools merge 失败" $EXIT_CODE
    fi
    
else
    # ========== Joint Calling 模式 ==========
    log "INFO" "执行 Joint Calling 管道..."
    
    # 构建 bcftools call 参数
    CALL_PARAMS="-m"  # multiallelic caller
    CALL_PARAMS="$CALL_PARAMS --ploidy $PLOIDY"
    
    if [ "$CALL_OUTPUT_ALL_SITES" = false ]; then
        CALL_PARAMS="$CALL_PARAMS -v"  # 只输出变异位点
        log "INFO" "模式: 仅输出变异位点 (variants only)"
    else
        log "INFO" "模式: 输出所有位点 (包括 0/0)"
    fi
    
    log "INFO" "倍性设置: ${PLOIDY}倍体"
    
    # 标准管道: merge -> call
    {
        bcftools merge \
            -l "$PROCESSED_LIST" \
            --threads "$THREADS_MERGE" \
            -O u \
        | bcftools call \
            $CALL_PARAMS \
            --threads "$THREADS_CALL" \
            -O z \
            -o "$TEMP_OUT"
    } 2>&1 | tee -a "$LOG_FILE"
    
    EXIT_CODE=${PIPESTATUS[0]}
    if [ $EXIT_CODE -ne 0 ]; then
        rm -f "$TEMP_OUT"
        error_exit "bcftools merge 失败" $EXIT_CODE
    fi
    
    EXIT_CODE=${PIPESTATUS[1]}
    if [ $EXIT_CODE -ne 0 ]; then
        rm -f "$TEMP_OUT"
        error_exit "bcftools call 失败" $EXIT_CODE
    fi
fi

# 移动临时文件到最终位置
log "INFO" "移动临时文件到最终位置..."
mv "$TEMP_OUT" "$OUT_FILE" || error_exit "无法移动输出文件"

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))
ELAPSED_MIN=$((ELAPSED / 60))
ELAPSED_SEC=$((ELAPSED % 60))

if [ "$USE_NAIVE_MERGE" = true ]; then
    log "INFO" "VCF 合并完成,耗时: ${ELAPSED_MIN}分${ELAPSED_SEC}秒"
else
    log "INFO" "Joint Calling 完成,耗时: ${ELAPSED_MIN}分${ELAPSED_SEC}秒"
fi

# =================== 索引和验证阶段 ===================

log "INFO" "创建输出文件索引..."
bcftools index -t "$OUT_FILE" --threads "$THREADS_INDEX" 2>&1 | tee -a "$LOG_FILE" \
    || error_exit "索引创建失败"

log "INFO" "索引创建完成"

# 验证输出文件
log "INFO" "验证输出文件..."
if bcftools view -h "$OUT_FILE" > /dev/null 2>&1; then
    log "INFO" "输出文件格式验证通过"
else
    error_exit "输出文件验证失败,可能已损坏"
fi

# =================== 统计信息 ===================

log "INFO" "生成统计信息..."

# 样本数
SAMPLE_COUNT=$(bcftools query -l "$OUT_FILE" 2>/dev/null | wc -l)
log "INFO" "样本数: $SAMPLE_COUNT"

# 变异位点数
log "INFO" "统计变异位点数 (可能需要一些时间)..."
VARIANT_COUNT=$(bcftools view -H "$OUT_FILE" 2>/dev/null | wc -l)
log "INFO" "变异位点数: $VARIANT_COUNT"

# 文件大小
FILE_SIZE=$(du -h "$OUT_FILE" | awk '{print $1}')
log "INFO" "输出文件大小: $FILE_SIZE"

# SNP vs Indel 统计
log "INFO" "统计变异类型分布..."
SNP_COUNT=$(bcftools view -v snps "$OUT_FILE" 2>/dev/null | bcftools view -H | wc -l)
INDEL_COUNT=$(bcftools view -v indels "$OUT_FILE" 2>/dev/null | bcftools view -H | wc -l)
log "INFO" "SNP 数量: $SNP_COUNT"
log "INFO" "Indel 数量: $INDEL_COUNT"

# 生成摘要文件
SUMMARY_FILE="${OUT_DIR}/joint_calling_summary.txt"

if [ "$USE_NAIVE_MERGE" = true ]; then
    PROCESS_TYPE="简单合并 (Naive Merge)"
else
    PROCESS_TYPE="Joint Calling"
fi

cat > "$SUMMARY_FILE" << EOF
VCF 处理任务摘要
=====================
完成时间: $(date)
处理类型: $PROCESS_TYPE
输入文件数: $FILE_COUNT
输出文件: $OUT_FILE
参考基因组: $REF_GENOME

统计信息:
---------
样本数: $SAMPLE_COUNT
总变异数: $VARIANT_COUNT
  - SNPs: $SNP_COUNT
  - Indels: $INDEL_COUNT
文件大小: $FILE_SIZE

资源使用:
---------
Merge 线程: $THREADS_MERGE
Call 线程: $THREADS_CALL
倍性设置: ${PLOIDY}倍体
总耗时: ${ELAPSED_MIN}分${ELAPSED_SEC}秒

日志文件: $LOG_FILE
EOF

log "INFO" "摘要已保存到: $SUMMARY_FILE"

# =================== 清理和完成 ===================

log "INFO" "清理临时文件..."
rm -rf "$TMP_DIR"

print_separator
log "INFO" "所有任务完成!"
log "INFO" "结果文件: $OUT_FILE"
log "INFO" "索引文件: ${OUT_FILE}.tbi"
log "INFO" "日志文件: $LOG_FILE"
log "INFO" "摘要文件: $SUMMARY_FILE"
print_separator

exit 0