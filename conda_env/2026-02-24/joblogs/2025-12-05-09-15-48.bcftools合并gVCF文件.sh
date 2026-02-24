#!/bin/bash

set -euo pipefail  # 遇到错误立即退出,未定义变量报错,管道命令失败时退出

# ================= 配置区域 =================
# 输入包含 g.vcf.gz 的文件夹路径
IN_DIR="/share/org/YZWL/yzwl_lixg/project/21.野生大豆群体/03.glnexus_joint/vcf"

# 输出文件夹路径
OUT_DIR="/share/org/YZWL/yzwl_lixg/project/21.野生大豆群体/03.第一批720个样品/03.bcftools_joint"

# 输出的合并文件名
OUT_FILE="${OUT_DIR}/all_samples_bcftools_joint.vcf.gz"

# 线程数
THREADS=88

# 日志文件
LOG_FILE="${OUT_DIR}/merge_$(date +%Y%m%d_%H%M%S).log"

# 临时目录
TMP_DIR="${OUT_DIR}/tmp"

# ===========================================

# 日志函数
log() {
    local level=$1
    shift
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [$level] $*" | tee -a "$LOG_FILE"
}

# 错误退出函数
error_exit() {
    log "ERROR" "$1"
    exit 1
}

# 开始记录
log "INFO" "========== VCF 合并任务开始 =========="
log "INFO" "输入目录: $IN_DIR"
log "INFO" "输出目录: $OUT_DIR"
log "INFO" "线程数: $THREADS"

# 1. 检查依赖工具
log "INFO" "检查必需工具..."
for tool in bcftools find; do
    if ! command -v $tool &> /dev/null; then
        error_exit "$tool 未安装或不在 PATH 中"
    fi
done
log "INFO" "工具检查完成"

# 2. 创建必要目录
log "INFO" "创建输出目录..."
mkdir -p "$OUT_DIR" || error_exit "无法创建输出目录"
mkdir -p "$TMP_DIR" || error_exit "无法创建临时目录"

# 3. 检查输入目录
if [ ! -d "$IN_DIR" ]; then
    error_exit "输入目录不存在: $IN_DIR"
fi

# 4. 生成文件列表
log "INFO" "搜索 g.vcf.gz 文件..."
VCF_LIST="${OUT_DIR}/vcf_file_list.txt"

find "$IN_DIR" -type f -name "*.g.vcf.gz" | sort > "$VCF_LIST"

FILE_COUNT=$(wc -l < "$VCF_LIST")
log "INFO" "找到 $FILE_COUNT 个 VCF 文件"

if [ "$FILE_COUNT" -eq 0 ]; then
    error_exit "在 ${IN_DIR} 未找到 .g.vcf.gz 文件"
fi

# 5. 检查索引文件
log "INFO" "检查索引文件..."
MISSING_INDEX=0
while IFS= read -r vcf_file; do
    if [ ! -f "${vcf_file}.tbi" ]; then
        log "WARN" "缺少索引: ${vcf_file}.tbi"
        ((MISSING_INDEX++))
    fi
done < "$VCF_LIST"

if [ $MISSING_INDEX -gt 0 ]; then
    log "WARN" "发现 $MISSING_INDEX 个文件缺少索引,将自动创建..."
    while IFS= read -r vcf_file; do
        if [ ! -f "${vcf_file}.tbi" ]; then
            log "INFO" "为 $(basename "$vcf_file") 创建索引..."
            bcftools index -t "$vcf_file" --threads 4 || log "ERROR" "索引创建失败: $vcf_file"
        fi
    done < "$VCF_LIST"
    log "INFO" "索引创建完成"
else
    log "INFO" "所有文件索引完整"
fi

# 6. 检查磁盘空间
log "INFO" "检查磁盘空间..."
AVAILABLE_SPACE=$(df -BG "$OUT_DIR" | awk 'NR==2 {print $4}' | sed 's/G//')
TOTAL_INPUT_SIZE=$(du -scBG $(cat "$VCF_LIST") | tail -n1 | awk '{print $1}' | sed 's/G//')
REQUIRED_SPACE=$((TOTAL_INPUT_SIZE * 2))  # 预估需要2倍输入大小

log "INFO" "可用空间: ${AVAILABLE_SPACE}G, 输入总大小: ${TOTAL_INPUT_SIZE}G, 预估需要: ${REQUIRED_SPACE}G"

if [ "$AVAILABLE_SPACE" -lt "$REQUIRED_SPACE" ]; then
    log "WARN" "磁盘空间可能不足,但继续执行"
else
    log "INFO" "磁盘空间充足"
fi

# 7. 执行合并
log "INFO" "开始合并 VCF 文件..."
log "INFO" "使用 $THREADS 线程"

# 备份旧文件
if [ -f "$OUT_FILE" ]; then
    BACKUP="${OUT_FILE}.backup_$(date +%Y%m%d_%H%M%S)"
    log "WARN" "输出文件已存在,备份至: $BACKUP"
    mv "$OUT_FILE" "$BACKUP"
    [ -f "${OUT_FILE}.tbi" ] && mv "${OUT_FILE}.tbi" "${BACKUP}.tbi"
fi

# 执行合并,使用临时文件避免中断时损坏输出
TEMP_OUT="${TMP_DIR}/$(basename "$OUT_FILE").tmp"

bcftools merge \
    -l "$VCF_LIST" \
    -O z \
    -o "$TEMP_OUT" \
    --threads "$THREADS" 2>&1 | tee -a "$LOG_FILE"

MERGE_EXIT_CODE=${PIPESTATUS[0]}

if [ $MERGE_EXIT_CODE -ne 0 ]; then
    error_exit "bcftools merge 失败,退出码: $MERGE_EXIT_CODE"
fi

# 移动临时文件到最终位置
mv "$TEMP_OUT" "$OUT_FILE" || error_exit "无法移动输出文件"
log "INFO" "合并完成,输出: $OUT_FILE"

# 8. 建立索引
log "INFO" "为输出文件创建索引..."
bcftools index -t "$OUT_FILE" --threads "$THREADS" 2>&1 | tee -a "$LOG_FILE" \
    || error_exit "索引创建失败"

log "INFO" "索引创建完成"

# 9. 验证输出
log "INFO" "验证输出文件..."
if bcftools view -H "$OUT_FILE" | head -n 1 > /dev/null 2>&1; then
    log "INFO" "输出文件验证通过"
else
    error_exit "输出文件验证失败"
fi

# 10. 统计信息
log "INFO" "生成统计信息..."
SAMPLE_COUNT=$(bcftools query -l "$OUT_FILE" | wc -l)
VARIANT_COUNT=$(bcftools view -H "$OUT_FILE" | wc -l)
FILE_SIZE=$(du -h "$OUT_FILE" | awk '{print $1}')

log "INFO" "样本数: $SAMPLE_COUNT"
log "INFO" "变异位点数: $VARIANT_COUNT"
log "INFO" "文件大小: $FILE_SIZE"

# 11. 清理临时文件
log "INFO" "清理临时文件..."
rm -rf "$TMP_DIR"

# 12. 完成
log "INFO" "========== 任务全部完成 =========="
log "INFO" "输出文件: $OUT_FILE"
log "INFO" "日志文件: $LOG_FILE"
log "INFO" "总耗时: $SECONDS 秒"

exit 0
