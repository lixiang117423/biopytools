#!/bin/bash

# ==============================================================================
# 🧬 Dsuite Dtrios 分析脚本 (增强版)
# ==============================================================================
# 功能: 基因渗入分析 (Introgression Analysis)
# 作者: Enhanced Version
# 日期: $(date +%Y-%m-%d)
# ==============================================================================

set -euo pipefail  # 严格模式: 遇到错误立即退出

# ==============================================================================
# 📋 配置部分
# ==============================================================================

# Dsuite 可执行程序路径
DSUITE_BIN="/share/org/YZWL/yzwl_lixg/software/Dsuite/Build/Dsuite"

# 输入文件路径
INPUT_VCF="/share/org/YZWL/yzwl_lixg/project/19.大豆疫霉菌/12.D检验/variation.filtered.snp.vcf.gz"
INPUT_SETS="/share/org/YZWL/yzwl_lixg/project/19.大豆疫霉菌/12.D检验/55.D检验的分组信息.txt"

# 输出设置
OUT_DIR="/share/org/YZWL/yzwl_lixg/project/19.大豆疫霉菌/12.D检验/output"
OUT_PREFIX="${OUT_DIR}/psoja"

# 日志文件
LOG_FILE="${OUT_DIR}/dsuite_analysis_$(date +%Y%m%d_%H%M%S).log"

# bcftools 命令
BCFTOOLS="bcftools"

# 分析参数
MIN_ALLELES=2          # 最小等位基因数
MAX_ALLELES=2          # 最大等位基因数 (双等位)
VARIANT_TYPE="snps"    # 变异类型: snps, indels, both
THREADS=64              # 线程数 (如果支持)

# ==============================================================================
# 🛠️ 工具函数
# ==============================================================================

# 日志函数
log() {
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    echo "[$timestamp] $1" | tee -a "$LOG_FILE"
}

# 错误处理函数
error_exit() {
    log "❌ 错误: $1"
    exit 1
}

# 成功信息函数
success() {
    log "✅ $1"
}

# 警告信息函数
warning() {
    log "⚠️ 警告: $1"
}

# 检查文件是否存在
check_file() {
    if [ ! -f "$1" ]; then
        error_exit "文件不存在: $1"
    fi
}

# 检查命令是否存在
check_command() {
    if ! command -v "$1" &> /dev/null; then
        error_exit "未找到命令: $1。请先安装或加载它"
    fi
}

# 打印分隔线
print_separator() {
    log "=================================================================="
}

# ==============================================================================
# 🚀 主程序开始
# ==============================================================================

main() {
    # 打印欢迎信息
    print_separator
    log "🧬 Dsuite Dtrios 分析脚本启动"
    log "分析项目: 大豆疫霉菌基因渗入检测"
    print_separator
    
    # 创建输出目录
    if [ ! -d "$OUT_DIR" ]; then
        log "📂 创建输出目录: $OUT_DIR"
        mkdir -p "$OUT_DIR" || error_exit "无法创建输出目录"
    fi
    
    # 初始化日志
    log "📝 日志文件: $LOG_FILE"
    
    # ==============================================================================
    # 1️⃣ 环境检查
    # ==============================================================================
    log ""
    log "🔍 第一步: 环境检查"
    print_separator
    
    # 检查必要的命令
    check_command "$BCFTOOLS"
    success "bcftools 可用: $(command -v $BCFTOOLS)"
    
    # 检查 Dsuite
    if [ ! -f "$DSUITE_BIN" ]; then
        error_exit "Dsuite 程序不存在: $DSUITE_BIN"
    fi
    if [ ! -x "$DSUITE_BIN" ]; then
        error_exit "Dsuite 程序不可执行: $DSUITE_BIN"
    fi
    success "Dsuite 可用: $DSUITE_BIN"
    
    # 检查输入文件
    check_file "$INPUT_VCF"
    success "VCF 文件存在: $INPUT_VCF"
    
    check_file "$INPUT_SETS"
    success "分组文件存在: $INPUT_SETS"
    
    # 检查 VCF 文件格式
    if [[ "$INPUT_VCF" == *.gz ]]; then
        if ! $BCFTOOLS view -h "$INPUT_VCF" &>/dev/null; then
            error_exit "VCF 文件格式错误或已损坏"
        fi
        success "VCF 文件格式验证通过"
    fi
    
    # 显示分组信息
    log ""
    log "📊 分组信息预览 (前10行):"
    head -n 10 "$INPUT_SETS" | while IFS= read -r line; do
        log "  $line"
    done
    
    # ==============================================================================
    # 2️⃣ VCF 统计信息
    # ==============================================================================
    log ""
    log "📈 第二步: VCF 文件统计"
    print_separator
    
    # 统计样本数
    SAMPLE_COUNT=$($BCFTOOLS query -l "$INPUT_VCF" | wc -l)
    log "样本数量: $SAMPLE_COUNT"
    
    # 统计总变异数
    TOTAL_VARIANTS=$($BCFTOOLS view -H "$INPUT_VCF" | wc -l)
    log "总变异数: $TOTAL_VARIANTS"
    
    # 统计染色体/scaffold 数量
    CHROM_COUNT=$($BCFTOOLS view -h "$INPUT_VCF" | grep "^##contig" | wc -l)
    log "染色体/Scaffold 数: $CHROM_COUNT"
    
    # ==============================================================================
    # 3️⃣ 计算过滤后的变异数
    # ==============================================================================
    log ""
    log "🔬 第三步: 应用过滤条件并统计"
    print_separator
    
    log "⚙️ 过滤参数:"
    log "  - 最小等位基因数: $MIN_ALLELES"
    log "  - 最大等位基因数: $MAX_ALLELES"
    log "  - 变异类型: $VARIANT_TYPE"
    
    log "⏳ 正在统计符合条件的变异数..."
    
    # 计算过滤后的行数
    NUMLINES=$($BCFTOOLS view \
        -m${MIN_ALLELES} \
        -M${MAX_ALLELES} \
        -v ${VARIANT_TYPE} \
        "$INPUT_VCF" | grep -v "^#" | wc -l)
    
    if [ "$NUMLINES" -eq 0 ]; then
        error_exit "过滤后没有变异位点！请检查过滤参数"
    fi
    
    success "符合条件的双等位 SNP 数: $NUMLINES"
    
    # 计算过滤比例
    FILTER_RATIO=$(awk "BEGIN {printf \"%.2f\", ($NUMLINES/$TOTAL_VARIANTS)*100}")
    log "过滤保留比例: ${FILTER_RATIO}%"
    
    # ==============================================================================
    # 4️⃣ 运行 Dsuite Dtrios
    # ==============================================================================
    log ""
    log "🏃 第四步: 运行 Dsuite Dtrios"
    print_separator
    
    log "💻 命令预览:"
    log "  $BCFTOOLS view -m${MIN_ALLELES} -M${MAX_ALLELES} -v ${VARIANT_TYPE} $INPUT_VCF | \\"
    log "  $DSUITE_BIN Dtrios -l $NUMLINES -o $OUT_PREFIX stdin $INPUT_SETS"
    
    log ""
    log "⏱️ 开始时间: $(date '+%Y-%m-%d %H:%M:%S')"
    START_TIME=$(date +%s)
    
    # 执行分析
    set +e  # 临时关闭严格模式以捕获退出码
    $BCFTOOLS view \
        -m${MIN_ALLELES} \
        -M${MAX_ALLELES} \
        -v ${VARIANT_TYPE} \
        "$INPUT_VCF" | \
    $DSUITE_BIN Dtrios \
        -l $NUMLINES \
        -o "$OUT_PREFIX" \
        stdin \
        "$INPUT_SETS" 2>&1 | tee -a "$LOG_FILE"
    
    EXIT_CODE=${PIPESTATUS[1]}  # 获取 Dsuite 的退出码
    set -e  # 恢复严格模式
    
    END_TIME=$(date +%s)
    ELAPSED=$((END_TIME - START_TIME))
    
    log "⏱️ 结束时间: $(date '+%Y-%m-%d %H:%M:%S')"
    log "⏱️ 总耗时: ${ELAPSED} 秒 ($(($ELAPSED/60)) 分钟)"
    
    # ==============================================================================
    # 5️⃣ 结果检查
    # ==============================================================================
    log ""
    log "📊 第五步: 结果验证"
    print_separator
    
    if [ $EXIT_CODE -eq 0 ]; then
        success "Dsuite 运行成功！"
        
        # 列出生成的文件
        log ""
        log "📁 生成的输出文件:"
        if ls "${OUT_PREFIX}"* 1> /dev/null 2>&1; then
            for file in "${OUT_PREFIX}"*; do
                if [ -f "$file" ]; then
                    FILE_SIZE=$(du -h "$file" | cut -f1)
                    log "  ✓ $(basename $file) (大小: $FILE_SIZE)"
                fi
            done
        else
            warning "未找到输出文件"
        fi
        
        # 检查主要结果文件
        MAIN_RESULT="${OUT_PREFIX}_BBAA.txt"
        if [ -f "$MAIN_RESULT" ]; then
            log ""
            log "📄 主要结果文件预览 (${MAIN_RESULT}):"
            log "前10行:"
            head -n 10 "$MAIN_RESULT" | while IFS= read -r line; do
                log "  $line"
            done
            
            # 统计三元组数量
            TRIO_COUNT=$(($(wc -l < "$MAIN_RESULT") - 1))
            log ""
            log "🔢 分析的三元组数量: $TRIO_COUNT"
        fi
        
    else
        error_exit "Dsuite 运行失败 (退出码: $EXIT_CODE)"
    fi
    
    # ==============================================================================
    # 6️⃣ 生成分析总结
    # ==============================================================================
    log ""
    log "📋 分析总结"
    print_separator
    log "输入文件: $INPUT_VCF"
    log "分组文件: $INPUT_SETS"
    log "输出目录: $OUT_DIR"
    log "输出前缀: $OUT_PREFIX"
    log "样本数量: $SAMPLE_COUNT"
    log "总变异数: $TOTAL_VARIANTS"
    log "过滤后变异数: $NUMLINES (${FILTER_RATIO}%)"
    log "运行时间: ${ELAPSED} 秒"
    log "日志文件: $LOG_FILE"
    
    # ==============================================================================
    # 7️⃣ 下一步提示
    # ==============================================================================
    log ""
    log "💡 下一步建议"
    print_separator
    log "1. 查看 D 统计结果: cat ${OUT_PREFIX}_BBAA.txt"
    log "2. 查看树文件: cat ${OUT_PREFIX}_tree.txt"
    log "3. 使用 Dsuite Dinvestigate 进行进一步分析"
    log "4. 使用 R 脚本可视化结果 (Dsuite 提供的 utils/)"
    
    print_separator
    log "🎉 分析流程完成！"
    print_separator
}

# ==============================================================================
# 执行主程序
# ==============================================================================
main "$@"