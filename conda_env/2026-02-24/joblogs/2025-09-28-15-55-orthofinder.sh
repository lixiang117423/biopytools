#!/bin/bash

# 🌿 叶绿体基因组OrthoFinder分析脚本
# 👨‍💻 作者: [Your Name]
# 📅 日期: $(date +%Y-%m-%d)
# 🔬 用途: 分析477个十字花科植物叶绿体基因组的直系同源关系

set -e  # 遇到错误时退出

# =============================================================================
# 参数设置
# =============================================================================

# 输入目录
INPUT_DIR="/share/org/YZWL/yzwl_lixg/project/15.十字花科泛基因组/04.叶绿体基因组/unique_genome/fasta/pep/by_sample"

# 输出目录
OUTPUT_DIR="$(pwd)/chloroplast_orthofinder_analysis"
PROCESSED_DIR="${OUTPUT_DIR}/processed_fasta"
RESULTS_DIR="${OUTPUT_DIR}/orthofinder_results"

# 线程数
THREADS=88

# 日志文件
LOG_FILE="${OUTPUT_DIR}/analysis.log"

# =============================================================================
# 函数定义
# =============================================================================

# 📝 日志记录函数
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$LOG_FILE"
}

# ❌ 错误处理函数
error_exit() {
    log "❌ 错误: $1"
    exit 1
}

# 🔍 检查依赖
check_dependencies() {
    log "🔍 检查依赖软件..."
    
    if ! command -v orthofinder &> /dev/null; then
        error_exit "OrthoFinder未找到，请确保已正确安装并添加到PATH"
    fi
    
    log "✅ 依赖检查完成"
}

# 📁 创建目录结构
setup_directories() {
    log "📁 创建目录结构..."
    
    mkdir -p "$OUTPUT_DIR"
    mkdir -p "$PROCESSED_DIR"
    mkdir -p "$RESULTS_DIR"
    
    log "✅ 目录创建完成"
}

# 🧬 处理FASTA文件，添加样本前缀避免基因名冲突
process_fasta_files() {
    log "🧬 开始处理FASTA文件，添加样本前缀..."
    
    local file_count=0
    local total_files=$(find "$INPUT_DIR" -name "*_pep.fasta" | wc -l)
    
    log "🔍 发现 $total_files 个蛋白质文件"
    
    for fasta_file in "$INPUT_DIR"/*_pep.fasta; do
        # 检查文件是否存在（避免通配符无匹配时的错误）
        if [[ ! -f "$fasta_file" ]]; then
            continue
        fi
        
        # 获取文件名并提取样本名称
        basename_file=$(basename "$fasta_file")
        # 去掉_pep.fasta后缀得到样本名称，如：Aethionema_arabicum_pep.fasta -> Aethionema_arabicum
        sample_name="${basename_file%_pep.fasta}"
        
        # 输出文件路径
        output_file="${PROCESSED_DIR}/${sample_name}_processed.faa"
        
        # 处理FASTA文件，给每个基因ID添加样本前缀
        awk -v sample="$sample_name" '
        /^>/ {
            # 提取基因名称（只取第一个字段，如：>accD Aethionema_arabicum | ... -> accD）
            split($0, parts, " ")
            gene_id = substr(parts[1], 2)  # 去掉>符号
            # 添加样本前缀
            print ">" sample "_" gene_id
            next
        }
        # 序列行直接输出
        {print}
        ' "$fasta_file" > "$output_file"
        
        ((file_count++))
        
        # 进度显示
        if ((file_count % 50 == 0)); then
            log "⏳ 已处理 $file_count / $total_files 个文件..."
        fi
    done
    
    log "✅ FASTA文件处理完成，共处理 $file_count 个文件"
    
    # 验证处理后的文件数量
    processed_count=$(find "$PROCESSED_DIR" -name "*.faa" | wc -l)
    log "📊 处理后的文件数量: $processed_count"
    
    if [[ $processed_count -eq 0 ]]; then
        error_exit "没有找到处理后的文件，请检查输入目录和文件格式"
    fi
}

# 🚀 运行OrthoFinder分析
run_orthofinder() {
    log "🚀 开始运行OrthoFinder分析..."
    log "💻 使用线程数: $THREADS"
    
    # OrthoFinder命令参数说明：
    # -f: 输入目录包含蛋白质序列文件
    # -t: 线程数
    # -a: 线程数用于BLAST搜索
    # -o: 输出目录
    # -S: 指定物种树推断方法 (diamond为快速模式)
    
    orthofinder \
        -f "$PROCESSED_DIR" \
        -t "$THREADS" \
        -a "$THREADS" \
        -o "$RESULTS_DIR" \
        -S diamond 2>&1 | tee -a "$LOG_FILE"
    
    if [[ ${PIPESTATUS[0]} -eq 0 ]]; then
        log "🎉 OrthoFinder分析成功完成"
    else
        error_exit "OrthoFinder分析失败"
    fi
}

# 📊 分析结果统计
analyze_results() {
    log "📊 分析结果统计..."
    
    # 查找结果目录
    result_subdir=$(find "$RESULTS_DIR" -name "Results_*" -type d | head -1)
    
    if [[ -z "$result_subdir" ]]; then
        log "⚠️  警告: 未找到Results目录"
        return 1
    fi
    
    log "📂 结果目录: $result_subdir"
    
    # 统计orthogroups
    orthogroups_file="$result_subdir/Orthogroups/Orthogroups.tsv"
    if [[ -f "$orthogroups_file" ]]; then
        orthogroup_count=$(tail -n +2 "$orthogroups_file" | wc -l)
        log "🔗 发现的Orthogroup数量: $orthogroup_count"
    fi
    
    # 统计单拷贝orthologs
    single_copy_file="$result_subdir/Single_Copy_Orthologue_Sequences"
    if [[ -d "$single_copy_file" ]]; then
        single_copy_count=$(ls "$single_copy_file"/*.fa 2>/dev/null | wc -l)
        log "🎯 单拷贝直系同源基因数量: $single_copy_count"
    fi
    
    # 统计基因树
    gene_trees_dir="$result_subdir/Gene_Trees"
    if [[ -d "$gene_trees_dir" ]]; then
        tree_count=$(ls "$gene_trees_dir"/*.txt 2>/dev/null | wc -l)
        log "🌳 构建的基因树数量: $tree_count"
    fi
    
    # 物种树
    species_tree="$result_subdir/Species_Tree/SpeciesTree_rooted.txt"
    if [[ -f "$species_tree" ]]; then
        log "🌲 物种树文件: $species_tree"
    fi
    
    log "📁 详细结果请查看目录: $result_subdir"
}

# 🧹 清理函数
cleanup() {
    log "🧹 清理临时文件..."
    # 如果需要，可以在这里添加清理代码
}

# =============================================================================
# 主程序
# =============================================================================

main() {
    log "🌿=========================================="
    log "🧬 叶绿体基因组OrthoFinder分析开始"
    log "🔬 十字花科植物477个样本直系同源分析"
    log "🌿=========================================="
    
    # 检查输入目录
    if [[ ! -d "$INPUT_DIR" ]]; then
        error_exit "输入目录不存在: $INPUT_DIR"
    fi
    
    log "📂 输入目录: $INPUT_DIR"
    log "📁 输出目录: $OUTPUT_DIR"
    
    # 执行分析步骤
    check_dependencies
    setup_directories
    process_fasta_files
    run_orthofinder
    analyze_results
    
    log "🌿=========================================="
    log "🎉 分析完成！"
    log "📊 结果保存在: $RESULTS_DIR"
    log "📝 日志文件: $LOG_FILE"
    log "🌿=========================================="
}

# 设置退出时的清理
trap cleanup EXIT

# 运行主程序
main "$@"
