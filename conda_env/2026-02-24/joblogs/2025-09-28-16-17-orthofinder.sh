# #!/bin/bash

# # 🌿 叶绿体基因组OrthoFinder分析脚本
# # 📅 日期: $(date +%Y-%m-%d)
# # 🧬 用途: 分析477个十字花科植物叶绿体基因组

# echo "🌿=========================================="
# echo "🧬 叶绿体基因组OrthoFinder分析开始"
# echo "🌿=========================================="

# # 参数设置
# INPUT_DIR="/share/org/YZWL/yzwl_lixg/project/15.十字花科泛基因组/04.叶绿体基因组/unique_genome/fasta/pep/by_sample"
# TIMESTAMP=$(date +%Y%m%d_%H%M%S)
# OUTPUT_DIR="$(pwd)/chloroplast_orthofinder_analysis_${TIMESTAMP}"
# PROCESSED_DIR="${OUTPUT_DIR}/processed_fasta"
# RESULTS_DIR="${OUTPUT_DIR}/orthofinder_results"
# THREADS=88
# LOG_FILE="${OUTPUT_DIR}/analysis.log"

# echo "📂 输入目录: $INPUT_DIR"
# echo "📁 输出目录: $OUTPUT_DIR"

# # 创建目录
# echo "📁 创建目录结构..."
# mkdir -p "$OUTPUT_DIR"
# mkdir -p "$PROCESSED_DIR"
# mkdir -p "$RESULTS_DIR"
# echo "✅ 目录创建完成" | tee -a "$LOG_FILE"

# # 检查OrthoFinder
# echo "🔍 检查OrthoFinder..."
# if ! command -v orthofinder &> /dev/null; then
#     echo "❌ 错误: OrthoFinder未找到" | tee -a "$LOG_FILE"
#     exit 1
# fi
# echo "✅ OrthoFinder检查通过" | tee -a "$LOG_FILE"

# # 检查输入目录
# if [[ ! -d "$INPUT_DIR" ]]; then
#     echo "❌ 错误: 输入目录不存在: $INPUT_DIR" | tee -a "$LOG_FILE"
#     exit 1
# fi

# # 统计文件
# file_count=$(find "$INPUT_DIR" -name "*_pep.fasta" | wc -l)
# echo "🔍 发现 $file_count 个蛋白质文件" | tee -a "$LOG_FILE"

# if [[ $file_count -eq 0 ]]; then
#     echo "❌ 错误: 没有找到匹配的文件" | tee -a "$LOG_FILE"
#     exit 1
# fi

# # 处理FASTA文件
# echo "🧬 开始处理FASTA文件..." | tee -a "$LOG_FILE"

# processed_count=0
# for fasta_file in "$INPUT_DIR"/*_pep.fasta; do
#     if [[ ! -f "$fasta_file" ]]; then
#         continue
#     fi
    
#     # 获取样本名称
#     basename_file=$(basename "$fasta_file")
#     sample_name="${basename_file%_pep.fasta}"
    
#     # 输出文件
#     output_file="${PROCESSED_DIR}/${sample_name}_processed.faa"
    
#     # 处理文件
#     awk -v sample="$sample_name" '
#     /^>/ {
#         split($0, parts, " ")
#         gene_id = substr(parts[1], 2)
#         print ">" sample "_" gene_id
#         next
#     }
#     {print}
#     ' "$fasta_file" > "$output_file"
    
#     ((processed_count++))
    
#     if ((processed_count % 50 == 0)); then
#         echo "⏳ 已处理 $processed_count / $file_count 个文件..." | tee -a "$LOG_FILE"
#     fi
# done

# echo "✅ FASTA文件处理完成，共处理 $processed_count 个文件" | tee -a "$LOG_FILE"

# # 验证处理结果
# final_count=$(find "$PROCESSED_DIR" -name "*.faa" | wc -l)
# echo "📊 最终处理文件数: $final_count" | tee -a "$LOG_FILE"

# if [[ $final_count -eq 0 ]]; then
#     echo "❌ 错误: 没有成功处理任何文件" | tee -a "$LOG_FILE"
#     exit 1
# fi

# # 运行OrthoFinder
# echo "🚀 开始运行OrthoFinder分析..." | tee -a "$LOG_FILE"
# echo "💻 使用线程数: $THREADS" | tee -a "$LOG_FILE"

# orthofinder \
#     -f "$PROCESSED_DIR" \
#     -t "$THREADS" \
#     -a "$THREADS" \
#     -o "$RESULTS_DIR" \
#     -S diamond 2>&1 | tee -a "$LOG_FILE"

# if [[ $? -eq 0 ]]; then
#     echo "🎉 OrthoFinder分析成功完成" | tee -a "$LOG_FILE"
# else
#     echo "❌ OrthoFinder分析失败" | tee -a "$LOG_FILE"
#     exit 1
# fi

# # 分析结果
# echo "📊 分析结果统计..." | tee -a "$LOG_FILE"

# result_dir=$(find "$RESULTS_DIR" -name "Results_*" -type d | head -1)
# if [[ -n "$result_dir" ]]; then
#     echo "📂 结果目录: $result_dir" | tee -a "$LOG_FILE"
    
#     # 统计orthogroups
#     orthogroups_file="$result_dir/Orthogroups/Orthogroups.tsv"
#     if [[ -f "$orthogroups_file" ]]; then
#         orthogroup_count=$(tail -n +2 "$orthogroups_file" | wc -l)
#         echo "🔗 Orthogroup数量: $orthogroup_count" | tee -a "$LOG_FILE"
#     fi
    
#     # 统计单拷贝基因
#     single_copy_dir="$result_dir/Single_Copy_Orthologue_Sequences"
#     if [[ -d "$single_copy_dir" ]]; then
#         single_copy_count=$(ls "$single_copy_dir"/*.fa 2>/dev/null | wc -l)
#         echo "🎯 单拷贝基因数量: $single_copy_count" | tee -a "$LOG_FILE"
#     fi
    
#     # 基因树
#     gene_trees_dir="$result_dir/Gene_Trees"
#     if [[ -d "$gene_trees_dir" ]]; then
#         tree_count=$(ls "$gene_trees_dir"/*.txt 2>/dev/null | wc -l)
#         echo "🌳 基因树数量: $tree_count" | tee -a "$LOG_FILE"
#     fi
    
#     # 物种树
#     species_tree="$result_dir/Species_Tree/SpeciesTree_rooted.txt"
#     if [[ -f "$species_tree" ]]; then
#         echo "🌲 物种树: $species_tree" | tee -a "$LOG_FILE"
#     fi
# fi

# echo "🌿=========================================="
# echo "🎉 分析完成！"
# echo "📊 结果保存在: $RESULTS_DIR"
# echo "📝 日志文件: $LOG_FILE"
# echo "🌿=========================================="

biopytools orthofinder -i processed_fasta -o ./ -s diamond 