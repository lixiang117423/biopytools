#!/bin/bash

# BAM文件覆盖度统计脚本
# 作者: Claude
# 日期: 2025-09-29

# 设置输入和输出目录
BAM_DIR="/share/org/YZWL/yzwl_lixg/project/16.荠菜/02.each/bam文件"
OUTPUT_DIR="/share/org/YZWL/yzwl_lixg/project/16.荠菜/02.each/coverage"

# 创建输出目录（如果不存在）
mkdir -p ${OUTPUT_DIR}

# 创建汇总文件
SUMMARY_FILE="${OUTPUT_DIR}/coverage_summary.txt"
echo -e "Sample\tAverage_Coverage\tCovered_Bases\tTotal_Bases\tCoverage_Rate" > ${SUMMARY_FILE}

# 检查BAM目录是否存在
if [ ! -d "${BAM_DIR}" ]; then
    echo "错误: BAM目录不存在: ${BAM_DIR}"
    exit 1
fi

# 统计BAM文件数量
BAM_COUNT=$(find ${BAM_DIR} -name "*.bam" | wc -l)
echo "找到 ${BAM_COUNT} 个BAM文件"

# 计数器
count=0

# 遍历所有BAM文件
for bam_file in ${BAM_DIR}/*.bam; do
    # 检查文件是否存在
    if [ ! -f "${bam_file}" ]; then
        echo "警告: 没有找到BAM文件"
        continue
    fi
    
    count=$((count + 1))
    
    # 获取文件名（不含路径和扩展名）
    sample_name=$(basename ${bam_file} .bam)
    
    echo "[$count/$BAM_COUNT] 正在处理: ${sample_name}"
    
    # 检查BAM索引文件是否存在，如果不存在则创建
    if [ ! -f "${bam_file}.bai" ]; then
        echo "  创建索引文件..."
        samtools index ${bam_file}
    fi
    
    # 使用samtools coverage计算覆盖度
    coverage_file="${OUTPUT_DIR}/${sample_name}_coverage.txt"
    echo "  计算覆盖度..."
    samtools coverage ${bam_file} > ${coverage_file}
    
    # 使用samtools depth计算详细的平均覆盖度（可选，更精确但更慢）
    # depth_file="${OUTPUT_DIR}/${sample_name}_depth.txt.gz"
    # echo "  计算深度..."
    # samtools depth -a ${bam_file} | gzip > ${depth_file}
    
    # 提取平均覆盖度信息并添加到汇总文件
    if [ -f "${coverage_file}" ]; then
        # 跳过标题行，计算所有染色体的加权平均覆盖度
        avg_cov=$(awk 'NR>1 {sum+=$7*$3; len+=$3} END {if(len>0) print sum/len; else print 0}' ${coverage_file})
        covered_bases=$(awk 'NR>1 {sum+=$5} END {print sum}' ${coverage_file})
        total_bases=$(awk 'NR>1 {sum+=$3} END {print sum}' ${coverage_file})
        coverage_rate=$(awk 'NR>1 {sum+=$6*$3; len+=$3} END {if(len>0) print sum/len; else print 0}' ${coverage_file})
        
        echo -e "${sample_name}\t${avg_cov}\t${covered_bases}\t${total_bases}\t${coverage_rate}" >> ${SUMMARY_FILE}
        
        echo "  平均覆盖度: ${avg_cov}x"
    fi
    
    echo "  完成!"
    echo ""
done

echo "所有文件处理完成!"
echo "详细结果保存在: ${OUTPUT_DIR}"
echo "汇总结果保存在: ${SUMMARY_FILE}"

# 显示汇总统计
echo ""
echo "========== 覆盖度汇总 =========="
column -t ${SUMMARY_FILE}
