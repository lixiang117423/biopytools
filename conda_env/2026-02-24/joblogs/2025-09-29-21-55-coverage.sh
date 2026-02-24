#!/bin/bash

# BAM文件覆盖度统计脚本（改进版）
# 作者: Claude
# 日期: 2025-09-29

# 设置输入和输出目录
BAM_DIR="/share/org/YZWL/yzwl_lixg/project/16.荠菜/02.each/bam"
OUTPUT_DIR="/share/org/YZWL/yzwl_lixg/project/16.荠菜/02.each/coverage"

# 创建输出目录（如果不存在）
mkdir -p ${OUTPUT_DIR}

# 创建日志文件
LOG_FILE="${OUTPUT_DIR}/coverage_log.txt"
echo "开始时间: $(date)" > ${LOG_FILE}

# 创建汇总文件
SUMMARY_FILE="${OUTPUT_DIR}/coverage_summary.txt"
echo -e "Sample\tTotal_Reads\tMapped_Reads\tAverage_Coverage\tCovered_Bases\tTotal_Bases\tCoverage_Rate" > ${SUMMARY_FILE}

# 检查BAM目录是否存在
if [ ! -d "${BAM_DIR}" ]; then
    echo "错误: BAM目录不存在: ${BAM_DIR}"
    exit 1
fi

# 检查samtools是否可用
if ! command -v samtools &> /dev/null; then
    echo "错误: samtools未安装或不在PATH中"
    exit 1
fi

echo "samtools版本: $(samtools --version | head -n1)"

# 统计BAM文件数量
BAM_COUNT=$(find ${BAM_DIR} -name "*.bam" | wc -l)
echo "找到 ${BAM_COUNT} 个BAM文件"
echo "找到 ${BAM_COUNT} 个BAM文件" >> ${LOG_FILE}

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
    
    echo ""
    echo "=========================================="
    echo "[$count/$BAM_COUNT] 正在处理: ${sample_name}"
    echo "=========================================="
    echo "" >> ${LOG_FILE}
    echo "[$count/$BAM_COUNT] 处理: ${sample_name}" >> ${LOG_FILE}
    
    # 检查BAM文件大小
    file_size=$(du -h ${bam_file} | cut -f1)
    echo "  文件大小: ${file_size}"
    echo "  文件大小: ${file_size}" >> ${LOG_FILE}
    
    # 使用samtools flagstat获取基本统计信息
    echo "  获取BAM文件统计信息..."
    flagstat_file="${OUTPUT_DIR}/${sample_name}_flagstat.txt"
    samtools flagstat ${bam_file} > ${flagstat_file}
    
    total_reads=$(grep "in total" ${flagstat_file} | awk '{print $1}')
    mapped_reads=$(grep "mapped (" ${flagstat_file} | head -n1 | awk '{print $1}')
    
    echo "  总reads数: ${total_reads}"
    echo "  已比对reads数: ${mapped_reads}"
    echo "  总reads数: ${total_reads}, 已比对: ${mapped_reads}" >> ${LOG_FILE}
    
    # 检查是否有比对的reads
    if [ "${mapped_reads}" = "0" ] || [ -z "${mapped_reads}" ]; then
        echo "  警告: 该BAM文件没有比对的reads！"
        echo "  警告: ${sample_name} 没有比对的reads" >> ${LOG_FILE}
        echo -e "${sample_name}\t${total_reads}\t0\t0\t0\t0\t0" >> ${SUMMARY_FILE}
        continue
    fi
    
    # 检查BAM索引文件是否存在，如果不存在则创建
    if [ ! -f "${bam_file}.bai" ]; then
        echo "  创建索引文件..."
        samtools index ${bam_file}
        if [ $? -ne 0 ]; then
            echo "  错误: 索引创建失败"
            echo "  错误: ${sample_name} 索引创建失败" >> ${LOG_FILE}
            continue
        fi
    fi
    
    # 使用samtools coverage计算覆盖度
    coverage_file="${OUTPUT_DIR}/${sample_name}_coverage.txt"
    echo "  计算覆盖度..."
    samtools coverage ${bam_file} > ${coverage_file}
    
    if [ $? -ne 0 ]; then
        echo "  错误: samtools coverage 执行失败"
        echo "  错误: ${sample_name} samtools coverage失败" >> ${LOG_FILE}
        continue
    fi
    
    # 显示coverage文件的前几行用于调试
    echo "  Coverage文件预览:"
    head -n 3 ${coverage_file}
    
    # 检查coverage文件是否有数据
    line_count=$(wc -l < ${coverage_file})
    echo "  Coverage文件行数: ${line_count}"
    
    if [ ${line_count} -le 1 ]; then
        echo "  警告: Coverage文件没有数据"
        echo "  警告: ${sample_name} coverage文件无数据" >> ${LOG_FILE}
        echo -e "${sample_name}\t${total_reads}\t${mapped_reads}\t0\t0\t0\t0" >> ${SUMMARY_FILE}
        continue
    fi
    
    # 提取覆盖度信息
    # samtools coverage输出格式：
    # #rname  startpos  endpos  numreads  covbases  coverage  meandepth  meanbaseq  meanmapq
    
    # 计算总的平均覆盖度（按长度加权）
    avg_cov=$(awk 'NR>1 {
        len = $3 - $2 + 1
        sum += $7 * len
        total_len += len
    } END {
        if(total_len > 0) printf "%.2f", sum/total_len
        else print "0"
    }' ${coverage_file})
    
    # 覆盖的碱基总数
    covered_bases=$(awk 'NR>1 {sum+=$5} END {print sum+0}' ${coverage_file})
    
    # 总碱基数
    total_bases=$(awk 'NR>1 {sum+=($3-$2+1)} END {print sum+0}' ${coverage_file})
    
    # 覆盖率（加权平均）
    coverage_rate=$(awk 'NR>1 {
        len = $3 - $2 + 1
        sum += $6 * len
        total_len += len
    } END {
        if(total_len > 0) printf "%.2f", sum/total_len
        else print "0"
    }' ${coverage_file})
    
    echo "  平均覆盖度: ${avg_cov}x"
    echo "  覆盖的碱基: ${covered_bases}"
    echo "  总碱基数: ${total_bases}"
    echo "  覆盖率: ${coverage_rate}%"
    
    # 写入汇总文件
    echo -e "${sample_name}\t${total_reads}\t${mapped_reads}\t${avg_cov}\t${covered_bases}\t${total_bases}\t${coverage_rate}" >> ${SUMMARY_FILE}
    
    echo "  完成!"
done

echo ""
echo "========================================"
echo "所有文件处理完成!"
echo "========================================"
echo "详细结果保存在: ${OUTPUT_DIR}"
echo "汇总结果保存在: ${SUMMARY_FILE}"
echo "日志文件: ${LOG_FILE}"

# 显示汇总统计
echo ""
echo "========== 覆盖度汇总 =========="
column -t ${SUMMARY_FILE}

echo "" >> ${LOG_FILE}
echo "结束时间: $(date)" >> ${LOG_FILE}