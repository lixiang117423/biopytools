#!/bin/bash

# 设置源目录和目标目录
SOURCE_DIR="/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/31.重测序数据比对到挂载的结果上/02.mapping/bam"
OUTPUT_DIR="$(pwd)"

# 线程数
THREADS=8

# 染色体名称
CHR="Chr12"

# 创建日志文件
LOG_FILE="${OUTPUT_DIR}/extract_Chr12_to_fastq.log"
echo "开始提取Chr12的reads并转换为FASTQ - $(date)" > "$LOG_FILE"

# 计数器
TOTAL=0
SUCCESS=0
FAILED=0

# 遍历源目录中的所有BAM文件
for bamfile in "${SOURCE_DIR}"/*.sorted.bam; do
    if [ -f "$bamfile" ]; then
        # 获取文件名并提取样品名
        filename=$(basename "$bamfile")
        # 从文件名中提取样品名，例如 OV8-277.sorted.bam -> OV8-277
        sample_name=$(echo "$filename" | sed 's/\.sorted\.bam//')
        output_file="${OUTPUT_DIR}/${sample_name}.fq.gz"

        TOTAL=$((TOTAL + 1))

        echo "[$TOTAL] 处理: $filename" | tee -a "$LOG_FILE"
        echo "    样品名: $sample_name" | tee -a "$LOG_FILE"
        echo "    提取染色体: $CHR" | tee -a "$LOG_FILE"
        echo "    输出: $output_file" | tee -a "$LOG_FILE"

        # 检查输出文件是否已存在
        if [ -f "$output_file" ]; then
            echo "    警告: 输出文件已存在，跳过" | tee -a "$LOG_FILE"
            continue
        fi

        # 创建临时文件
        temp_bam="${OUTPUT_DIR}/.${sample_name}_temp.bam"

        # 方法1: 先提取Chr12的reads到临时BAM文件，再转换为fastq
        if samtools view -@ "$THREADS" -b "$bamfile" "$CHR" -o "$temp_bam" 2>> "$LOG_FILE"; then

            # 统计reads数量
            read_count=$(samtools view -c "$temp_bam" 2>> "$LOG_FILE")

            # 转换为fastq并压缩
            if samtools fastq -@ "$THREADS" -n "$temp_bam" | gzip -c > "$output_file" 2>> "$LOG_FILE"; then

                # 获取文件大小
                file_size=$(du -h "$output_file" | cut -f1)

                SUCCESS=$((SUCCESS + 1))
                echo "    成功! Chr12 reads数: $read_count, 文件大小: $file_size" | tee -a "$LOG_FILE"
            else
                FAILED=$((FAILED + 1))
                echo "    失败! (转换fastq时出错)" | tee -a "$LOG_FILE"
                rm -f "$output_file"
            fi

            # 删除临时文件
            rm -f "$temp_bam"
        else
            FAILED=$((FAILED + 1))
            echo "    失败! (提取Chr12时出错)" | tee -a "$LOG_FILE"
        fi

        echo "----------------------------------------" | tee -a "$LOG_FILE"
    fi
done

# 输出统计信息
echo "" | tee -a "$LOG_FILE"
echo "转换完成 - $(date)" | tee -a "$LOG_FILE"
echo "总计文件数: $TOTAL" | tee -a "$LOG_FILE"
echo "成功: $SUCCESS" | tee -a "$LOG_FILE"
echo "失败: $FAILED" | tee -a "$LOG_FILE"

# 检查是否没有找到Chr12
if [ $SUCCESS -eq 0 ] && [ $TOTAL -gt 0 ]; then
    echo "" | tee -a "$LOG_FILE"
    echo "警告: 所有文件都提取失败!" | tee -a "$LOG_FILE"
    echo "可能的原因:" | tee -a "$LOG_FILE"
    echo "1. 染色体名称可能不是'$CHR',而是其他格式" | tee -a "$LOG_FILE"
    echo "" | tee -a "$LOG_FILE"
    echo "检查第一个BAM文件中的染色体名称:" | tee -a "$LOG_FILE"
    first_bam=$(ls "${SOURCE_DIR}"/*.sorted.bam | head -1)
    samtools view -H "$first_bam" | grep "^@SQ" | tee -a "$LOG_FILE"
fi
