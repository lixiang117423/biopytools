#!/bin/bash

# 遇到错误时立即退出
set -e

# --- 配置参数 ---
SRA_DIR="/share/org/YZWL/yzwl_lixg/project/95.gdl_metagenome/01.data/raw/sra_data"
OUTPUT_DIR="/share/org/YZWL/yzwl_lixg/project/95.gdl_metagenome/01.data/raw/fastq_data"
# --- 配置结束 ---

# 检查 fastq-dump 命令
if ! command -v fastq-dump &> /dev/null
then
    echo "错误: 'fastq-dump' 命令未找到。请安装 SRA Toolkit。"
    exit 1
fi

# 创建输出目录
echo "创建输出目录 (如果不存在): $OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR"

echo "开始处理 SRA 文件..."

# 计数器
file_count=0

# 循环处理指定目录下所有以 SRR 开头的文件
# 修改点：从 "$SRA_DIR"/*.sra 改为 "$SRA_DIR"/SRR*
for sra_file in "$SRA_DIR"/SRR*
do
  # 检查匹配到的是否为文件
  if [ -f "$sra_file" ]; then
    sra_filename=$(basename "$sra_file")
    
    echo "----------------------------------------"
    echo "正在处理: $sra_filename"
    
    # 执行 fastq-dump
    fastq-dump --split-files --gzip -O "$OUTPUT_DIR" "$sra_file"
    
    if [ $? -eq 0 ]; then
      echo "成功拆分 $sra_filename"
      file_count=$((file_count + 1))
    else
      echo "警告: 处理 $sra_filename 时遇到错误。"
    fi
  fi
done

echo "----------------------------------------"
if [ "$file_count" -eq 0 ]; then
  echo "警告: 未找到任何以 'SRR' 开头的文件进行处理。"
else
  echo "所有文件处理完成！共处理了 $file_count 个文件。"
fi
