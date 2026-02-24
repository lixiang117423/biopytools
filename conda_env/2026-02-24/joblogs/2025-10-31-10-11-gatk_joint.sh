#!/bin/bash

# --- 配置 ---
# 设置你的GVCF文件所在的目录
GVCF_DIR="/share/org/YZWL/yzwl_lixg/project/96.yabian/14.新的变异信息/01.gvcf"
# --- 配置结束 ---

# 1. 检查 tabix 命令是否存在
if ! command -v tabix &> /dev/null; then
    echo "错误：tabix 命令未找到。请先安装 htslib/samtools。"
    exit 1
fi

# 2. 检查目录是否存在
if [ ! -d "$GVCF_DIR" ]; then
    echo "错误：目录不存在: $GVCF_DIR"
    exit 1
fi

echo "============================================="
echo "开始扫描并创建GVCF索引"
echo "目标目录: $GVCF_DIR"
echo "============================================="

# 初始化计数器
total_files=0
indexed_files=0
skipped_files=0
failed_files=0

# 3. 遍历目录下所有以 .g.vcf.gz 结尾的文件
for gvcf_file in "$GVCF_DIR"/*.g.vcf.gz; do
    
    # 检查是否找到了文件，如果没有匹配的文件，循环会把 "*.g.vcf.gz" 作为文件名
    if [ ! -e "$gvcf_file" ]; then
        continue
    fi
    
    ((total_files++))
    index_file="${gvcf_file}.tbi"
    
    echo -n "处理中: $(basename "$gvcf_file") ... "

    # 4. 检查索引是否已存在
    if [ -f "$index_file" ]; then
        echo "已跳过 (索引已存在)"
        ((skipped_files++))
    else
        # 5. 创建索引
        tabix -p vcf "$gvcf_file"
        
        # 检查命令执行结果
        if [ $? -eq 0 ]; then
            echo "索引创建成功"
            ((indexed_files++))
        else
            echo "索引创建失败！"
            ((failed_files++))
        fi
    fi
done

echo "============================================="
echo "处理完成！"
echo "总共扫描文件: $total_files"
echo "成功创建索引: $indexed_files"
echo "跳过 (已存在): $skipped_files"
echo "失败: $failed_files"
echo "============================================="

# 如果有失败，返回非零退出码
if [ $failed_files -gt 0 ]; then
    exit 1
fi


biopytools gatk-joint \
    -i /share/org/YZWL/yzwl_lixg/project/96.yabian/14.新的变异信息/01.gvcf \
    -r /share/org/YZWL/yzwl_lixg/project/96.yabian/01.data/genome/genome.fa \
    -o /share/org/YZWL/yzwl_lixg/project/96.yabian/14.新的变异信息/02.gatk \
    -t 80 -m 200g
