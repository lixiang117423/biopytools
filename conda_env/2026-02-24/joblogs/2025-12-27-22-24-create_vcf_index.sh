#!/bin/bash

# VCF文件批量创建索引脚本
# 使用tabix为vcf.gz文件创建索引

# VCF文件所在目录
VCF_DIR="/share/org/YZWL/yzwl_lixg/project/21.野生大豆群体/02.mapping/vcf"

# 线程数
THREADS=64

# 检查tabix是否可用
if ! command -v tabix &> /dev/null; then
    echo "错误: tabix未安装或不在PATH中"
    exit 1
fi

# 切换到vcf目录
cd "$VCF_DIR" || exit 1

echo "开始为VCF文件创建索引..."
echo "目录: $VCF_DIR"
echo "线程数: $THREADS"
echo "----------------------------------------"

# 计数器
count=0
total=0

# 统计总文件数
total=$(ls *.g.vcf.gz 2>/dev/null | wc -l)

if [ $total -eq 0 ]; then
    echo "错误: 未找到 .g.vcf.gz 文件"
    exit 1
fi

echo "共找到 $total 个VCF文件"
echo "----------------------------------------"

# 遍历所有vcf.gz文件
for vcf_file in *.g.vcf.gz; do
    # 检查索引是否已存在
    if [ -f "${vcf_file}.tbi" ]; then
        echo "[$((count+1))/$total] 索引已存在: $vcf_file"
        ((count++))
        continue
    fi

    echo "[$((count+1))/$total] 正在创建索引: $vcf_file"

    # 使用tabix创建索引
    if tabix -p vcf -@ $THREADS "$vcf_file"; then
        echo "  ✓ 索引创建成功: ${vcf_file}.tbi"
    else
        echo "  ✗ 索引创建失败: $vcf_file"
    fi

    ((count++))
done

echo "----------------------------------------"
echo "索引创建完成！"
