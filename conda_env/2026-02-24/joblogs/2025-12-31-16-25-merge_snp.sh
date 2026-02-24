#!/bin/bash

# 合并野生大豆群体VCF文件脚本
# 用法: bash merge_vcf.sh

# 设置基础路径
BASE_DIR="/share/org/YZWL/yzwl_lixg/project/21.野生大豆群体/05.所有样品/03.按染色体过滤VCF文件"

# 设置输出目录和文件名
OUTPUT_DIR="${BASE_DIR}"
mkdir -p ${OUTPUT_DIR}
OUTPUT_FILE="${OUTPUT_DIR}/chr.snp.vcf.gz"

# 创建临时文件列表
VCF_LIST="${OUTPUT_DIR}/vcf_list_snp.txt"

echo "开始准备VCF文件列表..."

# 生成VCF文件列表
for i in $(seq -f "%02g" 1 20); do
    vcf_file="${BASE_DIR}/chr${i}/variation.filtered.snp.biallelic.vcf.gz"
    if [ -f "$vcf_file" ]; then
        echo "$vcf_file" >> ${VCF_LIST}
        echo "找到文件: chr${i}"
    else
        echo "警告: 文件不存在 - chr${i}"
    fi
done

# 检查是否有文件要合并
if [ ! -s ${VCF_LIST} ]; then
    echo "错误: 没有找到任何VCF文件!"
    exit 1
fi

echo ""
echo "找到 $(wc -l < ${VCF_LIST}) 个VCF文件"
echo ""
echo "开始合并VCF文件..."
echo "输出文件: ${OUTPUT_FILE}"
echo ""

# 使用bcftools合并VCF文件
# -O z: 输出压缩的VCF格式
# --threads: 使用多线程加速
bcftools concat \
    -O z \
    --threads 64 \
    -f ${VCF_LIST} \
    -o ${OUTPUT_FILE}

# 检查合并是否成功
if [ $? -eq 0 ]; then
    echo ""
    echo "VCF文件合并成功!"
    echo ""
    echo "创建索引文件..."
    # 创建索引
    bcftools index -t ${OUTPUT_FILE}
    
    if [ $? -eq 0 ]; then
        echo "索引创建成功!"
        echo ""
        echo "输出文件信息:"
        ls -lh ${OUTPUT_FILE}*
        echo ""
        echo "文件统计:"
        bcftools stats ${OUTPUT_FILE} | grep "^SN"
    else
        echo "警告: 索引创建失败"
    fi
else
    echo ""
    echo "错误: VCF文件合并失败!"
    exit 1
fi

# 清理临时文件
rm -f ${VCF_LIST}

echo ""
echo "完成!"