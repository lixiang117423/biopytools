#!/bin/bash

# 从VCF文件中随机抽取10%的SNP
# 作者: Claude
# 日期: 2026-01-03

# 设置路径
INPUT_VCF="/share/org/YZWL/yzwl_lixg/project/21.野生大豆群体/05.所有样品/03.按染色体过滤VCF文件/chr.snp.vcf.gz"
OUTPUT_DIR="/share/org/YZWL/yzwl_lixg/project/21.野生大豆群体/05.所有样品/08.10%的SNP的admixture"
OUTPUT_VCF="${OUTPUT_DIR}/chr.snp.10p.random.vcf.gz"

# 创建输出目录
mkdir -p "${OUTPUT_DIR}"

# 方法1: 使用bcftools (推荐，速度快)
echo "使用bcftools随机抽取10%的SNP..."
bcftools view "${INPUT_VCF}" \
    | bcftools +random \
    --threads 4 \
    --number 20% \
    -Oz -o "${OUTPUT_VCF}"

# 建立索引
tabix -p vcf "${OUTPUT_VCF}"

echo "完成！输出文件: ${OUTPUT_VCF}"
echo ""
echo "统计信息："
echo "原始SNP数量:"
zcat "${INPUT_VCF}" | grep -v "^#" | wc -l
echo "抽取后SNP数量:"
zcat "${OUTPUT_VCF}" | grep -v "^#" | wc -l

# 如果没有bcftools +random插件，可以用方法2（vcftools的--thin参数，但这是按位置均匀抽稀，不是完全随机）
# 注释掉的方法2，仅作为备用
# echo "使用vcftools抽取10%的SNP..."
# vcftools --gzvcf "${INPUT_VCF}" \
#          --thin 0.2 \
#          --recode \
#          --recode-INFO-all \
#          --out "${OUTPUT_DIR}/chr.snp.10p.thin"
#
# bgzip -c ${OUTPUT_DIR}/chr.snp.10p.thin.recode.vcf > ${OUTPUT_DIR}/chr.snp.10p.thin.vcf.gz
# tabix -p vcf ${OUTPUT_DIR}/chr.snp.10p.thin.vcf.gz
