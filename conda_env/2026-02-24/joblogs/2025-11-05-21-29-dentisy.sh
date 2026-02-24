#!/bin/bash

# --- 用户可配置变量 ---
SNP_VCF="variation.filtered.snp.vcf.gz"
INDEL_VCF="variation.filtered.indel.vcf.gz"
POP1_PREFIX="BS"
POP2_PREFIX="GC"
WINDOW_SIZE=100000

# --- 脚本主体 ---

# 检查VCF文件是否存在
# ... (这部分和原脚本一样，此处省略) ...

echo "步骤 1: 创建群体样本列表"
bcftools query -l ${SNP_VCF} | grep "^${POP1_PREFIX}" > ${POP1_PREFIX}_samples.txt
bcftools query -l ${SNP_VCF} | grep "^${POP2_PREFIX}" > ${POP2_PREFIX}_samples.txt
echo "样本列表创建完成。"
echo ""

# --- 核心修正步骤：为每个群体生成过滤后的VCF文件 ---

echo "步骤 2: 为群体 '${POP1_PREFIX}' 创建只包含其多态性位点的VCF"
vcftools --gzvcf ${SNP_VCF} \
    --keep ${POP1_PREFIX}_samples.txt \
    --mac 1 \
    --recode --recode-INFO-all \
    --out ${POP1_PREFIX}_poly_snps
# mv ${POP1_PREFIX}_poly_snps.recode.vcf ${POP1_PREFIX}_poly_snps.vcf # 可选：重命名
# bgzip ${POP1_PREFIX}_poly_snps.vcf && tabix -p vcf ${POP1_PREFIX}_poly_snps.vcf.gz # 可选：压缩并索引

echo "步骤 3: 为群体 '${POP2_PREFIX}' 创建只包含其多态性位点的VCF"
vcftools --gzvcf ${SNP_VCF} \
    --keep ${POP2_PREFIX}_samples.txt \
    --mac 1 \
    --recode --recode-INFO-all \
    --out ${POP2_PREFIX}_poly_snps
# mv ${POP2_PREFIX}_poly_snps.recode.vcf ${POP2_PREFIX}_poly_snps.vcf # 可选
# bgzip ${POP2_PREFIX}_poly_snps.vcf && tabix -p vcf ${POP2_PREFIX}_poly_snps.vcf.gz # 可选

# 对INDEL文件做同样处理
echo "步骤 4: 为群体 '${POP1_PREFIX}' 创建只包含其多态性位点的INDEL VCF"
vcftools --gzvcf ${INDEL_VCF} \
    --keep ${POP1_PREFIX}_samples.txt \
    --mac 1 \
    --recode --recode-INFO-all \
    --out ${POP1_PREFIX}_poly_indels

echo "步骤 5: 为群体 '${POP2_PREFIX}' 创建只包含其多态性位点的INDEL VCF"
vcftools --gzvcf ${INDEL_VCF} \
    --keep ${POP2_PREFIX}_samples.txt \
    --mac 1 \
    --recode --recode-INFO-all \
    --out ${POP2_PREFIX}_poly_indels

echo ""
echo "特定群体的VCF文件已生成，现在开始计算密度..."
echo ""

# --- 步骤 6: 在新的、过滤后的VCF文件上计算密度 ---

# 计算SNP密度
echo "计算 '${POP1_PREFIX}' 群体的SNP密度..."
vcftools --vcf ${POP1_PREFIX}_poly_snps.recode.vcf \
    --SNPdensity ${WINDOW_SIZE} \
    --out ${POP1_PREFIX}_snp_filtered

echo "计算 '${POP2_PREFIX}' 群体的SNP密度..."
vcftools --vcf ${POP2_PREFIX}_poly_snps.recode.vcf \
    --SNPdensity ${WINDOW_SIZE} \
    --out ${POP2_PREFIX}_snp_filtered

# 计算INDEL密度
echo "计算 '${POP1_PREFIX}' 群体的INDEL密度..."
vcftools --vcf ${POP1_PREFIX}_poly_indels.recode.vcf \
    --site-density ${WINDOW_SIZE} \
    --out ${POP1_PREFIX}_indel_filtered

echo "计算 '${POP2_PREFIX}' 群体的INDEL密度..."
vcftools --vcf ${POP2_PREFIX}_poly_indels.recode.vcf \
    --site-density ${WINDOW_SIZE} \
    --out ${POP2_PREFIX}_indel_filtered

echo ""
echo "所有计算已完成！现在的结果应该是不同的了。"
echo "请检查后缀为 '_filtered.snpden' 和 '_filtered.siteden' 的文件。"

# 清理中间文件（可选）
# rm *_poly_snps.recode.vcf *_poly_indels.recode.vcf *_poly_snps.log *_poly_indels.log
