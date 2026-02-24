#!/bin/bash

# --- 用户可配置变量 ---
SNP_VCF="variation.filtered.snp.vcf.gz"      # 输入的SNP VCF文件名
INDEL_VCF="variation.filtered.indel.vcf.gz"    # 输入的INDEL VCF文件名
POP1_PREFIX="BS"                             # 第一个群体的样本名前缀
POP2_PREFIX="GC"                             # 第二个群体的样本名前缀
WINDOW_SIZE=100000                           # 窗口大小 (100kb)

# --- 脚本主体 ---

# 检查VCF文件是否存在
if [ ! -f "$SNP_VCF" ]; then
    echo "错误: SNP VCF文件 '$SNP_VCF' 不存在。"
    exit 1
fi

if [ ! -f "$INDEL_VCF" ]; then
    echo "错误: INDEL VCF文件 '$INDEL_VCF' 不存在。"
    exit 1
fi

echo "步骤 1: 创建群体样本列表 (从SNP VCF文件中读取)"
# 从VCF文件中提取所有样本名，并根据前缀分组
# 假设SNP和INDEL VCF文件中的样本名是一致的，我们只需要创建一次列表
bcftools query -l ${SNP_VCF} | grep "^${POP1_PREFIX}" > ${POP1_PREFIX}_samples.txt
bcftools query -l ${SNP_VCF} | grep "^${POP2_PREFIX}" > ${POP2_PREFIX}_samples.txt

# 检查样本列表是否成功创建
if [ ! -s "${POP1_PREFIX}_samples.txt" ] || [ ! -s "${POP2_PREFIX}_samples.txt" ]; then
    echo "错误: 未能根据提供的前缀 '${POP1_PREFIX}' 和 '${POP2_PREFIX}' 创建一个或两个样本列表。"
    echo "请检查您的VCF文件中的样本名是否正确。"
    exit 1
fi

echo "为群体 '${POP1_PREFIX}' 创建了样本列表: ${POP1_PREFIX}_samples.txt"
echo "为群体 '${POP2_PREFIX}' 创建了样本列表: ${POP2_PREFIX}_samples.txt"
echo ""

# --- 计算SNP密度 ---
echo "步骤 2: 为群体 '${POP1_PREFIX}' 计算SNP密度"
vcftools --gzvcf ${SNP_VCF} \
    --keep ${POP1_PREFIX}_samples.txt \
    --SNPdensity ${WINDOW_SIZE} \
    --out ${POP1_PREFIX}_snp

echo "群体 '${POP1_PREFIX}' 的SNP密度计算完成。结果保存在 '${POP1_PREFIX}_snp.snpden' 文件中。"
echo ""

echo "步骤 3: 为群体 '${POP2_PREFIX}' 计算SNP密度"
vcftools --gzvcf ${SNP_VCF} \
    --keep ${POP2_PREFIX}_samples.txt \
    --SNPdensity ${WINDOW_SIZE} \
    --out ${POP2_PREFIX}_snp

echo "群体 '${POP2_PREFIX}' 的SNP密度计算完成。结果保存在 '${POP2_PREFIX}_snp.snpden' 文件中。"
echo ""


# --- 计算INDEL密度 ---
# VCFtools 没有专门的 --INDELdensity 选项，但 --site-density 会计算文件中的所有变异类型。
# 因为我们提供了只包含INDEL的VCF文件，所以这等同于计算INDEL密度。
echo "步骤 4: 为群体 '${POP1_PREFIX}' 计算INDEL密度"
vcftools --gzvcf ${INDEL_VCF} \
    --keep ${POP1_PREFIX}_samples.txt \
    --site-density ${WINDOW_SIZE} \
    --out ${POP1_PREFIX}_indel

echo "群体 '${POP1_PREFIX}' 的INDEL密度计算完成。结果保存在 '${POP1_PREFIX}_indel.siteden' 文件中。"
echo ""

echo "步骤 5: 为群体 '${POP2_PREFIX}' 计算INDEL密度"
vcftools --gzvcf ${INDEL_VCF} \
    --keep ${POP2_PREFIX}_samples.txt \
    --site-density ${WINDOW_SIZE} \
    --out ${POP2_PREFIX}_indel

echo "群体 '${POP2_PREFIX}' 的INDEL密度计算完成。结果保存在 '${POP2_PREFIX}_indel.siteden' 文件中。"
echo ""

echo "所有计算已完成！"
