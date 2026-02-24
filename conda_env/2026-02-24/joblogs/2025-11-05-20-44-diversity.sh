#!/bin/bash

# --- 用户可配置变量 ---
VCF_FILE="variation.filtered.snp.vcf.gz"  # 输入的VCF文件名
POP1_PREFIX="BS"                         # 第一个群体的样本名前缀
POP2_PREFIX="GC"                         # 第二个群体的样本名前缀
WINDOW_SIZE=500000                       # 滑窗大小 (500kb)
STEP_SIZE=100000                         # 步长 (100kb)

# --- 脚本主体 ---

# 检查VCF文件是否存在
if [ ! -f "$VCF_FILE" ]; then
    echo "错误: VCF文件 '$VCF_FILE' 不存在。"
    exit 1
fi

echo "步骤 1: 创建群体样本列表"
# 从VCF文件中提取所有样本名，并根据前缀分组
bcftools query -l ${VCF_FILE} | grep "^${POP1_PREFIX}" > ${POP1_PREFIX}_samples.txt
bcftools query -l ${VCF_FILE} | grep "^${POP2_PREFIX}" > ${POP2_PREFIX}_samples.txt

# 检查样本列表是否成功创建
if [ ! -s "${POP1_PREFIX}_samples.txt" ] || [ ! -s "${POP2_PREFIX}_samples.txt" ]; then
    echo "错误: 未能根据提供的前缀 '${POP1_PREFIX}' 和 '${POP2_PREFIX}' 创建一个或两个样本列表。"
    echo "请检查您的VCF文件中的样本名是否正确。"
    exit 1
fi

echo "为群体 '${POP1_PREFIX}' 创建了样本列表: ${POP1_PREFIX}_samples.txt"
echo "为群体 '${POP2_PREFIX}' 创建了样本列表: ${POP2_PREFIX}_samples.txt"
echo ""

# --- 计算核苷酸多样性 (π) ---
echo "步骤 2: 为群体 '${POP1_PREFIX}' 计算核苷酸多样性 (π)"
vcftools --gzvcf ${VCF_FILE} \
    --keep ${POP1_PREFIX}_samples.txt \
    --window-pi ${WINDOW_SIZE} \
    --window-pi-step ${STEP_SIZE} \
    --out ${POP1_PREFIX}

echo "群体 '${POP1_PREFIX}' 的 π 计算完成。结果保存在 '${POP1_PREFIX}.windowed.pi' 文件中。"
echo ""

echo "步骤 3: 为群体 '${POP2_PREFIX}' 计算核苷酸多样性 (π)"
vcftools --gzvcf ${VCF_FILE} \
    --keep ${POP2_PREFIX}_samples.txt \
    --window-pi ${WINDOW_SIZE} \
    --window-pi-step ${STEP_SIZE} \
    --out ${POP2_PREFIX}

echo "群体 '${POP2_PREFIX}' 的 π 计算完成。结果保存在 '${POP2_PREFIX}.windowed.pi' 文件中。"
echo ""

# --- 计算Tajima's D ---
# 注意: VCFtools的Tajima's D计算是基于固定大小的窗口（bins），而不是滑窗。
# 这里我们将窗口大小设置为与您的步长一致，这是一个常见的做法。
echo "步骤 4: 为群体 '${POP1_PREFIX}' 计算Tajima's D (窗口大小 ${STEP_SIZE}bp)"
vcftools --gzvcf ${VCF_FILE} \
    --keep ${POP1_PREFIX}_samples.txt \
    --TajimaD ${STEP_SIZE} \
    --out ${POP1_PREFIX}

echo "群体 '${POP1_PREFIX}' 的 Tajima's D 计算完成。结果保存在 '${POP1_PREFIX}.Tajima.D' 文件中。"
echo ""

echo "步骤 5: 为群体 '${POP2_PREFIX}' 计算Tajima's D (窗口大小 ${STEP_SIZE}bp)"
vcftools --gzvcf ${VCF_FILE} \
    --keep ${POP2_PREFIX}_samples.txt \
    --TajimaD ${STEP_SIZE} \
    --out ${POP2_PREFIX}

echo "群体 '${POP2_PREFIX}' 的 Tajima's D 计算完成。结果保存在 '${POP2_PREFIX}.Tajima.D' 文件中。"
echo ""

echo "所有计算已完成！"
