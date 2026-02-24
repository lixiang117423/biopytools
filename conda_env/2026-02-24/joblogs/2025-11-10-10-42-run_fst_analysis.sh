#!/bin/bash

# --- 设置变量 ---
# 输入的VCF文件
VCF_FILE="variation.filtered.snp.vcf.gz"

# 群体1的前缀
POP1_PREFIX="BS"

# 群体2的前缀
POP2_PREFIX="GC"

# 输出文件的前缀
OUT_PREFIX="my_analysis"

# 滑窗FST的参数
WINDOW_SIZE=100000  # 窗口大小，例如 100kb
STEP_SIZE=10000     # 步长，例如 10kb

# --- 脚本开始 ---
echo "FST 分析开始..."
echo "输入 VCF 文件: ${VCF_FILE}"

# --- 步骤 1: 创建群体文件 ---
# 使用 zcat 和 bcftools query（或 vcftools --list-individuals）来获取样本列表
# 这里使用 zcat 和 grep 更为直接
echo "步骤 1: 从VCF文件头部创建群体文件..."

# 从 VCF 文件头中提取所有样本名，然后用 grep 筛选
zcat ${VCF_FILE} | grep -m1 '^#CHROM' | cut -f 10- | tr '\t' '\n' > all_samples.txt

# 创建群体1的文件 (BS 开头的样本)
grep "^${POP1_PREFIX}" all_samples.txt > pop1.txt
echo "群体1 (${POP1_PREFIX}) 包含 $(wc -l < pop1.txt) 个样本，列表保存在 pop1.txt"

# 创建群体2的文件 (GC 开头的样本)
grep "^${POP2_PREFIX}" all_samples.txt > pop2.txt
echo "群体2 (${POP2_PREFIX}) 包含 $(wc -l < pop2.txt) 个样本，列表保存在 pop2.txt"

# 检查群体文件是否为空
if [ ! -s pop1.txt ] || [ ! -s pop2.txt ]; then
    echo "错误：一个或两个群体文件为空！请检查您的 VCF 文件和群体前缀是否正确。"
    exit 1
fi


# --- 步骤 2: 计算每个位点的 FST ---
echo "步骤 2: 计算每个位点的 Weir and Cockerham's FST..."
vcftools --gzvcf ${VCF_FILE} \
    --weir-fst-pop pop1.txt \
    --weir-fst-pop pop2.txt \
    --out ${OUT_PREFIX}

echo "逐位点 FST 计算完成。结果保存在 ${OUT_PREFIX}.weir.fst 文件中。"


# --- 步骤 3: 计算滑窗 FST ---
echo "步骤 3: 计算滑窗 FST (窗口大小: ${WINDOW_SIZE}bp, 步长: ${STEP_SIZE}bp)..."
vcftools --gzvcf ${VCF_FILE} \
    --weir-fst-pop pop1.txt \
    --weir-fst-pop pop2.txt \
    --fst-window-size ${WINDOW_SIZE} \
    --fst-window-step ${STEP_SIZE} \
    --out ${OUT_PREFIX}_windowed

echo "滑窗 FST 计算完成。结果保存在 ${OUT_PREFIX}_windowed.windowed.weir.fst 文件中。"


# --- 清理临时文件 ---
rm all_samples.txt
# 如果您想保留群体文件，可以注释掉下面两行
# rm pop1.txt
# rm pop2.txt

echo "分析流程结束！"
