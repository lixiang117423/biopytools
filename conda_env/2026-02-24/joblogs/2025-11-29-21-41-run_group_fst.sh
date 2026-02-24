#!/bin/bash

# ================= 设置 =================
INPUT_VCF="variation.filtered.merged.vcf.gz"
OUT_PREFIX="BS_vs_GC_all_sites"
GROUP1_TAG="BS"
GROUP2_TAG="GC"
# =======================================

echo ">>> Step 1: 准备分组文件..."

# 1. 获取所有样品名
zgrep -m 1 "^#CHROM" $INPUT_VCF | cut -f 10- | tr '\t' '\n' > all_samples.txt

# 2. 根据前缀生成分组列表
grep "^$GROUP1_TAG" all_samples.txt > group_BS.txt
grep "^$GROUP2_TAG" all_samples.txt > group_GC.txt

# 检查分组数量
N_BS=$(wc -l < group_BS.txt)
N_GC=$(wc -l < group_GC.txt)
echo "    - BS 组样品数: $N_BS"
echo "    - GC 组样品数: $N_GC"

if [[ $N_BS -eq 0 || $N_GC -eq 0 ]]; then
    echo "!!! 错误: 没找到对应的样品，请检查 VCF 头部的样品名是否以 BS 或 GC 开头。"
    exit 1
fi

echo ">>> Step 2: 运行 VCFtools 计算全基因组 Fst (包含 SNP 和 INDEL)..."

# --gzvcf: 读取 gz 格式 VCF
# --weir-fst-pop: 指定两个群体
# --out: 输出文件名
# 注意：这里没有加 --remove-indels，所以 INDEL 也会被计算
vcftools --gzvcf $INPUT_VCF \
         --weir-fst-pop group_BS.txt \
         --weir-fst-pop group_GC.txt \
         --out $OUT_PREFIX

echo ">>> Step 3: 整理结果..."

# 原始结果文件名为 OUT_PREFIX.weir.fst
# 这里我们做一个简单的清洗：
# 1. 去掉 NaN 的行 (NaN通常是因为某个位点在某个群体中完全缺失或无变异，无法计算Fst)
# 2. 加上 gzip 压缩，方便传输和保存

grep -v "NaN" ${OUT_PREFIX}.weir.fst | gzip > ${OUT_PREFIX}.clean.fst.gz

echo "========================================================"
echo "分析完成！"
echo ""
echo "原始输出文件 (含表头和NaN): ${OUT_PREFIX}.weir.fst"
echo "清洗压缩文件 (无NaN,可直接读取): ${OUT_PREFIX}.clean.fst.gz"
echo ""
echo "结果文件包含三列: CHROM(染色体)  POS(位置)  WEIR_AND_COCKERHAM_FST(Fst值)"
echo "========================================================"
