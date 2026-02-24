#!/bin/bash

# ================= 配置区域 =================
INPUT_VCF="variation.filtered.merged.vcf.gz"
GROUP1_TAG="BS"
GROUP2_TAG="GC"
FST_THRESHOLD=0.5
OUT_PREFIX="Fst0.5_PAV_Matrix"
# ===========================================

echo ">>> Step 1: 准备分组并计算 Fst..."

# 1. 提取样品名并分组
zgrep -m 1 "^#CHROM" $INPUT_VCF | cut -f 10- | tr '\t' '\n' > all_samples.txt
grep "^$GROUP1_TAG" all_samples.txt > group_BS.txt
grep "^$GROUP2_TAG" all_samples.txt > group_GC.txt

# 2. 计算 Fst (如果文件不存在)
if [ ! -f "${OUT_PREFIX}.weir.fst" ]; then
    vcftools --gzvcf $INPUT_VCF \
             --weir-fst-pop group_BS.txt \
             --weir-fst-pop group_GC.txt \
             --out $OUT_PREFIX
fi

echo ">>> Step 2: 筛选 Fst > $FST_THRESHOLD 的位点..."

# 筛选 Fst > 0.5 的行，提取 CHROM, POS, FST 三列
# 结果保存为临时文件 fst_info.tmp
awk -v th="$FST_THRESHOLD" 'NR>1 && $3 != "-nan" && $3 != "nan" && $3 > th {print $1 "\t" $2 "\t" $3}' ${OUT_PREFIX}.weir.fst > fst_info.tmp

# 提取只有 CHROM 和 POS 的列表，用于给 VCFtools 传参
cut -f 1,2 fst_info.tmp > fst_positions.tmp

NUM=$(wc -l < fst_positions.tmp)
echo "    共筛选出 $NUM 个显著差异位点。"

if [ "$NUM" -eq 0 ]; then
    echo "Error: 没有位点满足条件，请降低阈值。"
    rm *.tmp
    exit 1
fi

echo ">>> Step 3: 提取基因型并转换为 0/1 格式..."

# 1. 提取原始 GT 矩阵 (格式如 0/0, 0/1...)
vcftools --gzvcf $INPUT_VCF \
         --positions fst_positions.tmp \
         --extract-FORMAT-pg GT \
         --out raw_genotypes

# 2. 使用 sed 进行批量替换，生成 0/1 矩阵
# 逻辑：
# \t0/0 -> \t0 (Ref)
# \t1/1 -> \t1 (Alt)
# \t0/1 -> \t1 (Het -> 1, 视为有变异)
# \t1/0 -> \t1 (Het -> 1)
# \t./. -> \tNA (Missing)
# 注意：这里假设是双等位基因位点(Biallelic)，这是Fst分析最常见的情况
sed 's#\t0/0#\t0#g; s#\t0|0#\t0#g; s#\t1/1#\t1#g; s#\t1|1#\t1#g; s#\t0/1#\t1#g; s#\t0|1#\t1#g; s#\t1/0#\t1#g; s#\t1|0#\t1#g; s#\t\./\.#\tNA#g' raw_genotypes.GT.FORMAT > genotypes_01.tmp

echo ">>> Step 4: 合并 Fst 数值和 0/1 矩阵..."

# 此时我们有两个文件：
# 1. fst_info.tmp: [CHROM] [POS] [FST_VALUE]
# 2. genotypes_01.tmp: [CHROM] [POS] [Sample1] [Sample2]...

# 为了合并，我们需要去掉 genotypes_01.tmp 的前两列(CHROM和POS)，因为 fst_info.tmp 里面已经有了
# 但是要注意：fst_info.tmp 没有表头，而 genotypes_01.tmp 有表头。

# 4.1 处理表头
# 提取 genotypes_01.tmp 的第一行（样品名），去掉前两列，保留样品名
head -n 1 genotypes_01.tmp | cut -f 3- > header_samples.tmp
# 创建总表头: CHROM POS FST_VALUE [Sample_Names...]
echo -e "CHROM\tPOS\tFST_VALUE\t$(cat header_samples.tmp)" > ${OUT_PREFIX}_Final_Merged.txt

# 4.2 处理数据体
# 提取 genotype 的数据部分（去掉第1行表头，去掉前2列位置信息）
tail -n +2 genotypes_01.tmp | cut -f 3- > body_genotypes.tmp

# 4.3 将 [Fst信息] 与 [0/1矩阵] 左右拼接 (paste)
paste fst_info.tmp body_genotypes.tmp >> ${OUT_PREFIX}_Final_Merged.txt

# 清理临时文件
rm *.tmp raw_genotypes.GT.FORMAT raw_genotypes.log

echo "========================================================"
echo "处理完成！"
echo "最终输出文件: ${OUT_PREFIX}_Final_Merged.txt"
echo "格式说明: CHROM  POS  FST值  样品1(0/1)  样品2(0/1) ..."
echo "========================================================"
# 预览前5行
head -n 5 ${OUT_PREFIX}_Final_Merged.txt | cut -f 1-6
echo "..."
