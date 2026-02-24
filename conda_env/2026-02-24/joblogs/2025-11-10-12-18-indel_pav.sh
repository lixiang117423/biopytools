#!/bin/bash

# --- 设置变量 ---

# 输入的 Indel VCF 文件
VCF_FILE="variation.filtered.indel.vcf.gz"

# 输出的 0/1 矩阵文件名
OUTPUT_MATRIX="indel_presence_absence_matrix.01.tsv"


# --- 脚本开始 ---
echo "正在为 Indel 数据生成 0/1 存在/缺失矩阵..."
echo "输入 VCF 文件: ${VCF_FILE}"

# --- 步骤 1: 创建矩阵的表头 ---
# 首先，提取 VCF 中的样本名，并将它们格式化为制表符分隔的一行
# 然后，在样本名前加上固定的变异信息列名
# 最后，将完整的表头写入输出文件
(echo -n -e "CHROM\tPOS\tREF\tALT\t" && \
 bcftools query -l ${VCF_FILE} | tr '\n' '\t' | sed 's/\t$//g'
) > ${OUTPUT_MATRIX}


# --- 步骤 2: 生成 0/1 矩阵主体内容 ---
# 使用 bcftools query 提取每个位点的基因型 (GT)
# 然后通过管道(|)将输出实时传递给 awk 命令进行处理
# awk 会检查每个样本的基因型，并将其转换为 0 或 1
bcftools query -f '\n%CHROM\t%POS\t%REF\t%ALT[\t%GT]' ${VCF_FILE} | \
awk '
BEGIN { OFS="\t" } # 设置输出字段分隔符为制表符
{
    # 打印前四列 (CHROM, POS, REF, ALT)
    printf "%s\t%s\t%s\t%s", $1, $2, $3, $4;

    # 从第五列开始循环，处理每个样本的基因型
    for (i=5; i<=NF; i++) {
        # 如果基因型为杂合(0/1)或纯合变异(1/1)，则该变异存在，输出 1
        # 同时涵盖了有相位信息(phased)的基因型，如 0|1, 1|0, 1|1
        if ($i == "0/1" || $i == "1/1" || $i ~ /\|/ && $i != "0|0") {
            printf "\t1";
        } else {
        # 否则 (基因型为 0/0 或 ./.)，则该变异不存在，输出 0
            printf "\t0";
        }
    }
    # 每处理完 VCF 的一行后换行
    printf "\n";
}
' >> ${OUTPUT_MATRIX} # 将 awk 处理后的结果追加到输出文件中

echo "处理完成！"
echo "0/1 矩阵已保存到: ${OUTPUT_MATRIX}"
