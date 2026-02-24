#!/bin/bash

# 设置输入输出路径
VCF_FILE="/share/org/YZWL/yzwl_lixg/project/21.野生大豆群体/05.所有样品/03.按染色体过滤VCF文件/chr.snp.vcf.gz"
OUT_DIR="/share/org/YZWL/yzwl_lixg/project/21.野生大豆群体/05.所有样品/15.所有样品SNP数量统计"

# 创建输出文件夹
mkdir -p ${OUT_DIR}

# 输出文件
OUT_FILE="${OUT_DIR}/sample_snp_counts.txt"
SUMMARY_FILE="${OUT_DIR}/summary_statistics.txt"

# 检查输入文件是否存在
if [ ! -f "${VCF_FILE}" ]; then
    echo "错误: VCF文件不存在: ${VCF_FILE}"
    exit 1
fi

echo "开始统计每个样品的有效SNP数量..."
echo "VCF文件: ${VCF_FILE}"

# 提取样品列表
echo "提取样品列表..."
bcftools query -l ${VCF_FILE} > ${OUT_DIR}/sample_list.txt
SAMPLE_COUNT=$(wc -l < ${OUT_DIR}/sample_list.txt)
echo "共有 ${SAMPLE_COUNT} 个样品"

# 初始化输出文件
echo -e "Sample\tTotal_SNPs\tHet_SNPs\tHom_Alt_SNPs\tMissing_SNPs" > ${OUT_FILE}

# 统计每个样品的SNP
echo "统计每个样品的SNP数量（这可能需要一些时间）..."
bcftools query -f '[%SAMPLE\t%GT\n]' ${VCF_FILE} | \
awk '{
    sample=$1
    gt=$2
    
    # 跳过缺失基因型
    if (gt != "./." && gt != ".|.") {
        total[sample]++
        
        # 统计杂合和纯合
        if (gt == "0/1" || gt == "0|1" || gt == "1/0" || gt == "1|0") {
            het[sample]++
        } else if (gt == "1/1" || gt == "1|1") {
            hom_alt[sample]++
        }
    } else {
        missing[sample]++
    }
}
END {
    for (s in total) {
        print s, total[s], het[s]+0, hom_alt[s]+0, missing[s]+0
    }
}' | sort -k1,1 >> ${OUT_FILE}

echo "生成统计摘要..."

# 生成统计摘要
awk 'NR>1 {
    sum+=$2; het_sum+=$3; hom_sum+=$4; miss_sum+=$5; n++
    if (NR==2 || $2>max_snp) {max_snp=$2; max_sample=$1}
    if (NR==2 || $2<min_snp) {min_snp=$2; min_sample=$1}
}
END {
    print "=== SNP统计摘要 ==="
    print ""
    print "总样品数: " n
    print ""
    print "平均有效SNP数: " sprintf("%.0f", sum/n)
    print "平均杂合SNP数: " sprintf("%.0f", het_sum/n)
    print "平均纯合变异SNP数: " sprintf("%.0f", hom_sum/n)
    print "平均缺失位点数: " sprintf("%.0f", miss_sum/n)
    print ""
    print "最高SNP数: " max_snp " (" max_sample ")"
    print "最低SNP数: " min_snp " (" min_sample ")"
}' ${OUT_FILE} > ${SUMMARY_FILE}

echo ""
echo "========================================"
echo "统计完成！"
echo "========================================"
echo "结果保存在: ${OUT_FILE}"
echo "统计摘要保存在: ${SUMMARY_FILE}"
echo ""
echo "统计摘要预览:"
cat ${SUMMARY_FILE}
