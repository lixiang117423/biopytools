#!/bin/bash

# 设置输入输出路径
VCF_DIR="/share/org/YZWL/yzwl_lixg/project/21.野生大豆群体/05.所有样品/03.按染色体过滤VCF文件"
OUT_DIR="/share/org/YZWL/yzwl_lixg/project/21.野生大豆群体/05.所有样品/15.所有样品SNP数量统计"

# 创建输出文件夹
mkdir -p ${OUT_DIR}

# 输出文件
OUT_FILE="${OUT_DIR}/sample_snp_counts.txt"
SUMMARY_FILE="${OUT_DIR}/summary_statistics.txt"

echo "开始统计每个样品的有效SNP数量..."

# 初始化输出文件
echo -e "Sample\tTotal_SNPs\tHet_SNPs\tHom_Alt_SNPs\tMissing_SNPs" > ${OUT_FILE}

# 遍历所有染色体的VCF文件
for chr in {1..20}; do
    VCF_FILE="${VCF_DIR}/chr${chr}.snp.vcf.gz"
    
    if [ ! -f "${VCF_FILE}" ]; then
        echo "警告: ${VCF_FILE} 不存在，跳过"
        continue
    fi
    
    echo "处理 chr${chr}..."
    
    # 提取样品名称（仅在第一个染色体时执行）
    if [ ${chr} -eq 1 ]; then
        bcftools query -l ${VCF_FILE} > ${OUT_DIR}/sample_list.txt
    fi
    
    # 统计每个样品的SNP
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
    }' > ${OUT_DIR}/chr${chr}_counts.tmp
done

echo "合并所有染色体的统计结果..."

# 合并所有染色体的统计结果
awk '{
    sample=$1
    total[sample]+=$2
    het[sample]+=$3
    hom_alt[sample]+=$4
    missing[sample]+=$5
}
END {
    for (s in total) {
        print s, total[s], het[s], hom_alt[s], missing[s]
    }
}' ${OUT_DIR}/chr*_counts.tmp | sort -k1,1 >> ${OUT_FILE}

# 生成统计摘要
echo "生成统计摘要..."
awk 'NR>1 {
    sum+=$2; het_sum+=$3; hom_sum+=$4; miss_sum+=$5; n++
}
END {
    print "总样品数: " n
    print "平均SNP数: " sum/n
    print "平均杂合SNP数: " het_sum/n
    print "平均纯合SNP数: " hom_sum/n
    print "平均缺失SNP数: " miss_sum/n
    print "\n最高SNP数: " max_snp
    print "最低SNP数: " min_snp
}' ${OUT_FILE} > ${SUMMARY_FILE}

# 清理临时文件
rm -f ${OUT_DIR}/chr*_counts.tmp

echo "统计完成！"
echo "结果保存在: ${OUT_FILE}"
echo "统计摘要保存在: ${SUMMARY_FILE}"
