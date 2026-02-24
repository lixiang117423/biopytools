#!/bin/bash

# --- 用户需配置的变量 ---
REF_GENOME="/share/org/YZWL/yzwl_lixg/project/19.大豆疫霉菌/01.data/genome/Phytophthora_sojae_JS2.fa"
GVCF_DIR="/share/org/YZWL/yzwl_lixg/project/19.大豆疫霉菌/02.mapping/vcf"
OUTPUT_DIR="/share/org/YZWL/yzwl_lixg/project/19.大豆疫霉菌/03.glnexus_joint"
# !!! 关键：根据您的gVCF来源修改此配置 !!!
GLNEXUS_CONFIG="gatk" 


# --- 脚本开始 ---
echo "步骤 0: 创建输出目录"
mkdir -p ${OUTPUT_DIR}

# 定义输出文件路径
BED_FILE="${OUTPUT_DIR}/genome.bed"
BCF_OUT="${OUTPUT_DIR}/merged_output.bcf"
VCF_OUT="${OUTPUT_DIR}/merged_output.vcf.gz"

echo "步骤 1: 准备参考基因组 BED 文件"
# 检查是否存在fai索引，没有则创建
if [ ! -f "${REF_GENOME}.fai" ]; then
    echo "正在为参考基因组创建 .fai 索引..."
    samtools faidx ${REF_GENOME}
fi
# 从fai文件生成BED文件
awk -v OFS='\t' '{print $1, 0, $2}' ${REF_GENOME}.fai > ${BED_FILE}
echo "BED 文件已创建: ${BED_FILE}"

echo "步骤 2: 运行 glnexus_cli"
glnexus_cli \
  --config ${GLNEXUS_CONFIG} \
  --bed ${BED_FILE} \
  ${GVCF_DIR}/*.g.vcf.gz \
  --threads 88 \
  > ${BCF_OUT}
  
echo "GLnexus 运行完成，BCF 输出至: ${BCF_OUT}"

echo "步骤 3: 后处理 - 转换为 VCF 并索引"
bcftools view ${BCF_OUT} | bgzip -c > ${VCF_OUT}
bcftools index ${VCF_OUT}
echo "已成功转换为 VCF.gz 并创建索引: ${VCF_OUT}"

echo "所有步骤完成！"