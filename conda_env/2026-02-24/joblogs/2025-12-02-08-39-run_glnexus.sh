#!/bin/bash

# --- 脚本安全设置 ---
# 遇到错误立即停止脚本，防止后续步骤对空文件报错
set -e 

# --- 用户需配置的变量 ---
REF_GENOME="/share/org/YZWL/yzwl_lixg/project/21.野生大豆群体/01.data/genome/genome.fa"
GVCF_DIR="/share/org/YZWL/yzwl_lixg/project/21.野生大豆群体/03.glnexus_joint/vcf"
OUTPUT_DIR="/share/org/YZWL/yzwl_lixg/project/21.野生大豆群体/03.glnexus_joint"
GLNEXUS_CONFIG="gatk"

# 定义 GLnexus 数据库的临时目录路径 (建议放在输出目录下，方便管理)
DB_DIR="${OUTPUT_DIR}/GLnexus.DB"

# --- 脚本开始 ---
echo "步骤 0: 检查及创建输出目录"
mkdir -p ${OUTPUT_DIR}

# 方案一：使用for循环为 VCF 建索引
# (注意：如果文件非常多，这里可能会跑很久，如果已有索引建议注释掉)
echo "正在检查/创建 VCF 索引..."
for file in ${GVCF_DIR}/*.gz; do
    if [ ! -f "${file}.tbi" ]; then
        tabix -@ 64 -p vcf "$file"  # 稍微降低索引线程数，避免卡顿
    fi
done

# 定义输出文件路径
BED_FILE="${OUTPUT_DIR}/genome.bed"
BCF_OUT="${OUTPUT_DIR}/merged_output.bcf"
VCF_OUT="${OUTPUT_DIR}/merged_output.vcf.gz"

echo "步骤 1: 准备参考基因组 BED 文件"
if [ ! -f "${REF_GENOME}.fai" ]; then
    echo "正在为参考基因组创建 .fai 索引..."
    samtools faidx ${REF_GENOME}
fi
awk -v OFS='\t' '{print $1, 0, $2}' ${REF_GENOME}.fai > ${BED_FILE}
echo "BED 文件已创建: ${BED_FILE}"

# --- [关键修改] 清理旧数据库 ---
echo "检查是否存在旧的 GLnexus 数据库..."
if [ -d "${DB_DIR}" ]; then
    echo "发现旧数据库目录 ${DB_DIR}，正在删除..."
    rm -rf "${DB_DIR}"
elif [ -d "GLnexus.DB" ]; then
    # 防止你在当前目录下运行脚本，产生的遗留文件
    echo "发现当前目录下存在 GLnexus.DB，正在删除..."
    rm -rf "GLnexus.DB"
fi

echo "步骤 2: 运行 glnexus_cli"
# 注意：--threads 88 非常高，请确保您的服务器内存足够（通常 glnexus 需要大量内存）
# 增加了 --dir 参数明确指定数据库位置
glnexus_cli \
  --config ${GLNEXUS_CONFIG} \
  --bed ${BED_FILE} \
  --dir ${DB_DIR} \
  --threads 88 \
  ${GVCF_DIR}/*.g.vcf.gz \
  > ${BCF_OUT}

echo "GLnexus 运行完成，BCF 输出至: ${BCF_OUT}"

echo "步骤 3: 后处理 - 转换为 VCF 并索引"
# 只有上面的步骤成功了，才会运行这里
bcftools view --threads 64 ${BCF_OUT} | bgzip -c -@ 64 > ${VCF_OUT}
bcftools index ${VCF_OUT}

echo "所有步骤完成！输出文件：${VCF_OUT}"