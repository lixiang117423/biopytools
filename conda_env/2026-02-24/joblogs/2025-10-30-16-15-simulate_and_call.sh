#!/bin/bash

# --- 用户需配置的变量 ---
REF_GENOME="120.sorted.genome.fa"    # 参考基因组
QUERY_GENOME="119.sorted.genome.fa"  # 查询基因组 (从中模拟 reads)
PREFIX="119_sim_on_120"              # 输出文件前缀
THREADS=80                           # 使用的CPU线程数
COVERAGE=50                          # 模拟的测序深度 (X)
READ_LEN=150                         # 模拟的读长 (bp)
INSERT_SIZE=500                      # 模拟的插入片段大小 (bp)
MIN_INDEL_LENGTH=15                  # 最终筛选的INDEL最小长度

# --- 脚本主体 ---

# --- 步骤 1: 模拟测序数据 ---
echo "步骤 1: 使用 wgsim 从 ${QUERY_GENOME} 模拟 ${COVERAGE}X 的测序数据..."
# -N: 要生成的 read-pairs 数量。计算方法: (基因组大小 * 覆盖度) / (2 * 读长)
#     我们用一个 awk 小技巧来计算基因组大小
GENOME_SIZE=$(grep -v ">" ${QUERY_GENOME} | wc -c)
NUM_READS=$(awk -v size=${GENOME_SIZE} -v cov=${COVERAGE} -v len=${READ_LEN} 'BEGIN{printf "%.0f", size * cov / (2 * len)}')

# -e 0: 测序错误率为0 (关键！因为我们比较的是完美的基因组)
# -r 0: 突变率为0 (关键！)
# -R 0: indels比率为0 (关键！)
# -1 ${READ_LEN} -2 ${READ_LEN}: Paired-end reads长度
# -d ${INSERT_SIZE}: 插入片段大小
wgsim -e 0 -r 0 -R 0 -N ${NUM_READS} -1 ${READ_LEN} -2 ${READ_LEN} -d ${INSERT_SIZE} \
  ${QUERY_GENOME} ${PREFIX}.read1.fq.gz ${PREFIX}.read2.fq.gz

# --- 步骤 2: 为参考基因组建立索引 ---
echo "步骤 2: 使用 bwa-mem2 为参考基因组 ${REF_GENOME} 建立索引..."
if [ ! -f "${REF_GENOME}.bwt" ]; then
    bwa-mem2 index ${REF_GENOME}
fi

# --- 步骤 3: 比对模拟数据到参考基因组 ---
echo "步骤 3: 使用 bwa-mem2 将模拟数据比对到参考基因组..."
# -t: 线程数
# -R: 添加 Read Group 信息，是 GATK 等工具的最佳实践
bwa-mem2 mem -t ${THREADS} -R "@RG\tID:sim\tSM:sample119" ${REF_GENOME} \
  ${PREFIX}.read1.fq.gz ${PREFIX}.read2.fq.gz | \
  samtools view -bS - | \
  samtools sort -@ ${THREADS} -o ${PREFIX}.sorted.bam -

# --- 步骤 4: 为 BAM 文件建立索引 ---
echo "步骤 4: 使用 samtools index 为 BAM 文件建立索引..."
samtools index ${PREFIX}.sorted.bam

# --- 步骤 5: 变异检测 ---
echo "步骤 5: 使用 bcftools 进行变异检测..."
# -f: 指定参考基因组
# -O v: 输出未压缩的 VCF
bcftools mpileup -f ${REF_GENOME} ${PREFIX}.sorted.bam | \
  bcftools call -mv -O v -o ${PREFIX}.raw.vcf

# --- 步骤 6: 筛选大型 INDEL ---
echo "步骤 6: 筛选长度大于 ${MIN_INDEL_LENGTH}bp 的 INDEL..."
# ILEN 是 bcftools 自动计算的 INDEL 长度 (插入为正, 缺失为负)
bcftools filter -O v -o ${PREFIX}.large_indels.vcf \
  -i "ILEN > ${MIN_INDEL_LENGTH} || ILEN < -${MIN_INDEL_LENGTH}" \
  ${PREFIX}.raw.vcf

# --- 步骤 7: 格式化输出 ---
echo "步骤 7: 将 VCF 结果转换为表格..."
(echo -e "Chr\tStart\tEnd\tType\tLength\tRefSeq\tAltSeq" && \
 bcftools query -f '%CHROM\t%POS\t%END\t%TYPE\t%INFO/ILEN\t%REF\t%ALT\n' ${PREFIX}.large_indels.vcf \
) > large_indels_gt${MIN_INDEL_LENGTH}bp_from_sim.tsv

echo "完成！"
echo "最终结果保存在: large_indels_gt${MIN_INDEL_LENGTH}bp_from_sim.tsv"
echo "标准的VCF结果保存在: ${PREFIX}.large_indels.vcf"

# 清理中间文件 (可选)
# echo "正在清理中间文件..."
# rm ${PREFIX}.read1.fq.gz ${PREFIX}.read2.fq.gz ${PREFIX}.raw.vcf
# rm ${PREFIX}.sorted.bam ${PREFIX}.sorted.bam.bai
