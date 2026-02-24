#!/bin/bash

set -e # 任何命令失败，脚本立即退出

# --- 用户需配置的变量 ---
REF_GENOME="120.sorted.genome.fa"    # 参考基因组
QUERY_GENOME="119.sorted.genome.fa"  # 查询基因组 (从中模拟 reads)
PREFIX="119_sim_on_120_gatk"         # 输出文件前缀
THREADS=80                           # 使用的CPU线程数
MEMORY="160g"                        # 分配给GATK的Java内存

# --- 模拟参数 ---
COVERAGE=50                          # 模拟的测序深度 (X)
READ_LEN=150                         # 模拟的读长 (bp)
INSERT_SIZE=500                      # 模拟的插入片段大小 (bp)

# --- 筛选参数 ---
MIN_INDEL_LENGTH=15                  # 最终筛选的INDEL最小长度

# --- 脚本主体 ---

echo "==== 步骤 1: 模拟完美的Paired-End测序数据 ===="
GENOME_SIZE=$(grep -v ">" ${QUERY_GENOME} | wc -c)
NUM_READS=$(awk -v size=${GENOME_SIZE} -v cov=${COVERAGE} -v len=${READ_LEN} 'BEGIN{printf "%.0f", size * cov / (2 * len)}')
wgsim -e 0 -r 0 -R 0 -N ${NUM_READS} -1 ${READ_LEN} -2 ${READ_LEN} -d ${INSERT_SIZE} \
  ${QUERY_GENOME} ${PREFIX}.read1.fq.gz ${PREFIX}.read2.fq.gz

echo "==== 步骤 2: 为参考基因组建立 BWA 和 Samtools 索引 ===="
if [ ! -f "${REF_GENOME}.bwt" ]; then
    bwa index ${REF_GENOME}
fi
if [ ! -f "${REF_GENOME}.fai" ]; then
    samtools faidx ${REF_GENOME}
fi
# GATK需要一个.dict文件
REF_DICT=$(echo $REF_GENOME | sed 's/\.fa$/.dict/;s/\.fasta$/.dict/')
if [ ! -f "${REF_DICT}" ]; then
    gatk CreateSequenceDictionary -R ${REF_GENOME} -O ${REF_DICT}
fi

echo "==== 步骤 3: 比对、排序并添加Read Group ===="
# -R: Read Group信息是GATK流程的必需品
RG_INFO="@RG\tID:sim_group\tSM:sample119\tPL:ILLUMINA\tLB:lib1"
bwa mem -t ${THREADS} -R "${RG_INFO}" ${REF_GENOME} \
  ${PREFIX}.read1.fq.gz ${PREFIX}.read2.fq.gz | \
  samtools view -@ ${THREADS} -bS - | \
  samtools sort -@ ${THREADS} -o ${PREFIX}.sorted.bam -

# BAM文件已生成，可用于可视化
echo "==== BAM文件已生成: ${PREFIX}.sorted.bam ===="

echo "==== 步骤 4: 标记PCR重复 (GATK最佳实践) ===="
gatk --java-options "-Xmx${MEMORY}" MarkDuplicates \
  -I ${PREFIX}.sorted.bam \
  -O ${PREFIX}.dedup.bam \
  -M ${PREFIX}.dedup_metrics.txt

echo "==== 步骤 5: 为处理后的BAM文件建立索引 ===="
samtools index ${PREFIX}.dedup.bam

echo "==== 步骤 6: 使用GATK HaplotypeCaller进行变异检测 ===="
# -ploidy 1: 假设是单倍体基因组（如细菌），如果是二倍体请使用 -ploidy 2
gatk --java-options "-Xmx${MEMORY}" HaplotypeCaller \
  -R ${REF_GENOME} \
  -I ${PREFIX}.dedup.bam \
  -O ${PREFIX}.raw.vcf.gz \
  -ploidy 2

echo "==== 步骤 7: 筛选出INDEL变异 ===="
gatk SelectVariants \
  -R ${REF_GENOME} \
  -V ${PREFIX}.raw.vcf.gz \
  --select-type-to-include INDEL \
  -O ${PREFIX}.raw_indels.vcf.gz

echo "==== 步骤 8: 对INDEL进行质量过滤和长度筛选 ===="
# GATK VariantFiltration 添加过滤标签
gatk VariantFiltration \
  -R ${REF_GENOME} \
  -V ${PREFIX}.raw_indels.vcf.gz \
  -O ${PREFIX}.filtered_indels.vcf.gz \
  --filter-expression "QD < 2.0 || FS > 200.0 || SOR > 10.0 || MQ < 40.0" \
  --filter-name "StandardFilters"

# bcftools view 筛选出通过过滤(PASS)且长度符合要求的INDEL
# 'abs(ILEN)' 计算INDEL长度的绝对值
bcftools view \
  -f PASS \
  -i "abs(ILEN) > ${MIN_INDEL_LENGTH}" \
  -O v -o ${PREFIX}.final_large_indels.vcf \
  ${PREFIX}.filtered_indels.vcf.gz

echo "==== 步骤 9: 格式化输出最终结果 ===="
(echo -e "Chr\tStart\tEnd\tType\tLength\tRefSeq\tAltSeq\tQUAL\tFilterStatus" && \
 bcftools query -f '%CHROM\t%POS\t%END\t%TYPE\t%INFO/ILEN\t%REF\t%ALT\t%QUAL\t%FILTER\n' ${PREFIX}.final_large_indels.vcf \
) > large_indels_gt${MIN_INDEL_LENGTH}bp_from_gatk.tsv

echo "==== 清理中间文件 ===="
# 如果需要，可以取消下面的注释来删除大的中间文件
# rm ${PREFIX}.read1.fq.gz ${PREFIX}.read2.fq.gz
# rm ${PREFIX}.sorted.bam ${PREFIX}.dedup.bam ${PREFIX}.dedup.bam.bai
# rm ${PREFIX}.raw.vcf.gz* ${PREFIX}.raw_indels.vcf.gz* ${PREFIX}.filtered_indels.vcf.gz*

echo "==== 所有步骤完成！ ===="
echo "主要产出文件:"
echo "1. 可视化BAM文件: ${PREFIX}.sorted.bam (和它的索引 .bai)"
echo "2. 最终筛选出的VCF: ${PREFIX}.final_large_indels.vcf"
echo "3. 最终的表格结果: large_indels_gt${MIN_INDEL_LENGTH}bp_from_gatk.tsv"