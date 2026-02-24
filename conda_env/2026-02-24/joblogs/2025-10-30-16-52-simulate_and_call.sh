#!/bin/bash

# --- 用户需配置的变量 ---
REF_GENOME="120.sorted.genome.fa"    # 参考基因组
QUERY_GENOME="119.sorted.genome.fa"  # 查询基因组 (从中模拟 reads)
PREFIX="119_sim_on_120"              # 输出文件前缀
OUTPUT_BAM="${PREFIX}.sorted.bam"    # 最终输出的BAM文件名
THREADS=8                            # 使用的CPU线程数
COVERAGE=50                          # 模拟的测序深度 (X)，50X对于可视化已足够
READ_LEN=150                         # 模拟的读长 (bp)
INSERT_SIZE=500                      # 模拟的插入片段大小 (bp)

# --- 脚本主体 ---

echo "==== 步骤 1: 模拟完美的Paired-End测序数据 ===="
# -N: 计算需要生成的read-pairs总数以达到目标覆盖度
GENOME_SIZE=$(grep -v ">" ${QUERY_GENOME} | wc -c)
NUM_READS=$(awk -v size=${GENOME_SIZE} -v cov=${COVERAGE} -v len=${READ_LEN} 'BEGIN{printf "%.0f", size * cov / (2 * len)}')
echo "基因组大小: ${GENOME_SIZE} bp, 目标覆盖度: ${COVERAGE}X, 模拟Read-Pair数量: ${NUM_READS}"

# -e 0: 测序错误率为0 (关键！因为我们比较的是完美的基因组)
# -r 0: 突变率为0
# -R 0: indels比率为0
wgsim -e 0 -r 0 -R 0 -N ${NUM_READS} -1 ${READ_LEN} -2 ${READ_LEN} -d ${INSERT_SIZE} \
  ${QUERY_GENOME} ${PREFIX}.read1.fq.gz ${PREFIX}.read2.fq.gz

if [ $? -ne 0 ]; then echo "错误: wgsim 模拟数据失败。"; exit 1; fi

echo "==== 步骤 2: 为参考基因组建立 BWA 索引 ===="
if [ ! -f "${REF_GENOME}.bwt" ]; then
    echo "未找到索引，正在创建..."
    bwa index ${REF_GENOME}
else
    echo "索引已存在，跳过此步。"
fi

echo "==== 步骤 3: 使用 bwa mem 将模拟数据比对到参考基因组 ===="
# -t: 线程数
# -R: 添加Read Group信息，这是BAM文件的标准组成部分
bwa mem -t ${THREADS} -R "@RG\tID:sim\tSM:sample119\tPL:ILLUMINA" ${REF_GENOME} \
  ${PREFIX}.read1.fq.gz ${PREFIX}.read2.fq.gz | \
  samtools view -@ ${THREADS} -bS - | \
  samtools sort -@ ${THREADS} -o ${OUTPUT_BAM} -

if [ $? -ne 0 ]; then echo "错误: BWA比对或samtools排序失败。"; exit 1; fi

echo "==== 步骤 4: 为最终的 BAM 文件创建索引 ===="
samtools index -@ ${THREADS} ${OUTPUT_BAM}

if [ $? -ne 0 ]; then echo "错误: samtools index 失败。"; exit 1; fi

echo "==== 清理中间文件 ===="
# rm ${PREFIX}.read1.fq.gz ${PREFIX}.read2.fq.gz

echo "==== 完成！ ===="
echo "已成功生成用于可视化的比对文件: ${OUTPUT_BAM}"
echo "以及其索引文件: ${OUTPUT_BAM}.bai"
echo "现在您可以将这两个文件加载到 JBrowse 中查看。"