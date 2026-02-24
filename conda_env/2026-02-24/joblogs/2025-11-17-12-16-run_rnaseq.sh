#!/bin/bash

#####################################################################
# 三代全长转录组表达量分析流程
# 使用minimap2比对 + 统计reads数量计算表达量
#####################################################################

# 设置线程数
THREADS=88

# 设置路径
GENOME="/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/66.三代转录组比对到组装的基因组/rnaseq/Orychophragmus_violaceus_OV53_1_HiFi.fa"
SAMPLE1="/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/66.三代转录组比对到组装的基因组/rnaseq/genA.fastq"
SAMPLE2="/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/66.三代转录组比对到组装的基因组/rnaseq/jingyeA.fastq"
GTF="/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/66.三代转录组比对到组装的基因组/rnaseq/complete.genomic.gtf"

# 创建输出目录
OUTDIR="pacbio_expression_analysis"
mkdir -p ${OUTDIR}/{alignment,counts,logs}

echo "=========================================="
echo "三代全长转录组表达量分析流程"
echo "开始时间: $(date)"
echo "=========================================="

#####################################################################
# 步骤1: 使用minimap2进行转录组比对
#####################################################################
echo "[$(date)] 步骤1: 开始比对转录组数据到基因组..."

# 比对genA样本
echo "[$(date)] 正在比对 genA 样本..."
minimap2 -ax splice:hq -uf \
    -t ${THREADS} \
    ${GENOME} \
    ${SAMPLE1} 2> ${OUTDIR}/logs/genA_minimap2.log | \
    samtools view -bS -@ 20 - | \
    samtools sort -@ 20 -o ${OUTDIR}/alignment/genA.sorted.bam -
samtools index -@ ${THREADS} ${OUTDIR}/alignment/genA.sorted.bam

# 比对jingyeA样本
echo "[$(date)] 正在比对 jingyeA 样本..."
minimap2 -ax splice:hq -uf \
    -t ${THREADS} \
    ${GENOME} \
    ${SAMPLE2} 2> ${OUTDIR}/logs/jingyeA_minimap2.log | \
    samtools view -bS -@ 20 - | \
    samtools sort -@ 20 -o ${OUTDIR}/alignment/jingyeA.sorted.bam -
samtools index -@ ${THREADS} ${OUTDIR}/alignment/jingyeA.sorted.bam

echo "[$(date)] 比对完成!"

#####################################################################
# 步骤2: 生成比对统计信息
#####################################################################
echo "[$(date)] 步骤2: 生成比对统计信息..."

samtools flagstat ${OUTDIR}/alignment/genA.sorted.bam > ${OUTDIR}/logs/genA_flagstat.txt
samtools flagstat ${OUTDIR}/alignment/jingyeA.sorted.bam > ${OUTDIR}/logs/jingyeA_flagstat.txt

echo "[$(date)] 比对统计完成!"

#####################################################################
# 步骤3: 使用featureCounts计算基因表达量
#####################################################################
echo "[$(date)] 步骤3: 使用featureCounts计算基因表达量..."

featureCounts -L -T ${THREADS} \
    -a ${GTF} \
    -o ${OUTDIR}/counts/gene_counts.txt \
    ${OUTDIR}/alignment/genA.sorted.bam \
    ${OUTDIR}/alignment/jingyeA.sorted.bam

echo "[$(date)] featureCounts完成!"

#####################################################################
# 步骤4: 使用StringTie进行转录本定量(基于已有GTF注释)
#####################################################################
echo "[$(date)] 步骤4: 使用StringTie基于GTF注释进行定量..."

# genA样本 - 使用-e参数仅定量已有转录本
echo "[$(date)] StringTie定量 genA..."
stringtie -e -B \
    -p ${THREADS} \
    -G ${GTF} \
    -o ${OUTDIR}/counts/genA_stringtie.gtf \
    -A ${OUTDIR}/counts/genA_gene_abundance.txt \
    ${OUTDIR}/alignment/genA.sorted.bam

# jingyeA样本 - 使用-e参数仅定量已有转录本
echo "[$(date)] StringTie定量 jingyeA..."
stringtie -e -B \
    -p ${THREADS} \
    -G ${GTF} \
    -o ${OUTDIR}/counts/jingyeA_stringtie.gtf \
    -A ${OUTDIR}/counts/jingyeA_gene_abundance.txt \
    ${OUTDIR}/alignment/jingyeA.sorted.bam

echo "[$(date)] StringTie定量完成!"

#####################################################################
# 步骤5: 整理表达量结果
#####################################################################
echo "[$(date)] 步骤5: 整理表达量矩阵..."

# 创建Python脚本整理表达量矩阵
cat > ${OUTDIR}/merge_expression.py << 'EOF'
import pandas as pd
import sys

# 读取两个样本的表达量
genA = pd.read_csv('pacbio_expression_analysis/counts/genA_gene_abundance.txt', 
                   sep='\t', skiprows=0)
jingyeA = pd.read_csv('pacbio_expression_analysis/counts/jingyeA_gene_abundance.txt', 
                      sep='\t', skiprows=0)

# 提取基因ID和TPM/FPKM值
result = pd.DataFrame({
    'Gene_ID': genA['Gene ID'],
    'Gene_Name': genA['Gene Name'],
    'genA_TPM': genA['TPM'],
    'genA_FPKM': genA['FPKM'],
    'genA_Coverage': genA['Coverage'],
    'jingyeA_TPM': jingyeA['TPM'],
    'jingyeA_FPKM': jingyeA['FPKM'],
    'jingyeA_Coverage': jingyeA['Coverage']
})

# 保存结果
result.to_csv('pacbio_expression_analysis/counts/expression_matrix.txt', 
              sep='\t', index=False)

print(f"表达量矩阵已保存到: pacbio_expression_analysis/counts/expression_matrix.txt")
print(f"共有 {len(result)} 个基因")
print("\n前5行预览:")
print(result.head())

# 额外保存一个简化版本(仅TPM)
result_simple = pd.DataFrame({
    'Gene_ID': genA['Gene ID'],
    'Gene_Name': genA['Gene Name'],
    'genA_TPM': genA['TPM'],
    'jingyeA_TPM': jingyeA['TPM']
})
result_simple.to_csv('pacbio_expression_analysis/counts/expression_matrix_TPM.txt', 
                     sep='\t', index=False)
print(f"\nTPM表达量矩阵已保存到: pacbio_expression_analysis/counts/expression_matrix_TPM.txt")
EOF

python3 ${OUTDIR}/merge_expression.py

#####################################################################
# 完成
#####################################################################
echo "=========================================="
echo "分析完成!"
echo "结束时间: $(date)"
echo "=========================================="
echo ""
echo "输出文件说明:"
echo "1. 比对文件: ${OUTDIR}/alignment/*.sorted.bam"
echo "2. 比对统计: ${OUTDIR}/logs/*_flagstat.txt"
echo "3. featureCounts结果: ${OUTDIR}/counts/gene_counts.txt"
echo "4. StringTie转录本注释: ${OUTDIR}/counts/*_stringtie.gtf"
echo "5. StringTie基因表达量: ${OUTDIR}/counts/*_gene_abundance.txt"
echo "6. 表达量矩阵(完整): ${OUTDIR}/counts/expression_matrix.txt"
echo "7. 表达量矩阵(仅TPM): ${OUTDIR}/counts/expression_matrix_TPM.txt"
echo ""
echo "查看比对率:"
echo "cat ${OUTDIR}/logs/genA_flagstat.txt"
echo "cat ${OUTDIR}/logs/jingyeA_flagstat.txt"
echo ""
echo "查看表达量矩阵:"
echo "head -n 20 ${OUTDIR}/counts/expression_matrix.txt"