#!/bin/bash

# 这是一个为荠菜（异源四倍体）设计的，从FASTQ到VCF再到杂合度分析的完整流程脚本。
# 所有路径已根据用户提供的信息配置完毕。

# --- 安全设置 ---
# set -e: 如果任何命令失败，脚本将立即退出。
# set -o pipefail: 管道中的任何命令失败，整个管道都算作失败。
set -e
set -o pipefail

# --- 用户配置区 ---

# 1. 输入文件路径 (已配置)
FQ_DIR="/share/org/YZWL/yzwl_lixg/project/16.荠菜/01.data/clean"
REF_GENOME="/share/org/YZWL/yzwl_lixg/project/16.荠菜/01.data/genome/genome.fa"

# 2. 输出文件夹路径 (已配置)
OUTPUT_DIR="/share/org/YZWL/yzwl_lixg/project/16.荠菜/20.四倍体参数重新比对和检测变异"

# 3. 软件和资源配置
THREADS=80 # 根据您的服务器CPU核心数调整
MEM_GB=400  # GATK使用的最大内存(GB)，根据您的服务器内存调整

# 4. GATK 过滤参数 (这些是常用的硬过滤标准，可根据需要调整)
QD_FILTER="QD < 2.0"
FS_FILTER="FS > 60.0"
MQ_FILTER="MQ < 40.0"
SOR_FILTER="SOR > 3.0"
MQRankSum_FILTER="MQRankSum < -12.5"
ReadPosRankSum_FILTER="ReadPosRankSum < -8.0"

# --- 脚本主体 ---

echo "================================================="
echo "=== 四倍体荠菜 GWAS 上游分析流程开始 ==="
echo "================================================="
echo "输入 FASTQ 文件夹: $FQ_DIR"
echo "输出主文件夹: $OUTPUT_DIR"
echo "参考基因组: $REF_GENOME"
echo "================================================="

# --- 步骤 0: 准备工作 (创建目录和索引) ---
echo -e "\n--- 步骤 0: 创建目录结构和参考基因组索引 ---"
BAM_DIR="$OUTPUT_DIR/01_bam_files"
GVCF_DIR="$OUTPUT_DIR/02_gvcfs"
VCF_DIR="$OUTPUT_DIR/03_vcfs"
STATS_DIR="$OUTPUT_DIR/04_stats"
DB_DIR="$OUTPUT_DIR/db" # 用于GenomicsDBImport

mkdir -p $BAM_DIR $GVCF_DIR $VCF_DIR $STATS_DIR $DB_DIR

# 检查参考基因组是否存在
if [ ! -f "$REF_GENOME" ]; then
    echo "错误: 参考基因组文件未找到: $REF_GENOME"
    exit 1
fi

# 为参考基因组创建索引 (如果不存在)
if [ ! -f "${REF_GENOME}.bwt" ]; then
    echo "为 BWA 创建索引..."
    bwa index $REF_GENOME
fi
if [ ! -f "${REF_GENOME}.fai" ]; then
    echo "为 Samtools 创建 .fai 索引..."
    samtools faidx $REF_GENOME
fi
REF_DICT=$(echo $REF_GENOME | sed 's/\.fasta$/.dict/' | sed 's/\.fa$/.dict/')
if [ ! -f "$REF_DICT" ]; then
    echo "为 GATK 创建 .dict 字典..."
    gatk CreateSequenceDictionary -R $REF_GENOME -O $REF_DICT
fi
echo "准备工作完成。"


# --- 步骤 1, 2, 3: 循环处理每个样本 (比对, 排序, 去重, HaplotypeCaller) ---
echo -e "\n--- 步骤 1-3: 开始循环处理每个样本 ---"

# 创建一个数组来存储所有样本的GVCF文件路径
GVCF_LIST=()
GVCF_MAP_FILE="$OUTPUT_DIR/sample_map.txt" # GATK 4.2+ 推荐使用 map 文件
> $GVCF_MAP_FILE # 清空或创建 map 文件

for fq1 in $FQ_DIR/*_1.clean.fq.gz
do
    # 从文件名中提取样本名
    SAMPLE_NAME=$(basename $fq1 _1.clean.fq.gz)
    fq2=${fq1/_1.clean.fq.gz/_2.clean.fq.gz}

    echo -e "\n>>> 正在处理样本: $SAMPLE_NAME <<<"

    # 定义输出文件名
    bam_sorted="$BAM_DIR/${SAMPLE_NAME}.sorted.bam"
    bam_dedup="$BAM_DIR/${SAMPLE_NAME}.sorted.dedup.bam"
    gvcf="$GVCF_DIR/${SAMPLE_NAME}.g.vcf.gz"
    
    # 将样本名和gvcf路径写入map文件
    echo -e "${SAMPLE_NAME}\t${gvcf}" >> $GVCF_MAP_FILE
    
    # 【优化】如果最终的GVCF文件已存在，则跳过此样本，方便续跑
    if [ -f "$gvcf" ]; then
        echo "样本 $SAMPLE_NAME 的GVCF文件已存在，跳过前3步。"
        continue
    fi

    # 步骤 1: BWA-MEM 比对 + Samtools 排序
    # RG = Read Group, 包含样本信息，对下游GATK分析至关重要
    RG_INFO="@RG\\tID:${SAMPLE_NAME}\\tSM:${SAMPLE_NAME}\\tPL:ILLUMINA"
    echo "步骤 1: 正在比对和排序..."
    bwa mem -t $THREADS -R "$RG_INFO" $REF_GENOME $fq1 $fq2 | samtools sort -@ $THREADS -o $bam_sorted -

    # 步骤 2: GATK 标记PCR重复
    echo "步骤 2: 正在标记PCR重复..."
    gatk --java-options "-Xmx${MEM_GB}G" MarkDuplicates \
        -I $bam_sorted \
        -O $bam_dedup \
        -M "$BAM_DIR/${SAMPLE_NAME}.dedup_metrics.txt"
    
    # 为去重后的BAM文件创建索引
    samtools index $bam_dedup

    # 删除中间的 sorted bam 文件以节省空间
    rm $bam_sorted

    # 步骤 3: GATK HaplotypeCaller, 生成GVCF
    # 【【【核心参数: --ploidy 4】】】
    echo "步骤 3: 正在运行 HaplotypeCaller (Ploidy=4)..."
    gatk --java-options "-Xmx${MEM_GB}G" HaplotypeCaller \
        -R $REF_GENOME \
        -I $bam_dedup \
        -O $gvcf \
        -ERC GVCF \
        --ploidy 4

    echo "样本 $SAMPLE_NAME 处理完成。"
done


# --- 步骤 4: 合并GVCF并进行联合基因分型 ---
echo -e "\n--- 步骤 4: 合并所有样本并进行联合基因分型 ---"
DB_PATH="$DB_DIR/capsella_db"
RAW_VCF="$VCF_DIR/all_samples.raw.vcf.gz"

echo "步骤 4.1: 导入 GVCFs 到 GenomicsDB..."
gatk --java-options "-Xmx${MEM_GB}G" GenomicsDBImport \
    --genomicsdb-workspace-path $DB_PATH \
    --batch-size 50 \
    --sample-name-map $GVCF_MAP_FILE \
    --reader-threads $THREADS \
    --overwrite-existing-genomicsdb

echo "步骤 4.2: 运行 GenotypeGVCFs 进行联合基因分型..."
gatk --java-options "-Xmx${MEM_GB}G" GenotypeGVCFs \
    -R $REF_GENOME \
    -V "gendb://$DB_PATH" \
    -O $RAW_VCF \
    --ploidy 4

# --- 步骤 5: 变异过滤 ---
echo -e "\n--- 步骤 5: 对SNP进行硬过滤 ---"
FILTERED_VCF="$VCF_DIR/all_samples.filtered_snps.vcf.gz"
gatk VariantFiltration \
    -V $RAW_VCF \
    -O $FILTERED_VCF \
    --filter-expression "$QD_FILTER" --filter-name "QD_filter" \
    --filter-expression "$FS_FILTER" --filter-name "FS_filter" \
    --filter-expression "$MQ_FILTER" --filter-name "MQ_filter" \
    --filter-expression "$SOR_FILTER" --filter-name "SOR_filter" \
    --filter-expression "$MQRankSum_FILTER" --filter-name "MQRankSum_filter" \
    --filter-expression "$ReadPosRankSum_FILTER" --filter-name "ReadPosRankSum_filter"
echo "SNP过滤完成, 最终VCF文件为: $FILTERED_VCF"


# --- 步骤 6: 杂合度计算和Reads统计 ---
echo -e "\n--- 步骤 6: 计算杂合度并统计等位基因Reads数 ---"
FINAL_VCF_FOR_STATS="$FILTERED_VCF" # 使用过滤后的VCF进行统计

# 步骤 6.1: 计算每个样本的杂合度
echo "步骤 6.1: 计算每个样本的杂合度 (使用 vcftools)..."
vcftools --gzvcf $FINAL_VCF_FOR_STATS --het -c > "$STATS_DIR/heterozygosity_by_sample.txt"
# O(HOM)/N(sites) 是纯合子比例, (N_SITES - O(HOM))/N_SITES 是杂合率

# 步骤 6.2: 计算每个位点的杂合度 (Allele Frequency)
echo "步骤 6.2: 计算每个位点的等位基因频率 (使用 vcftools)..."
vcftools --gzvcf $FINAL_VCF_FOR_STATS --freq2 -c > "$STATS_DIR/allele_frequency_by_site.txt"
# HET = 1 - (p^2 + q^2 + ...) 可从频率计算

# 步骤 6.3: 提取杂合位点的REF和ALT reads数量 (AD字段)
echo "步骤 6.3: 提取杂合位点的等位基因深度 (使用 bcftools)..."
OUTPUT_AD_STATS="$STATS_DIR/heterozygous_allele_depths.tsv"
echo -e "CHROM\tPOS\tSAMPLE\tGENOTYPE\tREF_DEPTH\tALT_DEPTH" > $OUTPUT_AD_STATS
bcftools query -i 'GT="het"' -f '%CHROM\t%POS[\t%SAMPLE\t%GT\t%AD]\n' $FINAL_VCF_FOR_STATS | sed 's/,/\t/g' >> $OUTPUT_AD_STATS

echo "统计分析完成，结果保存在 $STATS_DIR 文件夹下。"
echo "================================================="
echo "=== 所有分析流程成功结束! ==="
echo "================================================="
