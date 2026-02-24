#!/bin/bash

# This script performs joint-calling on a set of GVCF files using the GATK best practices workflow.
# It is designed for a tetraploid organism.

# --- 安全设置 ---
# set -e: If any command fails, the script will exit immediately.
# set -o pipefail: If any command in a pipeline fails, the entire pipeline is considered failed.
set -e
set -o pipefail

# --- 用户配置区 ---

# 1. 输入文件路径
# This is the directory where you have stored all the individual sample GVCFs (e.g., sample1.g.vcf.gz)
GVCF_DIR="/share/org/YZWL/yzwl_lixg/project/16.荠菜/19.筛选唯一比对上的reads/vcf/gvcfs"
REF_GENOME="/share/org/YZWL/yzwl_lixg/project/16.荠菜/01.data/genome/genome.fa"

# 2. 输出文件夹路径
# All final outputs will be placed here
OUTPUT_DIR="/share/org/YZWL/yzwl_lixg/project/16.荠菜/19.筛选唯一比对上的reads/vcf"

# 3. 软件和资源配置
MEM_GB=550  # Max memory in GB for GATK Java options
THREADS=80  # Number of threads for parallel processes

# 4. GATK Filtering Parameters (Standard hard-filtering recommendations)
QD_FILTER="QD < 2.0"
FS_FILTER="FS > 60.0"
MQ_FILTER="MQ < 40.0"
SOR_FILTER="SOR > 3.0"
MQRankSum_FILTER="MQRankSum < -12.5"
ReadPosRankSum_FILTER="ReadPosRankSum < -8.0"

# --- 脚本主体 ---

echo "================================================="
echo "=== GATK Joint-Calling Workflow for Tetraploids ==="
echo "================================================="
echo "Input GVCF Directory: $GVCF_DIR"
echo "Output Directory:     $OUTPUT_DIR"
echo "Reference Genome:     $REF_GENOME"
echo "================================================="

# --- 步骤 0: 准备工作 ---
echo -e "\n--- Step 0: Preparing environment and creating directories ---"
DB_DIR="$OUTPUT_DIR/genomics_db" # Directory for the GenomicsDB
mkdir -p $DB_DIR

# Check if the reference genome and its dictionaries exist
if [ ! -f "${REF_GENOME}.fai" ]; then
    echo "Reference .fai index not found. Creating..."
    samtools faidx $REF_GENOME
fi
REF_DICT=$(echo $REF_GENOME | sed 's/\.fasta$/.dict/' | sed 's/\.fa$/.dict/')
if [ ! -f "$REF_DICT" ]; then
    echo "Reference .dict dictionary not found. Creating..."
    gatk CreateSequenceDictionary -R $REF_GENOME -O $REF_DICT
fi
echo "Preparation complete."

# --- 步骤 1: 创建样本映射文件 ---
echo -e "\n--- Step 1: Creating a sample map file for GenomicsDBImport ---"
GVCF_MAP_FILE="$OUTPUT_DIR/sample_map.txt"
> "$GVCF_MAP_FILE" # Create or clear the file

for gvcf in "$GVCF_DIR"/*.g.vcf.gz; do
    if [ -f "$gvcf" ]; then
        SAMPLE_NAME=$(basename "$gvcf" .g.vcf.gz)
        echo -e "${SAMPLE_NAME}\t${gvcf}" >> "$GVCF_MAP_FILE"
    fi
done
echo "Sample map file created at: $GVCF_MAP_FILE"

# --- 步骤 2: 合并/导入 GVCFs ---
echo -e "\n--- Step 2: Consolidating GVCFs using GenomicsDBImport ---"
# This tool creates a database from all GVCFs, which is much more efficient for the next step.
gatk --java-options "-Xmx${MEM_GB}G" GenomicsDBImport \
    --genomicsdb-workspace-path $DB_DIR \
    --batch-size 50 \
    --sample-name-map $GVCF_MAP_FILE \
    --reader-threads $THREADS \
    --overwrite-existing-genomicsdb

# --- 步骤 3: 联合基因分型 ---
echo -e "\n--- Step 3: Performing joint-genotyping with GenotypeGVCFs ---"
RAW_VCF="$OUTPUT_DIR/all_samples.raw.vcf.gz"
# This step reads the database and calls variants across all samples simultaneously.
# 【【【 CRITICAL PARAMETER: --ploidy 4 】】】
gatk --java-options "-Xmx${MEM_GB}G" GenotypeGVCFs \
    -R $REF_GENOME \
    -V "gendb://$DB_DIR" \
    -O $RAW_VCF \
    --ploidy 4

# --- 步骤 4: 变异过滤 ---
echo -e "\n--- Step 4: Applying hard filters to the raw VCF ---"
FILTERED_VCF="$OUTPUT_DIR/all_samples.filtered.vcf.gz"
# This step flags low-quality variants based on the annotation values.
gatk VariantFiltration \
    -V $RAW_VCF \
    -O $FILTERED_VCF \
    --filter-expression "$QD_FILTER" --filter-name "QD_filter" \
    --filter-expression "$FS_FILTER" --filter-name "FS_filter" \
    --filter-expression "$MQ_FILTER" --filter-name "MQ_filter" \
    --filter-expression "$SOR_FILTER" --filter-name "SOR_filter" \
    --filter-expression "$MQRankSum_FILTER" --filter-name "MQRankSum_filter" \
    --filter-expression "$ReadPosRankSum_FILTER" --filter-name "ReadPosRankSum_filter"

echo "SNP filtering complete."
echo "The final, analysis-ready VCF file is: $FILTERED_VCF"
echo "================================================="
echo "=== All analysis pipelines completed successfully! ==="
echo "================================================="
