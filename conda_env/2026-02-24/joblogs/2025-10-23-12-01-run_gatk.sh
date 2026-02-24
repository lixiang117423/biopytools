#!/bin/bash

# 这是一个完整的、遵循GATK最佳实践的变异检测流程脚本。
# 它包含了两个核心阶段：1. 为每个样本单独生成GVCF；2. 对所有GVCF进行联合基因分型。
# 脚本已为四倍体生物和旧版GATK进行了适配。

# --- 安全设置 ---
set -e
set -o pipefail

# --- 用户配置区 ---

# 1. 输入文件路径
BAM_DIR="/share/org/YZWL/yzwl_lixg/project/16.荠菜/19.筛选唯一比对上的reads/bam"
REF_GENOME="/share/org/YZWL/yzwl_lixg/project/16.荠菜/01.data/genome/genome.fa"

# 2. 输出文件夹路径
OUTPUT_DIR="/share/org/YZWL/yzwl_lixg/project/16.荠菜/19.筛选唯一比对上的reads/vcf"

# 3. 软件和资源配置
MEM_GB=550
THREADS=80

# 4. GATK Filtering Parameters
QD_FILTER="QD < 2.0"
FS_FILTER="FS > 60.0"
MQ_FILTER="MQ < 40.0"
SOR_FILTER="SOR > 3.0"
MQRankSum_FILTER="MQRankSum < -12.5"
ReadPosRankSum_FILTER="ReadPosRankSum < -8.0"

# --- 脚本主体 ---

echo "================================================="
echo "=== 完整GATK流程启动 (BAM -> VCF) ==="
echo "================================================="
echo "输入 BAM 文件夹: $BAM_DIR"
echo "输出主文件夹: $OUTPUT_DIR"
echo "参考基因组: $REF_GENOME"
echo "================================================="

# --- 准备工作: 创建目录和索引 ---
echo -e "\n--- 准备工作: 创建目录结构和准备参考基因组 ---"
GVCF_DIR="$OUTPUT_DIR/gvcfs"
DB_DIR="$OUTPUT_DIR/genomics_db"
mkdir -p $GVCF_DIR $DB_DIR

if [ ! -f "${REF_GENOME}.fai" ]; then samtools faidx $REF_GENOME; fi
REF_DICT=$(echo $REF_GENOME | sed 's/\.fasta$/.dict/' | sed 's/\.fa$/.dict/')
if [ ! -f "$REF_DICT" ]; then gatk CreateSequenceDictionary -R $REF_GENOME -O $REF_DICT; fi
echo "准备工作完成。"

# ===================================================================
# === 阶段一: 为每个样本单独 Calling 变异 (HaplotypeCaller) ===
# ===================================================================
echo -e "\n--- 阶段一: 开始为每个样本生成GVCF文件 (Ploidy=4) ---"

# 创建一个样本映射文件，供阶段二使用
GVCF_MAP_FILE="$OUTPUT_DIR/sample_map.txt" 
> $GVCF_MAP_FILE 

for input_bam in "$BAM_DIR"/*.bam
do
    SAMPLE_NAME=$(basename "$input_bam" .bam)
    echo -e "\n>>> 正在处理样本: $SAMPLE_NAME <<<"

    gvcf="$GVCF_DIR/${SAMPLE_NAME}.g.vcf.gz"
    echo -e "${SAMPLE_NAME}\t${gvcf}" >> $GVCF_MAP_FILE
    
    if [ -f "$gvcf" ]; then
        echo "样本 $SAMPLE_NAME 的GVCF文件已存在，跳过 HaplotypeCaller。"
        continue
    fi

    if [ ! -f "${input_bam}.bai" ]; then
        echo "为 $input_bam 创建索引..."
        samtools index "$input_bam"
    fi

    # 运行 HaplotypeCaller 生成GVCF
    # 【【【核心参数: --ploidy 4】】】
    gatk --java-options "-Xmx${MEM_GB}G" HaplotypeCaller \
        -R $REF_GENOME \
        -I "$input_bam" \
        -O "$gvcf" \
        -ERC GVCF \
        --ploidy 2

    echo "样本 $SAMPLE_NAME 的GVCF生成完毕。"
done
echo -e "\n--- 阶段一完成：所有样本的GVCF文件已生成 ---"


# ===================================================================
# === 阶段二: 联合基因分型 (Joint-Calling) ===
# ===================================================================
echo -e "\n--- 阶段二: 开始进行联合基因分型 ---"

# 步骤 2.1: 合并/导入 GVCFs
echo -e "\n--- 步骤 2.1: 使用GenomicsDBImport合并所有GVCF ---"
# 【版本兼容性修正】: 手动删除已存在的数据库
echo "INFO: 正在手动删除旧的GenomicsDB工作区以确保干净运行..."
rm -rf "$DB_DIR"
gatk --java-options "-Xmx${MEM_GB}G" GenomicsDBImport \
    --genomicsdb-workspace-path $DB_DIR \
    --batch-size 50 \
    --sample-name-map $GVCF_MAP_FILE \
    --reader-threads $THREADS

# 步骤 2.2: 联合基因分型
echo -e "\n--- 步骤 2.2: 使用GenotypeGVCFs进行联合基因分型 (Ploidy=4) ---"
RAW_VCF="$OUTPUT_DIR/all_samples.raw.vcf.gz"
gatk --java-options "-Xmx${MEM_GB}G" GenotypeGVCFs \
    -R $REF_GENOME \
    -V "gendb://$DB_DIR" \
    -O $RAW_VCF \
    --ploidy 2

# 步骤 2.3: 变异过滤
echo -e "\n--- 步骤 2.3: 对SNP进行硬过滤 ---"
FILTERED_VCF="$OUTPUT_DIR/all_samples.filtered.vcf.gz"
gatk VariantFiltration \
    -V $RAW_VCF \
    -O $FILTERED_VCF \
    --filter-expression "$QD_FILTER" --filter-name "QD_filter" \
    --filter-expression "$FS_FILTER" --filter-name "FS_filter" \
    --filter-expression "$MQ_FILTER" --filter-name "MQ_filter" \
    --filter-expression "$SOR_FILTER" --filter-name "SOR_filter" \
    --filter-expression "$MQRankSum_FILTER" --filter-name "MQRankSum_filter" \
    --filter-expression "$ReadPosRankSum_FILTER" --filter-name "ReadPosRankSum_filter"

echo "SNP过滤完成。"
echo "最终分析就绪的VCF文件是: $FILTERED_VCF"
echo "================================================="
echo "=== 所有分析流程成功结束! ==="
echo "================================================="