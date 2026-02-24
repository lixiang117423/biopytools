#!/bin/bash
# ------------------------------------------------------------------
# 脚本名称: run_variant_calling_pipeline.sh
# 描述: 蔊菜重测序流程：GATK Joint Call -> 过滤 -> 基因型提取 -> 统计
# ------------------------------------------------------------------

# ================= 配置区域 (Config) =================
# 1. 基础设置
set -e  # 遇到错误立即退出
START_TIME=$(date +%s)

# ---【修改这里】项目名称/物种前缀 ---
SPECIES="hancai"  # <--- 这里改成你当前项目的物种名，例如 "arabidopsis" 或 "rice"
# -----------------------------------

# 2. 资源设置
THREADS=64          # 常规步骤线程数
BCF_THREADS=64      # bcftools并行压缩线程数 (根据服务器负载调整)

# 3. 文件与路径设置
PROJECT_DIR=$(pwd)  # 当前目录作为项目根目录
REF_GENOME="01.data/genome/genome.fa"

# 输入目录
GVCF_DIR="02.mapping/02.gvcf"
BAM_DIR="02.mapping/01.bam"

# 输出目录
OUT_MAPPING="02.mapping"
OUT_JOINT="03.gatk_joint"
OUT_FILTER="04.filtered_snp_indel"
OUT_GENOTYPE="05.snp_indel_genotype"
OUT_STATS="${OUT_MAPPING}/04.bam_stat"

# 脚本路径
BAM_STAT_SCRIPT="${HOME}/software/scripts/23.bam_stats_reporter.py"

echo ">>> 流程开始运行..."
echo ">>> 项目名称(前缀): ${SPECIES}"
echo ">>> 工作目录: ${PROJECT_DIR}"

# ================= 1. 创建目录结构 =================
echo "[1/6] 创建必要的输出目录..."
mkdir -p ${OUT_STATS}
mkdir -p ${OUT_JOINT}
mkdir -p ${OUT_FILTER}
mkdir -p ${OUT_GENOTYPE}/snp
mkdir -p ${OUT_GENOTYPE}/indel

# ================= 2. Mapping (可选/已注释) =================
# 如果需要运行 Mapping，取消下方注释
# echo "[Mapping] 开始运行 BWA + GATK HaplotypeCaller..."
# bash ${BWA_SCRIPT} \
#     -i 01.data/clean \
#     -r ${REF_GENOME} \
#     -o ${OUT_MAPPING}

# ================= 3. GATK Joint Calling =================
echo "[2/6] 开始 GATK 联合变异检测 (Joint Calling)..."
biopytools gatk-joint \
    -i ${GVCF_DIR} \
    -r ${REF_GENOME} \
    -o ${OUT_JOINT}

# ================= 4. 过滤 SNP 和 INDEL =================
echo "[3/6] 过滤 SNP 和 INDEL..."
# 这里的输入文件名根据 gatk-joint 的实际输出可能需要微调
biopytools filter-snp-indel \
    -i ${OUT_JOINT}/joint_genotyping_merged.vcf.gz \
    -o ${OUT_FILTER} \
    -t ${THREADS}

# ================= 5. MAF 0.05 筛选 =================
echo "[4/6] 进行 MAF > 0.05 筛选..."
INPUT_SNP_VCF="${OUT_FILTER}/variation.filtered.snp.vcf.gz"
OUTPUT_MAF_VCF="${OUT_FILTER}/variation.filtered.snp.maf_05.vcf.gz"

if [ -f "${INPUT_SNP_VCF}" ]; then
    bcftools view \
        -i 'MAF > 0.05' \
        --threads ${BCF_THREADS} \
        -Oz \
        -o ${OUTPUT_MAF_VCF} \
        ${INPUT_SNP_VCF}
    
    bcftools index -t ${OUTPUT_MAF_VCF}
else
    echo "错误: 未找到 SNP 文件 ${INPUT_SNP_VCF}，跳过 MAF 筛选。"
    exit 1
fi

# ================= 6. 提取基因型 (CSV格式) =================
echo "[5/6] 提取基因型矩阵..."

# 6.1 提取 SNP
echo "  -> 正在提取 SNP 基因型 (前缀: ${SPECIES})..."
biopytools vcf-genotype \
    -i ${OUTPUT_MAF_VCF} \
    -o ${SPECIES} \
    --output-dir ${OUT_GENOTYPE}/snp \
    --each yes \
    -t csv

# 6.2 提取 INDEL
echo "  -> 正在提取 INDEL 基因型 (前缀: ${SPECIES})..."
INDEL_VCF="${OUT_FILTER}/variation.filtered.indel.vcf.gz"

if [ -f "${INDEL_VCF}" ]; then
    biopytools vcf-genotype \
        -i ${INDEL_VCF} \
        -o ${SPECIES} \
        --output-dir ${OUT_GENOTYPE}/indel \
        --each yes \
        -t csv
else
    echo "警告: 未找到 INDEL 文件 ${INDEL_VCF}，跳过 INDEL 提取。"
fi

# ================= 7. 比对结果统计 =================
echo "[6/6] 统计 BAM 比对结果..."
if [ -f "${BAM_STAT_SCRIPT}" ]; then
    python3 ${BAM_STAT_SCRIPT} \
        -i ${BAM_DIR} \
        -o ${OUT_STATS}/bam_stat.xlsx
else
    echo "警告: 找不到统计脚本，跳过。"
fi

# ================= 结束 =================
END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))
echo ">>> 所有任务完成！耗时: ${DURATION} 秒。"