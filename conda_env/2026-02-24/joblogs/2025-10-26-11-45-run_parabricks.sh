# cp /share/org/YZWL/yzwl_lixg/project/18.小花糖芥/01.data/clean/EcAA_* /share/org/YZWL/yzwl_lixg/project/18.小花糖芥/02.小花糖芥抗感标记/01.data
# cp /share/org/YZWL/yzwl_lixg/project/18.小花糖芥/01.data/clean/KostanzA_* /share/org/YZWL/yzwl_lixg/project/18.小花糖芥/02.小花糖芥抗感标记/01.data

# biopytools parabricks \
#   -i /share/org/YZWL/yzwl_lixg/project/18.小花糖芥/02.小花糖芥抗感标记/01.data \
#   -o /share/org/YZWL/yzwl_lixg/project/18.小花糖芥/02.小花糖芥抗感标记/02.mapping \
#   -r /share/org/YZWL/yzwl_lixg/project/18.小花糖芥/02.小花糖芥抗感标记/01.data/Erysimum_cheiranthoides_chromosomes.fna \
#   -t 64

# biopytools gatk-joint \
#     -i /share/org/YZWL/yzwl_lixg/project/18.小花糖芥/02.小花糖芥抗感标记/02.mapping \
#     -o /share/org/YZWL/yzwl_lixg/project/18.小花糖芥/02.小花糖芥抗感标记/03.gatk_joint \
#     -r /share/org/YZWL/yzwl_lixg/project/18.小花糖芥/02.小花糖芥抗感标记/01.data/Erysimum_cheiranthoides_chromosomes.fna

# biopytools filter-snp-indel \
#     -i /share/org/YZWL/yzwl_lixg/project/18.小花糖芥/02.小花糖芥抗感标记/03.gatk_joint/joint_genotyping_raw.vcf.gz \
#     -o /share/org/YZWL/yzwl_lixg/project/18.小花糖芥/02.小花糖芥抗感标记/04.filtered_snp_indel \
#     -t 88

#!/bin/bash
set -e # 脚本中任何命令执行失败，则立即退出

# =================================================================
#               配置部分 (Configuration Section)
#       只需要在这里修改路径和参数，下游代码会自动更新
# =================================================================

# --- 1. 基础路径 ---
PROJECT_BASE="/share/org/YZWL/yzwl_lixg/project/18.小花糖芥"
ANALYSIS_BASE="${PROJECT_BASE}/02.小花糖芥抗感标记"

# --- 2. 输入/输出目录 ---
# 原始数据目录
RAW_DATA_DIR="${PROJECT_BASE}/01.data/clean"

# 本次分析的各个步骤目录
DATA_DIR="${ANALYSIS_BASE}/01.data"
MAPPING_DIR="${ANALYSIS_BASE}/02.mapping"
GATK_DIR="${ANALYSIS_BASE}/03.gatk_joint"
FILTER_DIR="${ANALYSIS_BASE}/04.filtered_snp_indel"

# --- 3. 关键文件 ---
# 参考基因组路径
REF_GENOME="${DATA_DIR}/Erysimum_cheiranthoides_chromosomes.fa"

# --- 4. 工具参数 ---
THREADS_MAPPING=64
THREADS_FILTER=88


# =================================================================
#                   分析流程 (Analysis Pipeline)
#         这部分代码通常不需要修改，除非流程本身发生变化
# =================================================================

echo ">>> Starting SNP/Indel calling pipeline for Erysimum cheiranthoides..."

# --- 步骤 0: 准备工作 ---
echo ">>> Step 0: Creating output directories..."
# 使用 mkdir -p 一次性创建所有需要的输出目录，如果目录已存在也不会报错
mkdir -p "${DATA_DIR}" "${MAPPING_DIR}" "${GATK_DIR}" "${FILTER_DIR}"

echo ">>> Step 0: Copying clean reads..."
cp "${RAW_DATA_DIR}/EcAA_"* "${DATA_DIR}/"
cp "${RAW_DATA_DIR}/KostanzA_"* "${DATA_DIR}/"
echo "--- Data preparation complete."

# 构建索引
bwa index "${REF_GENOME}"

# --- 步骤 1: 测序数据比对 (Mapping) ---
echo ">>> Step 1: Running mapping with biopytools parabricks..."
biopytools parabricks \
  -i "${DATA_DIR}" \
  -o "${MAPPING_DIR}" \
  -r "${REF_GENOME}" \
  -t "${THREADS_MAPPING}" \
  --no-joint-calling
echo "--- Mapping complete."


# --- 步骤 2: 联合变异检测 (Joint Genotyping) ---
echo ">>> Step 2: Running joint genotyping with biopytools gatk-joint..."
biopytools gatk-joint \
    -i "${MAPPING_DIR}"/vcf \
    -o "${GATK_DIR}" \
    -r "${REF_GENOME}"
echo "--- Joint genotyping complete."


# --- 步骤 3: 变异过滤 (Variant Filtering) ---
echo ">>> Step 3: Filtering variants with biopytools filter-snp-indel..."
# 定义输入VCF文件的完整路径
RAW_VCF="${GATK_DIR}/joint_genotyping_raw.vcf.gz"

biopytools filter-snp-indel \
    -i "${RAW_VCF}" \
    -o "${FILTER_DIR}" \
    -t "${THREADS_FILTER}"
echo "--- Filtering complete."


echo "====================================================="
echo ">>> Pipeline finished successfully!"
echo ">>> Final filtered VCF files can be found in: ${FILTER_DIR}"
echo "====================================================="
