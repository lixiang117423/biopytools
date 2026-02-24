#!/bin/bash
# ------------------------------------------------------------------
# 脚本名称: run_full_pipeline.sh
# 功能: 完整流程 Mapping -> Joint Call -> Filter -> Genotype -> Stat
# ------------------------------------------------------------------

# ================= 1. 配置区域 (只需修改这里) =================
set -e  # 遇到错误立即停止，避免产生错误文件

# --- [1] 项目参数 ---
SPECIES="hancai"       # 结果文件前缀 (例如: hancai, rice, soybean)
THREADS=64             # 通用线程数 (GATK, biopytools)
MAPPING_THREADS=64     # 比对步骤专用线程数 (BWA通常可以开高一点)
BCF_THREADS=64         # bcftools压缩线程数

# --- [2] 环境设置 ---
# 添加 GATK 到环境变量
export PATH="/share/org/YZWL/yzwl_lixg/miniforge3/envs/GATK_v.4.6.2.0/bin/:$PATH"

# --- [3] 外部脚本路径 (请确认这些脚本真实存在) ---
# Mapping 脚本路径
MAPPING_SCRIPT="${HOME}/software/scripts/32.run_bwa_gatk.sh"
# 统计脚本路径
BAM_STAT_SCRIPT="${HOME}/software/scripts/23.bam_stats_reporter.py"

# ================= 2. 自动构建绝对路径 =================
# 获取当前所在目录的绝对路径
PROJECT_DIR=$(pwd)

echo ">>> [初始化] 正在构建绝对路径..."
echo ">>> 项目根目录: ${PROJECT_DIR}"
echo ">>> 物种前缀: ${SPECIES}"

# 定义关键文件的【绝对路径】
REF_GENOME="${PROJECT_DIR}/01.data/genome/genome.fa"
CLEAN_DATA_DIR="${PROJECT_DIR}/01.data/clean"

# 定义输出目录的【绝对路径】
OUT_MAPPING="${PROJECT_DIR}/02.mapping"
OUT_BAM="${OUT_MAPPING}/01.bam"
OUT_GVCF="${OUT_MAPPING}/02.gvcf"
OUT_JOINT="${PROJECT_DIR}/03.gatk_joint"
OUT_FILTER="${PROJECT_DIR}/04.filtered_snp_indel"
OUT_GENOTYPE="${PROJECT_DIR}/05.snp_indel_genotype"
OUT_STATS="${OUT_MAPPING}/04.bam_stat"

# ================= 3. 环境与文件检查 =================
echo ">>> [检查] 验证输入文件..."

if [ ! -f "${REF_GENOME}" ]; then
    echo "❌ 错误: 参考基因组不存在: ${REF_GENOME}"
    exit 1
fi

if [ ! -d "${CLEAN_DATA_DIR}" ]; then
    echo "❌ 错误: 测序数据目录不存在: ${CLEAN_DATA_DIR}"
    exit 1
fi

# 创建输出目录结构
mkdir -p "${OUT_MAPPING}"
mkdir -p "${OUT_STATS}"
mkdir -p "${OUT_JOINT}"
mkdir -p "${OUT_FILTER}"
mkdir -p "${OUT_GENOTYPE}/snp"
mkdir -p "${OUT_GENOTYPE}/indel"

# ================= 4. Mapping & Variant Calling (单样本) =================
echo ""
echo "========================================================"
echo "[1/6] Mapping & HaplotypeCaller (BWA + GATK)..."
echo "========================================================"
# 调用外部脚本进行比对，传入绝对路径
# 注意：这里假设 32.run_bwa_gatk.sh 接受 -t 参数来控制线程，如果没有，请手动修改脚本内的线程
# bash "${MAPPING_SCRIPT}" \
#     -i "${CLEAN_DATA_DIR}" \
#     -r "${REF_GENOME}" \
#     -o "${OUT_MAPPING}" 

# 检查 GVCF 是否生成
if [ -z "$(ls -A ${OUT_GVCF}/*.g.vcf.gz 2>/dev/null)" ]; then
    echo "❌ 错误: Mapping 完成，但在 ${OUT_GVCF} 下未检测到 GVCF 文件。"
    exit 1
fi

# ================= 5. GATK Joint Calling (多样本联合) =================
echo ""
echo "========================================================"
echo "[2/6] GATK Joint Calling (联合变异检测)..."
echo "========================================================"

# 清理旧目录，防止路径嵌套 (03.gatk_joint/03.gatk_joint)
if [ -d "${OUT_JOINT}" ]; then
    echo "   -> 清理旧的 Joint Call 输出目录..."
    rm -rf "${OUT_JOINT:?}"/*
fi

# 运行 biopytools (传入绝对路径)
biopytools gatk-joint \
    -i "${OUT_GVCF}" \
    -r "${REF_GENOME}" \
    -o "${OUT_JOINT}"

# 检查联合变异结果
JOINT_VCF="${OUT_JOINT}/joint_genotyping_merged.vcf.gz"
if [ ! -f "${JOINT_VCF}" ]; then
    echo "❌ 错误: Joint Calling 失败，未生成文件: ${JOINT_VCF}"
    exit 1
fi

# ================= 6. 过滤 SNP 和 INDEL =================
echo ""
echo "========================================================"
echo "[3/6] 过滤 SNP 和 INDEL..."
echo "========================================================"

biopytools filter-snp-indel \
    -i "${JOINT_VCF}" \
    -o "${OUT_FILTER}" \
    -t ${THREADS}

# ================= 7. MAF 0.05 筛选 =================
echo ""
echo "========================================================"
echo "[4/6] 执行 MAF > 0.05 筛选..."
echo "========================================================"

INPUT_SNP_VCF="${OUT_FILTER}/variation.filtered.snp.vcf.gz"
OUTPUT_MAF_VCF="${OUT_FILTER}/variation.filtered.snp.maf_05.vcf.gz"

if [ -f "${INPUT_SNP_VCF}" ]; then
    bcftools view \
        -i 'MAF > 0.05' \
        --threads ${BCF_THREADS} \
        -Oz \
        -o "${OUTPUT_MAF_VCF}" \
        "${INPUT_SNP_VCF}"
    
    echo "   -> 建立索引..."
    bcftools index -t "${OUTPUT_MAF_VCF}"
else
    echo "❌ 错误: 未找到 SNP 文件 ${INPUT_SNP_VCF}，无法进行 MAF 筛选。"
    exit 1
fi

# ================= 8. 提取基因型 (CSV) =================
echo ""
echo "========================================================"
echo "[5/6] 提取基因型矩阵 (前缀: ${SPECIES})..."
echo "========================================================"

# 8.1 SNP
echo "   -> 正在处理 SNP..."
biopytools vcf-genotype \
    -i "${OUTPUT_MAF_VCF}" \
    -o "${SPECIES}" \
    --output-dir "${OUT_GENOTYPE}/snp" \
    --each yes \
    -t csv

# 8.2 INDEL
echo "   -> 正在处理 INDEL..."
INDEL_VCF="${OUT_FILTER}/variation.filtered.indel.vcf.gz"
if [ -f "${INDEL_VCF}" ]; then
    biopytools vcf-genotype \
        -i "${INDEL_VCF}" \
        -o "${SPECIES}" \
        --output-dir "${OUT_GENOTYPE}/indel" \
        --each yes \
        -t csv
else
    echo "⚠️ 警告: 未找到 INDEL 文件，跳过此步。"
fi

# ================= 9. 比对结果统计 =================
echo ""
echo "========================================================"
echo "[6/6] 统计 BAM 比对率..."
echo "========================================================"

if [ -f "${BAM_STAT_SCRIPT}" ]; then
    python3 "${BAM_STAT_SCRIPT}" \
        -i "${OUT_BAM}" \
        -o "${OUT_STATS}/bam_stat.xlsx"
else
    echo "⚠️ 警告: 找不到统计脚本 ${BAM_STAT_SCRIPT}，跳过统计。"
fi

# ================= 结束 =================
echo ""
echo ">>> ✅ 所有任务已完成！"