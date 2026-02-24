#!/bin/bash
# =============================================================================
#       Juicer "单机版" 分析流程 - 超算提交脚本 (CSUB) - 最终自动化版
# =============================================================================

# =============================================================================
# --- 计算节点环境与用户配置 ---
# =============================================================================
echo "INFO: 作业开始于: $(date)"
echo "INFO: 运行于计算节点: $(hostname)"

# 1. 环境激活 (!! 非常重要 !!)
#    请根据您的实际情况，取消注释并修改以下行，以确保
#    java, bwa, samtools, python 等命令可用。
#    -----------------------------------------------------------------
#    # >> 方案一: 如果您使用 Conda 环境
#    source /share/org/YZWL/yzwl_lixg/miniforge3/etc/profile.d/conda.sh
#    conda activate your_juicer_env_name # <--- 请替换为您的环境名
#
#    # >> 方案二: 如果您使用 module 系统
#    # module load java/1.8
#    # module load bwa/0.7.17
#    # module load samtools/1.12
#    # module load python/3.8
#    -----------------------------------------------------------------


# 2. 用户配置 (已根据您的信息填写)
JUICER_DIR="/share/org/YZWL/yzwl_lixg/software/juicer"
TOP_DIR="/share/org/YZWL/yzwl_lixg/project/22.拟南芥hic/02.小试/juicer"
GENOME_ID="Est1" # 为您的基因组起一个唯一的ID
RESTRICTION_ENZYME="MboI"
N_THREADS=64 # 必须与上面的 #CSUB -n 核心数保持一致

# =============================================================================
# --- (C) 自动化文件准备 ---
# =============================================================================
echo "INFO: --- 开始自动化文件准备 ---"
# 脚本将在 TOP_DIR (您的工作路径) 下运行
cd "${TOP_DIR}"

# 1. 定义 Juicer 需要的各种文件和目录的路径
GENOME_FA_FILE="${TOP_DIR}/genome.fa"
CHROM_SIZES_FILE="${TOP_DIR}/references/${GENOME_ID}.chrom.sizes"
SITE_FILE="${TOP_DIR}/restriction_sites/${GENOME_ID}_${RESTRICTION_ENZYME}.txt"
BWA_INDEX_CHECK_FILE="${GENOME_FA_FILE}.bwt"

# 2. 创建 Juicer 必需的子目录
mkdir -p references restriction_sites aligned splits
echo "INFO: 目录结构检查完毕。"

# 3. 检查并构建 BWA 索引
if [ -f "${BWA_INDEX_CHECK_FILE}" ]; then
    echo "INFO: BWA 索引已存在, 跳过构建。"
else
    echo "INFO: BWA 索引未找到, 正在构建..."
    bwa index "${GENOME_FA_FILE}"
    echo "INFO: BWA 索引构建完成。"
fi

# 4. 准备染色体大小文件
if [ -f "${CHROM_SIZES_FILE}" ]; then
    echo "INFO: 染色体大小文件已存在, 跳过。"
else
    echo "INFO: 正在生成染色体大小文件..."
    samtools faidx "${GENOME_FA_FILE}"
    awk -v OFS='\t' '{print $1, $2}' "${GENOME_FA_FILE}.fai" > "${CHROM_SIZES_FILE}"
    echo "INFO: 染色体大小文件已生成。"
fi

# 5. 准备酶切位点文件
if [ -f "${SITE_FILE}" ]; then
    echo "INFO: 酶切位点文件已存在, 跳过。"
else
    echo "INFO: 正在生成酶切位点文件..."
    python "${JUICER_DIR}/misc/generate_site_positions.py" ${RESTRICTION_ENZYME} ${GENOME_ID} "${GENOME_FA_FILE}" | sort -k1,1V -k2,2n > "${SITE_FILE}"
    echo "INFO: 酶切位点文件已生成。"
fi
echo "INFO: --- 文件准备工作完成 ---"

# =============================================================================
# --- (D) 启动 Juicer 主流程 ---
# =============================================================================
echo "INFO: --- 准备启动 Juicer (单作业模式) ---"
JUICER_SH="${JUICER_DIR}/scripts/juicer.sh"

if [ ! -f "${JUICER_SH}" ]; then
    echo "错误: Juicer主脚本未找到: ${JUICER_SH}"
    exit 1
fi

echo "INFO: 即将以 ${N_THREADS} 个线程运行 Juicer..."

# 运行 Juicer 主脚本
# 所有必需的路径都已通过参数明确提供
bash ${JUICER_SH} \
    -g ${GENOME_ID} \
    -d ${TOP_DIR} \
    -s ${RESTRICTION_ENZYME} \
    -p ${CHROM_SIZES_FILE} \
    -y ${SITE_FILE} \
    -z ${GENOME_FA_FILE} \
    -D ${JUICER_DIR} \
    -t ${N_THREADS}

echo "INFO: Juicer 流程执行完毕。"
echo "INFO: 作业结束于: $(date)"
