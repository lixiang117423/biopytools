#!/bin/bash
# =============================================================================
#       HiC-Pro 质量控制分析 - 超算提交脚本 (CSUB) - Singularity 最终版
# =============================================================================

# =============================================================================
# --- 用户配置 (这是唯一需要您修改的地方) ---
# =============================================================================

# 1. Singularity 镜像文件的绝对路径 (已根据您的信息更新)
IMAGE_PATH="/share/org/YZWL/yzwl_lixg/software/singularity/hicpro_latest.sif"

# 2. 您的主分析目录的绝对路径
#    !! 重要: 所有输入(genome, fastq)和输出都应在此目录或其子目录中 !!
ANALYSIS_DIR="/share/org/YZWL/yzwl_lixg/project/22.拟南芥hic/02.小试/hicpro"

# 3. CPU 核心数 (必须与上面的 #CSUB -n 核心数保持一致)
N_CPU=80

# =============================================================================
# --- (C) Singularity 执行命令 ---
# =============================================================================
echo "INFO: 作业开始于: $(date)"
echo "INFO: 运行于计算节点: $(hostname)"
echo "INFO: Singularity 镜像: ${IMAGE_PATH}"
echo "INFO: 分析目录: ${ANALYSIS_DIR}"

# 检查镜像文件是否存在
if [ ! -f "${IMAGE_PATH}" ]; then
    echo "错误: Singularity 镜像文件未找到: ${IMAGE_PATH}"
    exit 1
fi

# 'singularity exec' 是核心命令
# --bind 参数是关键: 它将您的主机分析目录 "映射" 到容器内部的同一个路径下
# bash -c '...' 允许我们在容器内执行一长串的命令
singularity exec --bind ${ANALYSIS_DIR}:${ANALYSIS_DIR} ${IMAGE_PATH} bash -c '

set -e
set -o pipefail

# --- 从这里开始，所有命令都在 Singularity 容器内部执行 ---
echo "INFO: (内部) 成功进入 Singularity 容器环境。"

# --- 用户配置 (在容器内部使用) ---
GENOME_FA="'${ANALYSIS_DIR}'/genome.fa"
RAW_DATA_DIR_REL="./hic_fastq"
ANALYSIS_DIR_INTERNAL="'${ANALYSIS_DIR}'" # 容器内的路径
ENZYME_NAME="MboI"
RE_SEQUENCE="^GATC"
LIGATION_SITE="GATCGATC"

# --- 脚本主体 (与之前的版本类似，但现在在容器内运行) ---
echo "INFO: (内部) 切换到工作目录: ${ANALYSIS_DIR_INTERNAL}"
cd ${ANALYSIS_DIR_INTERNAL}

REF_DIR_REL="./reference"
HIC_RESULTS_DIR_REL="./hic_results"
TMP_DIR_REL="./tmp"
CONFIG_FILE_REL="./config-hicpro.txt"

mkdir -p "${REF_DIR_REL}" "${HIC_RESULTS_DIR_REL}" "${TMP_DIR_REL}"

echo "INFO: (内部) --- 开始准备参考基因组文件 ---"
GENOME_PREFIX="${REF_DIR_REL}/genome"
if [ -f "${GENOME_PREFIX}.1.bt2" ]; then
    echo "INFO: (内部) Bowtie2索引已存在, 跳过。"
else
    echo "INFO: (内部) 正在构建Bowtie2索引..."
    bowtie2-build --threads '${N_CPU}' "${GENOME_FA}" "${GENOME_PREFIX}"
fi

CHROM_SIZES="${REF_DIR_REL}/genome.chrom.sizes"
if [ ! -f "${CHROM_SIZES}" ]; then
    samtools faidx "${GENOME_FA}"
    awk -v OFS='\''\t'\'' '\''{print $1, $2}'\'' "${GENOME_FA}.fai" > "${CHROM_SIZES}"
fi

GENOME_FRAGMENTS="${REF_DIR_REL}/genome_${ENZYME_NAME}.bed"
if [ ! -f "${GENOME_FRAGMENTS}" ]; then
    digest_genome.py -g "${GENOME_FA}" -r "'${RE_SEQUENCE}'" -o "${GENOME_FRAGMENTS}"
fi
echo "INFO: (内部) --- 参考基因组文件准备完毕 ---"

echo "INFO: (内部) 正在自动生成HiC-Pro配置文件: ${CONFIG_FILE_REL}"
cat << EOF > "${CONFIG_FILE_REL}"
N_CPU = '${N_CPU}'
LOGFILE = ./hicpro_main.log
TMP_DIR = ${TMP_DIR_REL}
PAIR1_EXT = _1.fastq.gz
PAIR2_EXT = _2.fastq.gz
MIN_MAPQ = 10
BOWTIE2_IDX_PATH = ${REF_DIR_REL}
REFERENCE_GENOME = genome
GENOME_SIZE = ${CHROM_SIZES}
GENOME_FRAGMENT = ${GENOME_FRAGMENTS}
LIGATION_SITE = '${LIGATION_SITE}'
GET_ALL_INTERACTION_CLASSES = 1
BIN_SIZE = 100000 50000
MATRIX_FORMAT = upper
MAX_ITER = 100
FILTER_LOW_COUNT_PERC = 0.02
EOF
echo "INFO: (内部) 配置文件创建成功。"

echo "INFO: (内部) --- 开始运行 HiC-Pro 分析流程 ---"
echo "INFO: (内部) 开始时间: $(date)"

# 使用相对路径调用 HiC-Pro
HiC-Pro -i "${RAW_DATA_DIR_REL}" \
        -o "${HIC_RESULTS_DIR_REL}" \
        -c "${CONFIG_FILE_REL}"

echo "INFO: (内部) --- HiC-Pro 分析流程结束 ---"
echo "INFO: (内部) 结束时间: $(date)"

' # "bash -c" 命令结束

# 检查 Singularity 命令的退出状态
EXIT_CODE=$?
if [ ${EXIT_CODE} -ne 0 ]; then
    echo "错误: Singularity 容器内的命令执行失败，退出码: ${EXIT_CODE}"
    exit ${EXIT_CODE}
fi

echo "INFO: 作业成功结束于: $(date)"