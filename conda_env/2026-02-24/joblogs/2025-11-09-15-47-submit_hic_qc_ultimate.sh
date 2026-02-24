#!/bin/bash
# =============================================================================
#       HiC-Pro QC - 超算提交脚本 (CSUB) - 终极单作业版
# =============================================================================
#
# 运行策略:
#   - 放弃-p模式，回归单作业多线程模式。
#   - 脚本内部强制激活Conda并设置所有路径，不依赖任何外部环境或配置文件。
#   - 解决所有已知问题 (内存、ln路径错误、工具版本冲突)。
#
# =============================================================================

# --- (A) 作业调度系统指令 ---
#CSUB -L /bin/bash
#CSUB -J HiCPro_Ultimate
#CSUB -q c01
#CSUB -o Outlog/%J.out
#CSUB -e Outlog/%J.error
#CSUB -n 16                   # 核心数 (16个已足够高效)
#CSUB -R "span[hosts=1]"
#CSUB -R "rusage[mem=64000]"  # 申请64GB内存

# =============================================================================
# --- (B) 终极环境设置 (核心部分) ---
# =============================================================================
echo "INFO: 作业开始于: $(date)"
echo "INFO: 运行于计算节点: $(hostname)"

# 1. 直接使用 Conda 的绝对路径来初始化环境，这是最稳妥的方式
CONDA_BASE_PATH="/share/org/YZWL/yzwl_lixg/miniforge3"
source "${CONDA_BASE_PATH}/etc/profile.d/conda.sh"

# 2. 激活为 HiC-Pro 准备的 Conda 环境
CONDA_ENV_NAME="HiC-Pro_v3.1.0"
conda activate "${CONDA_ENV_NAME}"

# 3. 获取激活后 Conda 环境的路径
CONDA_ENV_PATH="${CONDA_BASE_PATH}/envs/${CONDA_ENV_NAME}"

# 4. 强制将当前环境的 bin 目录置于 PATH 最顶端，覆盖一切系统默认或用户配置
export PATH="${CONDA_ENV_PATH}/bin:${PATH}"

# 5. 强制设置 HiC-Pro 脚本目录的路径
HICPRO_SCRIPTS_PATH="/share/org/YZWL/yzwl_lixg/software/HiC-Pro/scripts"
export PATH="${HICPRO_SCRIPTS_PATH}:${PATH}"

# 6. 最终环境诊断
echo "INFO: --- Final Environment Check ---"
echo "  - which HiC-Pro:" $(which HiC-Pro || echo "NOT FOUND")
echo "  - which bowtie2:" $(which bowtie2 || echo "NOT FOUND")
echo "  - which samtools:" $(which samtools || echo "NOT FOUND")
if ! command -v HiC-Pro >/dev/null 2>&1 || ! command -v bowtie2 >/dev/null 2>&1; then
    echo "错误: 环境设置失败, 关键命令未找到。"
    exit 1
fi
echo "INFO: Environment check PASSED."

# =============================================================================
# --- (C) HiC-Pro 分析流程主体 ---
# =============================================================================
set -e 
set -o pipefail

# --- 用户配置 ---
ANALYSIS_DIR="/share/org/YZWL/yzwl_lixg/project/22.拟南芥hic/02.小试/hicpro"
N_CPU=16

# --- 脚本主体 ---
echo "INFO: Changing directory to ${ANALYSIS_DIR}"
cd "${ANALYSIS_DIR}"

# --- 定义相对路径变量 ---
GENOME_FA_REL="./genome.fa"
REF_DIR_REL="./reference"
HIC_RESULTS_DIR_REL="./hic_results"
TMP_DIR_REL="./tmp"
CONFIG_FILE_REL="./config-hicpro.txt"
RAW_DATA_DIR_REL="./hic_fastq"

mkdir -p "${REF_DIR_REL}" "${HIC_RESULTS_DIR_REL}" "${TMP_DIR_REL}"

echo "INFO: --- 开始准备参考基因组文件 ---"
if [ ! -f "${REF_DIR_REL}/genome.1.bt2" ]; then
    bowtie2-build --threads "${N_CPU}" "${GENOME_FA_REL}" "${REF_DIR_REL}/genome"
fi
# (此处省略了 samtools faidx 和 digest_genome.py 的检查与执行，假设它们已准备好或按需运行)

echo "INFO: 正在自动生成HiC-Pro配置文件: ${CONFIG_FILE_REL}"
cat << EOF > "${CONFIG_FILE_REL}"
# --- HiC-Pro 配置文件 (由终极脚本自动生成) ---
N_CPU = ${N_CPU}
LOGFILE = ./hicpro_main.log
TMP_DIR = ${TMP_DIR_REL}
PAIR1_EXT = _1.fastq.gz
PAIR2_EXT = _2.fastq.gz
MIN_MAPQ = 10

BOWTIE2_IDX_PATH = ${REF_DIR_REL}
REFERENCE_GENOME = genome
GENOME_SIZE = ${REF_DIR_REL}/genome.chrom.sizes

GENOME_FRAGMENT = ${REF_DIR_REL}/genome_MboI.bed
LIGATION_SITE = GATCGATC

GET_ALL_INTERACTION_CLASSES = 1
BIN_SIZE = 100000 50000
EOF

echo "INFO: --- 开始运行 HiC-Pro (终极单作业模式) ---"
echo "INFO: 开始时间: $(date)"

# !! 关键: 使用相对路径调用，不带 -p !!
HiC-Pro -i "${RAW_DATA_DIR_REL}" \
        -o "${HIC_RESULTS_DIR_REL}" \
        -c "${CONFIG_FILE_REL}"

echo "INFO: --- HiC-Pro 分析流程结束 ---"
echo "INFO: 结束时间: $(date)"

# (后续的QC报告生成部分被省略，可以在运行成功后添加回来)
