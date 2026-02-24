#!/bin/bash
# =============================================================================
#           HiC-Pro 质量控制分析 - 超算提交脚本 (CSUB) - 完整一体化版
# =============================================================================
#
# 使用方法:
# 1. 确认下面的 #CSUB 指令和用户配置部分的参数正确无误。
# 2. 在当前目录下创建 'Outlog' 文件夹: mkdir Outlog
# 3. 检查原始数据是否在 sample 子目录中, 例如: hic_fastq/AT_sample_1/
# 4. 使用命令提交此脚本: csub submit_hic_qc_final.sh
#
# =============================================================================

# --- (A) 作业调度系统指令 (CSUB) ---
#CSUB -L /bin/bash
#CSUB -J HiCPro_QC_Test      # 作业名称
#CSUB -q c01                 # 作业队列
#CSUB -o Outlog/%J.out       # 标准输出日志
#CSUB -e Outlog/%J.error     # 错误输出日志
#CSUB -n 80                  # CPU核心数
#CSUB -R "span[hosts=1]"     # 确保作业在一个节点上运行

# =============================================================================
# --- (B) 计算节点环境设置与诊断 ---
# =============================================================================
echo "INFO: 作业开始于: $(date)"
echo "INFO: 运行于计算节点: $(hostname)"

# # 1. 严格地激活 Shell 环境，这是确保 'conda' 命令可用的关键
# if [ -f "${HOME}/.bashrc" ]; then
#     echo "INFO: Sourcing ${HOME}/.bashrc"
#     source "${HOME}/.bashrc"
# else
#     echo "警告: ~/.bashrc 文件未找到。"
# fi

# # 2. 激活 HiC-Pro 所在的 Conda 环境
# # !! 重要：请确认您的环境名称是 'HiC-Pro_v3.1.0' !!
# CONDA_ENV_NAME="HiC-Pro_v3.1.0"
# echo "INFO: Activating conda environment: ${CONDA_ENV_NAME}"
# conda activate "${CONDA_ENV_NAME}"

# 3. 严格地设置 HiC-Pro 的路径
HICPRO_INSTALL_DIR="/share/org/YZWL/yzwl_lixg/software/HiC-Pro_3.1.0"
export PATH="${HICPRO_INSTALL_DIR}/bin:${HICPRO_INSTALL_DIR}/bin/utils:${PATH}"

# 4. 最终环境诊断，如果失败则提前退出
echo "INFO: --- Final Environment Check ---"
if ! command -v HiC-Pro &> /dev/null || ! command -v digest_genome.py &> /dev/null; then
    echo "错误: HiC-Pro 或其辅助工具未在 PATH 中找到。"
    echo "  - which HiC-Pro:" $(which HiC-Pro || echo "NOT FOUND")
    echo "  - which digest_genome.py:" $(which digest_genome.py || echo "NOT FOUND")
    echo "  - Current PATH: ${PATH}"
    exit 1
fi
echo "INFO: Environment check PASSED. All commands are available."

# =============================================================================
# --- (C) HiC-Pro 分析流程主体 ---
# =============================================================================

set -e # 若有任何命令执行失败，则立即退出脚本
set -o pipefail # 管道中任何一个命令失败，都视为整个管道失败

# --- 用户配置 ---
GENOME_FA="/share/org/YZWL/yzwl_lixg/project/22.拟南芥hic/02.小试/hicpro/genome.fa"
RAW_DATA_DIR="/share/org/YZWL/yzwl_lixg/project/22.拟南芥hic/02.小试/hicpro/hic_fastq"
ANALYSIS_DIR="/share/org/YZWL/yzwl_lixg/project/22.拟南芥hic/02.小试/hicpro"
ENZYME_NAME="MboI"
RE_SEQUENCE="^GATC"
LIGATION_SITE="GATCGATC"
N_CPU=80 # !! 必须与上面的 #CSUB -n 核心数保持一致 !!

# --- 脚本主体 ---
echo "INFO: Entering analysis directory: ${ANALYSIS_DIR}"
cd "${ANALYSIS_DIR}"

REF_DIR="./reference"
HIC_RESULTS_DIR="./hic_results"
TMP_DIR="./tmp"
CONFIG_FILE="./config-hicpro.txt"

mkdir -p "${REF_DIR}" "${HIC_RESULTS_DIR}" "${TMP_DIR}"

GENOME_FA=$(realpath "${GENOME_FA}")
RAW_DATA_DIR=$(realpath "${RAW_DATA_DIR}")

echo "INFO: --- 开始准备参考基因组文件 ---"
GENOME_PREFIX="${REF_DIR}/genome"
if [ -f "${GENOME_PREFIX}.1.bt2" ]; then
    echo "INFO: Bowtie2索引已存在, 跳过。"
else
    echo "INFO: 正在构建Bowtie2索引..."
    bowtie2-build --threads "${N_CPU}" "${GENOME_FA}" "${GENOME_PREFIX}"
fi

CHROM_SIZES="${REF_DIR}/genome.chrom.sizes"
if [ -f "${CHROM_SIZES}" ]; then
    echo "INFO: 染色体大小文件已存在, 跳过。"
else
    echo "INFO: 正在生成染色体大小文件..."
    samtools faidx "${GENOME_FA}"
    awk -v OFS='\t' '{print $1, $2}' "${GENOME_FA}.fai" > "${CHROM_SIZES}"
fi

GENOME_FRAGMENTS="${REF_DIR}/genome_${ENZYME_NAME}.bed"
if [ -f "${GENOME_FRAGMENTS}" ]; then
    echo "INFO: 基因组酶切片段文件已存在, 跳过。"
else
    echo "INFO: 正在使用 '${ENZYME_NAME}' 酶切基因组..."
    digest_genome.py -r "${RE_SEQUENCE}" -o "${GENOME_FRAGMENTS}" "${GENOME_FA}"
fi
echo "INFO: --- 参考基因组文件准备完毕 ---"

echo "INFO: 正在自动生成HiC-Pro配置文件: ${CONFIG_FILE}"
cat << EOF > "${CONFIG_FILE}"
# ==========================================================
#        HiC-Pro 配置文件 (由脚本自动生成)
# ==========================================================
TMP_DIR = ${TMP_DIR}
LOGS_DIR = ${HIC_RESULTS_DIR}/logs
BOWTIE2_OUTPUT_DIR = ${HIC_RESULTS_DIR}/bowtie_results
MAPC_OUTPUT = ${HIC_RESULTS_DIR}
RAW_DIR = ${RAW_DATA_DIR}

N_CPU = ${N_CPU}
LOGFILE = ${ANALYSIS_DIR}/hicpro_main.log

PAIR1_EXT = _1.fastq.gz
PAIR2_EXT = _2.fastq.gz

MIN_MAPQ = 10
BOWTIE2_IDX_PATH = ${REF_DIR}
BOWTIE2_GLOBAL_OPTIONS = --very-sensitive -L 30 --score-min L,-0.6,-0.2 --end-to-end --reorder
BOWTIE2_LOCAL_OPTIONS =  --very-sensitive -L 20 --score-min L,-0.6,-0.2 --end-to-end --reorder

REFERENCE_GENOME = genome
GENOME_SIZE = ${CHROM_SIZES}

GENOME_FRAGMENT = ${GENOME_FRAGMENTS}
LIGATION_SITE = ${LIGATION_SITE}

GET_ALL_INTERACTION_CLASSES = 1
RM_SINGLETON = 1
RM_MULTI = 1
RM_DUP = 1

BIN_SIZE = 100000 50000
MATRIX_FORMAT = upper

MAX_ITER = 100
FILTER_LOW_COUNT_PERC = 0.02
FILTER_HIGH_COUNT_PERC = 0
EPS = 0.1
EOF
echo "INFO: 配置文件创建成功。"

echo "INFO: --- 开始运行 HiC-Pro 分析流程 (单作业，多线程模式) ---"
echo "INFO: 原始数据目录: ${RAW_DATA_DIR}"
echo "INFO: 结果输出目录: ${HIC_RESULTS_DIR}"
echo "INFO: 开始时间: $(date)"

# !! 关键: 确保这里没有 -p 参数 !!
HiC-Pro -i "${RAW_DATA_DIR}" \
        -o "${HIC_RESULTS_DIR}" \
        -c "${CONFIG_FILE}"

echo "INFO: --- HiC-Pro 分析流程结束 ---"
echo "INFO: 结束时间: $(date)"


# =============================================================================
# --- (D) 解析结果并生成最终QC报告 ---
# =============================================================================
echo "INFO: --- 正在生成QC总结报告 ---"

MAP_STATS_FILE=$(find "${HIC_RESULTS_DIR}/hic_logs/" -name "*.mstat" | head -n 1)
BOWTIE_SUMMARY_FILE="${HIC_RESULTS_DIR}/bowtie_results/bwt2_global_summary.txt"

if [ ! -f "${MAP_STATS_FILE}" ] || [ ! -f "${BOWTIE_SUMMARY_FILE}" ]; then
    echo "错误: 关键统计文件未找到, HiC-Pro可能运行失败。"
    exit 1
fi

TOTAL_READS_BWT2=$(grep "Total pairs" "${BOWTIE_SUMMARY_FILE}" | awk '{print $4}')
MAPPED_PERCENT=$(grep "Overall alignment rate" "${BOWTIE_SUMMARY_FILE}" | awk '{print $1}' | tr -d '%')
TOTAL_PAIRS_PROC=$(grep "total_pairs_processed" "${MAP_STATS_FILE}" | awk '{print $2}')
VALID_INTERACTIONS=$(grep "valid_interaction" "${MAP_STATS_FILE}" | awk '{print $2}')
VALID_INTERACTIONS_NODUPS=$(grep "valid_interaction_rmdup" "${MAP_STATS_FILE}" | awk '{print $2}')
CIS_INTERACTIONS=$(grep "cis_interaction" "${MAP_STATS_FILE}" | awk '{print $2}')
TRANS_INTERACTIONS=$(grep "trans_interaction" "${MAP_STATS_FILE}" | awk '{print $2}')
CIS_SHORT=$(grep "cis_shortRange" "${MAP_STATS_FILE}" | awk '{print $2}')
CIS_LONG=$(grep "cis_longRange" "${MAP_STATS_FILE}" | awk '{print $2}')

VALID_PERCENT=$(awk "BEGIN {if ($TOTAL_PAIRS_PROC>0) printf \"%.2f\", $VALID_INTERACTIONS*100/$TOTAL_PAIRS_PROC; else print 0}")
DUPLICATION_PERCENT=$(awk "BEGIN {if ($VALID_INTERACTIONS>0) printf \"%.2f\", (1 - $VALID_INTERACTIONS_NODUPS/$VALID_INTERACTIONS)*100; else print 0}")
TOTAL_CIS_TRANS=$(awk "BEGIN {if ($CIS_INTERACTIONS+$TRANS_INTERACTIONS}")
CIS_PERCENT=$(awk "BEGIN {if ($TOTAL_CIS_TRANS>0) printf \"%.2f\", $CIS_INTERACTIONS*100/$TOTAL_CIS_TRANS; else print 0}")
TOTAL_CIS=$(awk "BEGIN {print $CIS_SHORT+$CIS_LONG}")
CIS_LONG_PERCENT=$(awk "BEGIN {if ($TOTAL_CIS>0) printf \"%.2f\", $CIS_LONG*100/$TOTAL_CIS; else print 0}")

get_eval() {
    value=$(echo "$1" | awk '{print int($1)}'); good_thr=$2; ok_thr=$3;
    if (( value >= good_thr )); then echo "优秀 (EXCELLENT)";
    elif (( value >= ok_thr )); then echo "可接受 (ACCEPTABLE)";
    else echo "较差 (POOR) - 需重点关注"; fi
}

MAPPED_EVAL=$(get_eval "$MAPPED_PERCENT" 90 80)
VALID_EVAL=$(get_eval "$VALID_PERCENT" 20 10)
CIS_EVAL=$(get_eval "$CIS_PERCENT" 80 60)
CIS_LONG_EVAL=$(get_eval "$CIS_LONG_PERCENT" 40 20)

REPORT_FILE="${ANALYSIS_DIR}/QC_SUMMARY_REPORT.txt"
{
    echo "================================================================";
    echo "             Hi-C 小试数据质量控制 (QC) 总结报告                ";
    echo "================================================================";
    echo "报告生成时间: $(date)";
    echo "原始数据目录: ${RAW_DATA_DIR}";
    echo "----------------------------------------------------------------"; echo "";
    echo "[1] 比对统计 (Mapping Statistics)";
    echo "    - 总 Reads 对数:      ${TOTAL_READS_BWT2}";
    echo "    - 整体比对率:         ${MAPPED_PERCENT}%  [评估: ${MAPPED_EVAL}]";
    echo "      (评估标准: >80% 可接受, >90% 优秀)"; echo "";
    echo "[2] 有效交互统计 (Valid Interaction Statistics)";
    echo "    - 进入Hi-C分析流程的总 Reads 对数: ${TOTAL_PAIRS_PROC}";
    echo "    - 去重后最终有效交互对 (Valid Pairs): ${VALID_INTERACTIONS_NODUPS}";
    echo "    - 有效交互率 (Valid Interaction Rate): ${VALID_PERCENT}%  [评估: ${VALID_EVAL}]";
    echo "      (【核心指标】评估标准: >10% 可接受, >20% 优秀)"; echo "";
    echo "[3] 文库质量评估 (Library Quality Metrics)";
    echo "    - PCR 重复率 (Duplication Rate):    ${DUPLICATION_PERCENT}%";
    echo "      (此值越低越好，代表文库复杂度越高)"; echo "";
    echo "    - Cis 交互占比 (顺式，同一染色体内部): ${CIS_PERCENT}%  [评估: ${CIS_EVAL}]";
    echo "      (Cis应占主导。评估标准: >60% 可接受, >80% 优秀)"; echo "";
    echo "    - Cis 长程交互占比 (Cis > 20kb): ${CIS_LONG_PERCENT}%  [评估: ${CIS_LONG_EVAL}]";
    echo "      (此值越高，捕获远距离互作能力越强。评估标准: >20% 可接受, >40% 优秀)"; echo "";
    echo "----------------------------------------------------------------";
    echo "最终决策建议:";
    if (( $(echo "$VALID_PERCENT >= 10" | bc -l) )) && (( $(echo "$CIS_PERCENT >= 60" | bc -l) )); then
        echo "  [建议: GO] 数据质量“可接受”或“优秀”。可以安排后续测序。";
    elif (( $(echo "$VALID_PERCENT >= 5" | bc -l) )); then
        echo "  [建议: DISCUSS] 数据质量处于“临界状态”。建议与测序公司讨论。";
    else
        echo "  [建议: NO-GO] 数据质量“较差”(<5%)，【不建议】进行大规模测序。";
    fi;
    echo "================================================================";
    echo "完整的 HiC-Pro 结果位于: ${HIC_RESULTS_DIR}";
} | tee "${REPORT_FILE}"

echo "INFO: QC总结报告已保存至 ${REPORT_FILE}"
echo "INFO: 流程运行完毕。"
echo "INFO: 作业结束于: $(date)"