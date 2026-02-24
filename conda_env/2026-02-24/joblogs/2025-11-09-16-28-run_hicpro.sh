#!/bin/bash
# =============================================================================
#           HiC-Pro QC - 最终版 (原生脚本 + Singularity 核心)
# =============================================================================

set -e # 若有任何命令执行失败，则立即退出脚本
set -o pipefail # 管道中任何一个命令失败，都视为整个管道失败

# --- (B) 用户配置 ---
# Singularity 镜像文件的绝对路径
IMAGE_PATH="/share/org/YZWL/yzwl_lixg/software/singularity/hicpro_3.1.0_ubuntu.img"
# 主分析目录的绝对路径
ANALYSIS_DIR="/share/org/YZWL/yzwl_lixg/project/22.拟南芥hic/02.小试/hicpro"
# CPU 核心数
N_CPU=80

# --- (C) 脚本主体 ---

echo "INFO: Changing directory to ${ANALYSIS_DIR}"
cd "${ANALYSIS_DIR}"

# --- 定义相对路径 ---
REF_DIR_REL="./reference"
HIC_RESULTS_DIR_REL="./hic_results"
TMP_DIR_REL="./tmp"
CONFIG_FILE_REL="./config-hicpro.txt"
RAW_DATA_DIR_REL="./hic_fastq"

mkdir -p "${REF_DIR_REL}" "${TMP_DIR_REL}"

# --- 将所有核心计算步骤封装在一个 'singularity exec' 命令中 ---
echo "INFO: --- Entering Singularity container to perform core computations ---"

singularity exec --bind ${ANALYSIS_DIR}:${ANALYSIS_DIR} ${IMAGE_PATH} bash -c '
set -e
set -o pipefail

# --- 从这里开始，所有命令都在 Singularity 容器内部执行 ---
echo "INFO: (内部) 成功进入 Singularity 容器环境。"
echo "INFO: (内部) 切换到工作目录: '${ANALYSIS_DIR}'"
cd "'${ANALYSIS_DIR}'"

# --- 定义内部变量 ---
GENOME_FA_REL="./genome.fa"
REF_DIR_REL="./reference"
HIC_RESULTS_DIR_REL="./hic_results"
TMP_DIR_REL="./tmp"
CONFIG_FILE_REL="./config-hicpro.txt"
RAW_DATA_DIR_REL="./hic_fastq"

echo "INFO: (内部) --- 开始准备参考基因组文件 ---"
if [ -f "${REF_DIR_REL}/genome.1.bt2" ]; then
    echo "INFO: (内部) Bowtie2索引已存在, 跳过。"
else
    echo "INFO: (内部) 正在构建Bowtie2索引..."
    bowtie2-build --threads '${N_CPU}' "${GENOME_FA_REL}" "${REF_DIR_REL}/genome"
fi
if [ ! -f "${REF_DIR_REL}/genome.chrom.sizes" ]; then
    samtools faidx "${GENOME_FA_REL}"
    awk -v OFS='\''\t'\'' '\''{print $1, $2}'\'' "${GENOME_FA_REL}.fai" > "${REF_DIR_REL}/genome.chrom.sizes"
fi
if [ ! -f "${REF_DIR_REL}/genome_MboI.bed" ]; then
    digest_genome.py -g "${GENOME_FA_REL}" -r "^GATC" -o "${REF_DIR_REL}/genome_MboI.bed"
fi
echo "INFO: (内部) --- 参考基因组文件准备完毕 ---"

echo "INFO: (内部) 正在自动生成HiC-Pro配置文件: ${CONFIG_FILE_REL}"
cat << EOF > "${CONFIG_FILE_REL}"
# Hi-C Pro 配置文件 (由脚本自动生成)
TMP_DIR = '${TMP_DIR_REL}'
LOGS_DIR = '${HIC_RESULTS_DIR_REL}'/logs
BOWTIE2_OUTPUT_DIR = '${HIC_RESULTS_DIR_REL}'/bowtie_results
MAPC_OUTPUT = '${HIC_RESULTS_DIR_REL}'
RAW_DIR = '${RAW_DATA_DIR_REL}'
N_CPU = '${N_CPU}'
LOGFILE = ./hicpro_main.log
PAIR1_EXT = _1
PAIR2_EXT = _2
MIN_MAPQ = 10
BOWTIE2_IDX_PATH = '${REF_DIR_REL}'
REFERENCE_GENOME = genome
GENOME_SIZE = '${REF_DIR_REL}'/genome.chrom.sizes
GENOME_FRAGMENT = '${REF_DIR_REL}'/genome_MboI.bed
LIGATION_SITE = GATCGATC
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
echo "INFO: (内部) 配置文件创建成功。"

echo "INFO: (内部) --- 开始运行 HiC-Pro 分析流程 ---"
HiC-Pro -i "${RAW_DATA_DIR_REL}" -o "${HIC_RESULTS_DIR_REL}" -c "${CONFIG_FILE_REL}"
echo "INFO: (内部) --- HiC-Pro 分析流程结束 ---"

' # "bash -c" 命令结束

echo "INFO: --- Core computations finished, exiting Singularity container ---"

# =============================================================================
# --- (D) 解析结果并生成最终QC报告 (在原生环境运行) ---
# =============================================================================
echo "INFO: --- 正在生成QC总结报告 ---"

# QC报告部分现在也自然地在当前目录下工作
MAP_STATS_FILE=$(find "${HIC_RESULTS_DIR_REL}/hic_logs/" -name "*.mstat" | head -n 1)
BOWTIE_SUMMARY_FILE="${HIC_RESULTS_DIR_REL}/bowtie_results/bwt2_global_summary.txt"

# ... (此处省略和您脚本中完全一样的QC报告生成代码) ...
if [ ! -f "${MAP_STATS_FILE}" ] || [ ! -f "${BOWTIE_SUMMARY_FILE}" ]; then echo "错误: 关键统计文件未找到, HiC-Pro可能运行失败。"; exit 1; fi; TOTAL_READS_BWT2=$(grep "Total pairs" "${BOWTIE_SUMMARY_FILE}" | awk '{print $4}'); MAPPED_PERCENT=$(grep "Overall alignment rate" "${BOWTIE_SUMMARY_FILE}" | awk '{print $1}' | tr -d '%'); TOTAL_PAIRS_PROC=$(grep "total_pairs_processed" "${MAP_STATS_FILE}" | awk '{print $2}'); VALID_INTERACTIONS=$(grep "valid_interaction" "${MAP_STATS_FILE}" | awk '{print $2}'); VALID_INTERACTIONS_NODUPS=$(grep "valid_interaction_rmdup" "${MAP_STATS_FILE}" | awk '{print $2}'); CIS_INTERACTIONS=$(grep "cis_interaction" "${MAP_STATS_FILE}" | awk '{print $2}'); TRANS_INTERACTIONS=$(grep "trans_interaction" "${MAP_STATS_FILE}" | awk '{print $2}'); CIS_SHORT=$(grep "cis_shortRange" "${MAP_STATS_FILE}" | awk '{print $2}'); CIS_LONG=$(grep "cis_longRange" "${MAP_STATS_FILE}" | awk '{print $2}'); VALID_PERCENT=$(awk "BEGIN {if ($TOTAL_PAIRS_PROC>0) printf \"%.2f\", $VALID_INTERACTIONS*100/$TOTAL_PAIRS_PROC; else print 0}"); DUPLICATION_PERCENT=$(awk "BEGIN {if ($VALID_INTERACTIONS>0) printf \"%.2f\", (1 - $VALID_INTERACTIONS_NODUPS/$VALID_INTERACTIONS)*100; else print 0}"); TOTAL_CIS_TRANS=$(awk "BEGIN {print $CIS_INTERACTIONS+$TRANS_INTERACTIONS}"); CIS_PERCENT=$(awk "BEGIN {if ($TOTAL_CIS_TRANS>0) printf \"%.2f\", $CIS_INTERACTIONS*100/$TOTAL_CIS_TRANS; else print 0}"); TOTAL_CIS=$(awk "BEGIN {print $CIS_SHORT+$CIS_LONG}"); CIS_LONG_PERCENT=$(awk "BEGIN {if ($TOTAL_CIS>0) printf \"%.2f\", $CIS_LONG*100/$TOTAL_CIS; else print 0}"); get_eval() { value=$(echo "$1" | awk '{print int($1)}'); good_thr=$2; ok_thr=$3; if (( value >= good_thr )); then echo "优秀 (EXCELLENT)"; elif (( value >= ok_thr )); then echo "可接受 (ACCEPTABLE)"; else echo "较差 (POOR) - 需重点关注"; fi }; MAPPED_EVAL=$(get_eval "$MAPPED_PERCENT" 90 80); VALID_EVAL=$(get_eval "$VALID_PERCENT" 20 10); CIS_EVAL=$(get_eval "$CIS_PERCENT" 80 60); CIS_LONG_EVAL=$(get_eval "$CIS_LONG_PERCENT" 40 20); REPORT_FILE="./QC_SUMMARY_REPORT.txt"; { echo "================================================================"; echo "             Hi-C 小试数据质量控制 (QC) 总结报告                "; echo "================================================================"; echo "报告生成时间: $(date)"; echo "原始数据目录: ${ANALYSIS_DIR}/hic_fastq"; echo "----------------------------------------------------------------"; echo ""; echo "[1] 比对统计 (Mapping Statistics)"; echo "    - 总 Reads 对数:      ${TOTAL_READS_BWT2}"; echo "    - 整体比对率:         ${MAPPED_PERCENT}%  [评估: ${MAPPED_EVAL}]"; echo "      (评估标准: >80% 可接受, >90% 优秀)"; echo ""; echo "[2] 有效交互统计 (Valid Interaction Statistics)"; echo "    - 进入Hi-C分析流程的总 Reads 对数: ${TOTAL_PAIRS_PROC}"; echo "    - 去重后最终有效交互对 (Valid Pairs): ${VALID_INTERACTIONS_NODUPS}"; echo "    - 有效交互率 (Valid Interaction Rate): ${VALID_PERCENT}%  [评估: ${VALID_EVAL}]"; echo "      (【核心指标】评估标准: >10% 可接受, >20% 优秀)"; echo ""; echo "[3] 文库质量评估 (Library Quality Metrics)"; echo "    - PCR 重复率 (Duplication Rate):    ${DUPLICATION_PERCENT}%"; echo "      (此值越低越好，代表文库复杂度越高)"; echo ""; echo "    - Cis 交互占比 (顺式，同一染色体内部): ${CIS_PERCENT}%  [评估: ${CIS_EVAL}]"; echo "      (Cis应占主导。评估标准: >60% 可接受, >80% 优秀)"; echo ""; echo "    - Cis 长程交互占比 (Cis > 20kb): ${CIS_LONG_PERCENT}%  [评估: ${CIS_LONG_EVAL}]"; echo "      (此值越高，捕获远距离互作能力越强。评估标准: >20% 可接受, >40% 优秀)"; echo ""; echo "----------------------------------------------------------------"; echo "最终决策建议:"; if (( $(echo "$VALID_PERCENT >= 10" | bc -l) )) && (( $(echo "$CIS_PERCENT >= 60" | bc -l) )); then echo "  [建议: GO] 数据质量“可接受”或“优秀”。可以安排后续测序。"; elif (( $(echo "$VALID_PERCENT >= 5" | bc -l) )); then echo "  [建议: DISCUSS] 数据质量处于“临界状态”。建议与测序公司讨论。"; else echo "  [建议: NO-GO] 数据质量“较差”(<5%)，【不建议】进行大规模测序。"; fi; echo "================================================================"; echo "完整的 HiC-Pro 结果位于: ${ANALYSIS_DIR}/hic_results"; } | tee "${REPORT_FILE}";

echo "INFO: QC总结报告已保存至 ${REPORT_FILE}"
echo "INFO: 流程运行完毕。"
echo "INFO: 作业结束于: $(date)"