#!/bin/bash
# =============================================================================
#       Juicer "单机版" 分析流程 - 超算提交脚本 (CSUB) - 最终自动化版
# =============================================================================

set -e
set -o pipefail

CONDA_ENV_BIN_PATH="/share/org/YZWL/yzwl_lixg/miniforge3/envs/juicer_v.1.6/bin"
export PATH="${CONDA_ENV_BIN_PATH}:${PATH}"

# =============================================================================
# ---  计算节点环境与用户配置 ---
# =============================================================================
echo "INFO: 作业开始于: $(date)"
echo "INFO: 运行于计算节点: $(hostname)"

# 检查 'aligned' 目录是否存在。如果存在，则将其重命名为带有时间戳的备份
if [ -d "aligned" ]; then
    TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
    echo "警告: 'aligned' 目录已存在。正在将其重命名为 'aligned_backup_${TIMESTAMP}'"
    mv aligned "aligned_backup_${TIMESTAMP}"
fi

# 2. 用户配置
JUICER_DIR="/share/org/YZWL/yzwl_lixg/software/juicer"
TOP_DIR="/share/org/YZWL/yzwl_lixg/project/22.Est1_hic/02.test/juicer"
GENOME_ID="Est1"
RESTRICTION_ENZYME="MboI"
N_THREADS=8

# =============================================================================
# --- (C) 自动化文件准备 (保持不变) ---
# =============================================================================
echo "INFO: --- 开始自动化文件准备 ---"
cd "${TOP_DIR}"

GENOME_FA_FILE="./genome.fa"
CHROM_SIZES_FILE="./references/${GENOME_ID}.chrom.sizes"
SITE_FILE="./restriction_sites/${GENOME_ID}_${RESTRICTION_ENZYME}.txt"
BWA_INDEX_CHECK_FILE="${GENOME_FA_FILE}.bwt"

mkdir -p references restriction_sites splits

if [ -f "${BWA_INDEX_CHECK_FILE}" ]; then
    echo "INFO: BWA 索引已存在, 跳过构建。"
else
    echo "INFO: BWA 索引未找到, 正在构建..."
    bwa index "${GENOME_FA_FILE}"
    echo "INFO: BWA 索引构建完成。"
fi

if [ -f "${CHROM_SIZES_FILE}" ]; then
    echo "INFO: 染色体大小文件已存在, 跳过。"
else
    samtools faidx "${GENOME_FA_FILE}"
    awk -v OFS='\t' '{print $1, $2}' "${GENOME_FA_FILE}.fai" > "${CHROM_SIZES_FILE}"
    echo "INFO: 染色体大小文件已生成。"
fi

if [ -f "${SITE_FILE}" ]; then
    echo "INFO: 酶切位点文件已存在, 跳过。"
else
    python "${JUICER_DIR}/misc/generate_site_positions.py" ${RESTRICTION_ENZYME} ${GENOME_ID} "${GENOME_FA_FILE}" | sort -k1,1V -k2,2n > "${SITE_FILE}"
    echo "INFO: 酶切位点文件已生成。"
fi
echo "INFO: --- 文件准备工作完成 ---"

# =============================================================================
# --- (D) 启动 Juicer 主流程 (保持不变) ---
# =============================================================================
echo "INFO: --- 准备启动 Juicer (单作业模式) ---"
JUICER_SH="${JUICER_DIR}/scripts/juicer.sh"

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

# =============================================================================
# --- (E) 自动生成质量评估报告 (新增功能) ---
# =============================================================================
echo "INFO: --- 开始生成质量评估报告 ---"

# 1. 定义关键文件路径
STATS_FILE="${TOP_DIR}/aligned/inter_30.txt"
HIC_FILE="${TOP_DIR}/aligned/inter_30.hic"
REPORT_FILE="${TOP_DIR}/aligned/QC_SUMMARY_REPORT.md"

if [ ! -f "${STATS_FILE}" ]; then
    echo "错误: Juicer 核心统计文件未找到: ${STATS_FILE}"
    echo "报告生成失败。"
    exit 1
fi

# 2. 从 inter_30.txt 中提取核心统计数据
#    Juicer的统计数据格式是 "指标: 值 (百分比1, 百分比2)"
TOTAL_READS=$(grep "Sequenced Read Pairs" "${STATS_FILE}" | awk '{print $NF}')
UNMAPPED_READS=$(grep "Unmapped" "${STATS_FILE}" | awk '{print $5}' | tr -d '()%,')
CHIMERIC_PAIRS=$(grep "Chimeric Paired" "${STATS_FILE}" | awk '{print $5}' | tr -d '()%,')
NON_CHIMERIC_AMBIGUOUS=$(grep "Chimeric Ambiguous" "${STATS_FILE}" | awk '{print $5}' | tr -d '()%,')
TOTAL_ALIGNED=$(awk -v total=$TOTAL_READS -v unmapped=$UNMAPPED_READS 'BEGIN {printf "%.2f", 100 - unmapped}')
LIGATION_MOTIF=$(grep "Ligation Motif Present" "${STATS_FILE}" | awk '{print $5}' | tr -d '()%,')

# 提取最关键的 Unique 和 Duplicates 数据
UNIQUE_READS_LINE=$(grep "Total Unique" "${STATS_FILE}")
TOTAL_UNIQUE=$(echo ${UNIQUE_READS_LINE} | awk '{print $3}')
UNIQUE_PERCENT_OF_TOTAL=$(echo ${UNIQUE_READS_LINE} | awk '{print $5}' | tr -d '()%,') # 这是黄金指标

DUPS_READS_LINE=$(grep "Total Duplicates" "${STATS_FILE}")
TOTAL_DUPS=$(echo ${DUPS_READS_LINE} | awk '{print $3}')
DUPS_PERCENT_OF_TOTAL=$(echo ${DUPS_READS_LINE} | awk '{print $5}' | tr -d '()%,')

# 3. 生成 Markdown 格式的报告文件
{
    echo "# Juicer Hi-C 数据质量评估报告"
    echo ""
    echo "**报告生成时间:** $(date)"
    echo "**分析目录:** \`${TOP_DIR}\`"
    echo ""
    echo "---"
    echo ""
    echo "## 一、核心结论"
    echo ""
    if (( $(echo "$UNIQUE_PERCENT_OF_TOTAL >= 10" | bc -l) )); then
        echo "> **数据质量“可接受”或“优秀”。** 核心指标——最终有效交互率达到了 **${UNIQUE_PERCENT_OF_TOTAL}%**，超过了进行大规模测序的推荐阈值 (10%)。可以安排后续测序。"
    else
        echo "> **数据质量“较差”，需谨慎评估。** 核心指标——最终有效交互率仅为 **${UNIQUE_PERCENT_OF_TOTAL}%**，低于推荐阈值 (10%)。不建议立即进行大规模测序，请与测序公司讨论实验细节。"
    fi
    echo ""
    echo "---"
    echo ""
    echo "## 二、关键产出文件"
    echo ""
    echo "- **核心统计文件 (MAPQ>=30):** \`aligned/inter_30.txt\`"
    echo "- **Juicebox 可视化文件 (MAPQ>=30):** \`aligned/inter_30.hic\`"
    echo "  - *您可以下载此文件，并使用 Juicebox 桌面软件打开，以交互式地浏览互作热图。*"
    echo ""
    echo "---"
    echo ""
    echo "## 三、详细统计数据解读 (基于 inter_30.txt)"
    echo ""
    echo "| 指标 (Metric) | 数量 / 比例 | 解读 |"
    echo "|:---|:---|:---|"
    echo "| **总测序 Reads 对数 (Sequenced Read Pairs)** | \`${TOTAL_READS}\` | 原始数据总量。 |"
    echo "| **比对率 (Alignment Rate)** | \`${TOTAL_ALIGNED}%\` | Reads 能够成功比对到参考基因组的比例，越高越好。 |"
    echo "| **嵌合体 Reads 比例 (Chimeric Reads)** | \`${CHIMERIC_PAIRS}%\` (Paired) / \`${NON_CHIMERIC_AMBIGUOUS}%\` (Ambiguous) | Hi-C 文库构建过程中的正常副产物，比例越低越好。 |"
    echo "| **含连接基序比例 (Ligation Motif Present)** | \`${LIGATION_MOTIF}%\` | 证明 Hi-C 邻近连接实验成功的直接证据，越高越好。 |"
    echo "| **PCR 重复 Reads 比例 (Total Duplicates)** | \`${DUPS_PERCENT_OF_TOTAL}%\` (\`${TOTAL_DUPS}\` 对) | 文库的 PCR 扩增冗余，越低说明文库复杂度越高。 |"
    echo "| **<span style='color:red'>最终有效交互率 (Total Unique)</span>** | **\`${UNIQUE_PERCENT_OF_TOTAL}%\`** (\`${TOTAL_UNIQUE}\` 对) | **【黄金指标】** 最终用于分析的“精华数据”占总数据的比例。**这是评估文库质量最核心的依据。** |"
    echo ""
} > "${REPORT_FILE}"

echo "INFO: --- 质量评估报告生成完毕 ---"
echo "INFO: 报告文件路径: ${REPORT_FILE}"
echo "INFO: 您可以查看此文件以获取详细的评估结果。"
echo ""
echo "INFO: 作业结束于: $(date)"