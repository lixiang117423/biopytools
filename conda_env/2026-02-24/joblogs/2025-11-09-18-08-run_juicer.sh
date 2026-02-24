#!/bin/bash
# =============================================================================
#       Juicer v1.6 "单机版" 分析流程 - 超算提交脚本 (CSUB) - 终极纯净精确版
# =============================================================================

# =============================================================================
# --- 计算节点环境与用户配置 ---
# =============================================================================
echo "INFO: 作业开始于: $(date)"
set -e 
set -o pipefail

# 1. 环境激活 (!! 关键步骤 !!)
#    我们将把您指定的 Java 路径，以及其他工具所在的 Conda 环境 bin 目录
#    都添加到 PATH 的最顶端。
echo "INFO: 正在设置精确的运行环境..."
CONDA_ENV_BIN_PATH="/share/org/YZWL/yzwl_lixg/miniforge3/envs/juicer_v.1.6/bin"
export PATH="${CONDA_ENV_BIN_PATH}:${PATH}"

# 验证环境
echo "INFO: 当前使用的 Java 版本是:"
java -version
echo "INFO: 当前使用的 bwa 是:"
which bwa
echo "INFO: 当前使用的 python2 是:"
which python2

# 2. 用户配置 (已根据您的最新信息完全更新)
# Juicer v1.6 的主安装目录
JUICER_DIR="/share/org/YZWL/yzwl_lixg/software/juicer"
# 您的纯英文项目工作目录
TOP_DIR="/share/org/YZWL/yzwl_lixg/project/22.Est1_hic/02.test/juicer"
# 基因组 ID
GENOME_ID="Est1"
# 限制性内切酶
RESTRICTION_ENZYME="MboI"
# CPU 线程数
N_THREADS=64

# =============================================================================
# --- (C) 自动化文件准备 ---
# =============================================================================
echo "INFO: --- 开始自动化文件准备 ---"
cd "${TOP_DIR}"

# 定义文件路径 (现在都在纯英文目录下)
GENOME_FA_FILE="./genome.fa"
CHROM_SIZES_FILE="./references/${GENOME_ID}.chrom.sizes"
SITE_FILE="./restriction_sites/${GENOME_ID}_${RESTRICTION_ENZYME}.txt"

# 创建目录
mkdir -p references restriction_sites

# 准备BWA索引
if [ ! -f "${GENOME_FA_FILE}.bwt" ]; then
    bwa index "${GENOME_FA_FILE}"
fi
# 准备染色体大小文件
if [ ! -f "${CHROM_SIZES_FILE}" ]; then
    samtools faidx "${GENOME_FA_FILE}"
    awk -v OFS='\t' '{print $1, $2}' "${GENOME_FA_FILE}.fai" > "${CHROM_SIZES_FILE}"
fi
# 准备酶切位点文件 (使用python2)
if [ -f "${SITE_FILE}" ] && [ -s "${SITE_FILE}" ]; then
    echo "INFO: 酶切位点文件已存在且不为空, 跳过。"
else
    # v1.6 的脚本在 CPU/misc 目录下
    python2 "${JUICER_DIR}/CPU/misc/generate_site_positions.py" ${GENOME_ID} "${GENOME_FA_FILE}" > "${SITE_FILE}"
fi
echo "INFO: --- 所有必需的输入文件准备工作完成 ---"

# =============================================================================
# --- (D) 启动 Juicer v1.6 主流程 ---
# =============================================================================
echo "INFO: --- 准备启动 Juicer v1.6 ---"
JUICER_SH="${JUICER_DIR}/CPU/juicer.sh"

# 清理旧的输出和临时目录
rm -rf aligned splits

# 运行 Juicer 主脚本
bash ${JUICER_SH} \
    -g ${GENOME_ID} \
    -d ${TOP_DIR} \
    -s ${RESTRICTION_ENZYME} \
    -p ${CHROM_SIZES_FILE} \
    -y ${SITE_FILE} \
    -z ${GENOME_FA_FILE} \
    -D ${JUICER_DIR}/CPU \
    -t ${N_THREADS}

echo "INFO: Juicer 流程执行完毕。"

# =============================================================================
# --- (E) 在标准输出日志中打印最终QC报告 ---
# =============================================================================
echo ""
echo "================================================================"
echo "          Juicer Hi-C 数据质量评估报告 (基于 inter_30.txt)        "
echo "================================================================"
echo ""

STATS_FILE="${TOP_DIR}/aligned/inter_30.txt"

if [ ! -f "${STATS_FILE}" ]; then
    echo "错误: Juicer 核心统计文件 'inter_30.txt' 未找到！"
    exit 1
fi

# 直接使用 cat 和 grep 打印格式化的报告
TOTAL_READS=$(grep "Sequenced Read Pairs" "${STATS_FILE}")
UNIQUE_READS=$(grep "Total Unique" "${STATS_FILE}")
DUPS_READS=$(grep "Total Duplicates" "${STATS_FILE}")
UNIQUE_PERCENT_OF_TOTAL=$(echo ${UNIQUE_READS} | awk -F'[,)]' '{print $2}' | sed 's/%//')

echo "## 一、核心结论"
echo ""
if (( $(echo "$UNIQUE_PERCENT_OF_TOTAL >= 10" | bc -l) )); then
    echo "> **数据质量“可接受”或“优秀”。** 核心指标——最终有效交互率达到了 **${UNIQUE_PERCENT_OF_TOTAL}%**，可以安排后续测序。"
else
    echo "> **数据质量“较差”，需谨慎评估。** 核心指标——最终有效交互率仅为 **${UNIQUE_PERCENT_OF_TOTAL}%**，不建议立即进行大规模测序。"
fi
echo ""
echo "---"
echo ""
echo "## 二、详细统计数据"
echo ""
grep "Sequenced Read Pairs" "${STATS_FILE}"
grep "Unmapped" "${STATS_FILE}"
grep "Ligation Motif Present" "${STATS_FILE}"
grep "Total Duplicates" "${STATS_FILE}"
echo ""
echo "**黄金指标 -->**"
grep "Total Unique" "${STATS_FILE}"
echo ""
echo "---"
echo ""

# 检查 .hic 文件是否生成
if [ -f "${TOP_DIR}/aligned/inter_30.hic" ]; then
    echo "## 三、产出文件"
    echo "- **成功:** \`.hic\` 可视化文件已生成于 \`aligned/inter_30.hic\`"
else
    echo "## 三、产出文件"
    echo "- **警告:** \`.hic\` 文件未能生成。对于小数据量这是正常现象，评估结果不受影响。"
fi
echo "================================================================"
echo ""

echo "INFO: 作业结束于: $(date)"