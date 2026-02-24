#!/bin/bash
# =============================================================================
#       Juicer "单机版" 分析流程 - 超算提交脚本 (CSUB) - 最终版 v14
# =============================================================================
#
# 最终修正 v14 (拨乱反正):
#   - 严格遵循 Juicer 的命名规范，确保所有基因组相关文件名与 GENOME_ID 匹配。
#   - 脚本会自动创建名为 <GENOME_ID>.fasta 的软链接，并基于此链接进行所有操作。
#   - 强制使用 python2 执行 Juicer 的 misc 脚本，解决版本兼容性问题。
#   - 这是基于所有失败经验总结出的、最稳健的最终方案。
#
# =============================================================================

# --- (A) 作业调度系统指令 ---
#CSUB -L /bin/bash
#CSUB -J Juicer_Final_v14
#CSUB -q c01
#CSUB -o Outlog/%J.out
#CSUB -e Outlog/%J.error
#CSUB -n 64
#CSUB -R "span[hosts=1]"
#CSUB -R "rusage[mem=640000]" 

# =============================================================================
# --- (B) 计算节点环境与用户配置 ---
# =============================================================================
echo "INFO: 作业开始于: $(date)"
echo "INFO: 运行于计算节点: $(hostname)"

set -e 
set -o pipefail

# 1. 环境激活 (!! 非常重要 !!)
#    确保您的环境中包含 java, bwa, samtools, 以及一个可用的 python2
#    -----------------------------------------------------------------
#    # >> 方案一: 如果您使用 Conda 环境
#    source /share/org/YZWL/yzwl_lixg/miniforge3/etc/profile.d/conda.sh
#    conda activate your_juicer_env_name # <--- 这个环境里需要有 python2
#
#    # >> 方案二: 如果您使用 module 系统
#    # module load python/2.7
#    # module load java/1.8
#    # module load bwa/0.7.17
#    # module load samtools/1.12
#    -----------------------------------------------------------------


# 2. 用户配置
JUICER_DIR="/share/org/YZWL/yzwl_lixg/software/juicer"
TOP_DIR="/share/org/YZWL/yzwl_lixg/project/22.Est1_hic/02.test/juicer"
GENOME_ID="Est1"
RESTRICTION_ENZYME="MboI"
N_THREADS=64

# =============================================================================
# --- (C) 自动化文件准备 (最严格、最标准的命名规范) ---
# =============================================================================
echo "INFO: --- 开始自动化文件准备 ---"
cd "${TOP_DIR}"

# --- 定义文件路径 ---
# 原始基因组文件
ORIGINAL_GENOME_FA="${TOP_DIR}/genome.fa"
# !! 关键：所有后续操作都将基于这个与GENOME_ID同名的、位于references目录下的文件 !!
GENOME_FA_FILE_STD="${TOP_DIR}/references/${GENOME_ID}.fasta" 
CHROM_SIZES_FILE="${TOP_DIR}/references/${GENOME_ID}.chrom.sizes"
SITE_FILE="${TOP_DIR}/restriction_sites/${GENOME_ID}_${RESTRICTION_ENZYME}.txt"

# --- 1. 创建目录并建立标准化的基因组链接 ---
mkdir -p references restriction_sites
# 创建一个指向原始基因组的、符合 Juicer 命名规范的软链接
ln -sf "${ORIGINAL_GENOME_FA}" "${GENOME_FA_FILE_STD}"
echo "INFO: 已创建标准基因组链接: ${GENOME_FA_FILE_STD}"

# --- 2. 准备 BWA 索引 (基于标准链接名) ---
if [ ! -f "${GENOME_FA_FILE_STD}.bwt" ]; then
    echo "INFO: BWA 索引未找到, 正在为 ${GENOME_FA_FILE_STD} 构建..."
    bwa index "${GENOME_FA_FILE_STD}"
    echo "INFO: BWA 索引构建完成。"
else
    echo "INFO: BWA 索引已存在。"
fi

# --- 3. 准备染色体大小文件 (基于标准链接名) ---
if [ ! -f "${CHROM_SIZES_FILE}" ]; then
    echo "INFO: 正在生成染色体大小文件..."
    samtools faidx "${GENOME_FA_FILE_STD}"
    awk -v OFS='\t' '{print $1, $2}' "${GENOME_FA_FILE_STD}.fai" > "${CHROM_SIZES_FILE}"
    echo "INFO: 染色体大小文件已生成。"
else
    echo "INFO: 染色体大小文件已存在。"
fi

# --- 4. 准备酶切位点文件 (基于标准链接名和GENOME_ID) ---
if [ -f "${SITE_FILE}" ] && [ -s "${SITE_FILE}" ]; then
    echo "INFO: 酶切位点文件已存在且不为空, 跳过。"
else
    echo "INFO: 正在以 Python 2 的方式生成酶切位点文件..."
    # 调用脚本时，第二个参数<genome>现在是'Est1'，脚本会在'references/'下找到'Est1.fasta'
    python2 "${JUICER_DIR}/misc/generate_site_positions.py" ${RESTRICTION_ENZYME} ${GENOME_ID}
    
    # 检查脚本是否在当前目录生成了文件
    if [ ! -s "./${GENOME_ID}_${RESTRICTION_ENZYME}.txt" ]; then
        echo "错误: generate_site_positions.py 未能正确生成酶切位点文件！"
        exit 1
    fi
    # 将生成的文件排序并移动到正确的位置
    sort -k1,1V -k2,2n "./${GENOME_ID}_${RESTRICTION_ENZYME}.txt" > "${SITE_FILE}"
    rm "./${GENOME_ID}_${RESTRICTION_ENZYME}.txt"
    echo "INFO: 酶切位点文件已成功生成: ${SITE_FILE}"
fi
echo "INFO: --- 所有必需的输入文件准备工作完成 ---"

# =============================================================================
# --- (D) 启动 Juicer 主流程 ---
# =============================================================================
echo "INFO: --- 准备启动 Juicer ---"
JUICER_SH="${JUICER_DIR}/scripts/juicer.sh"

rm -rf aligned splits

echo "INFO: 即将以 ${N_THREADS} 个线程运行 Juicer..."

# 运行 Juicer 主脚本
# -z 参数现在指向标准链接名，与 BWA 索引名完全匹配
bash ${JUICER_SH} \
    -g ${GENOME_ID} \
    -d ${TOP_DIR} \
    -s ${RESTRICTION_ENZYME} \
    -p ${CHROM_SIZES_FILE} \
    -y ${SITE_FILE} \
    -z ${GENOME_FA_FILE_STD} \
    -D ${JUICER_DIR} \
    -t ${N_THREADS}

echo "INFO: Juicer 流程执行完毕。"

# =============================================================================
# --- (E) 自动生成质量评估报告 (完整版) ---
# =============================================================================
echo "INFO: --- 开始生成质量评估报告 ---"
STATS_FILE="${TOP_DIR}/aligned/inter_30.txt"
REPORT_FILE="${TOP_DIR}/aligned/QC_SUMMARY_REPORT.md"
if [ ! -f "${STATS_FILE}" ]; then
    echo "错误: Juicer 核心统计文件未找到: ${STATS_FILE}"; exit 1;
fi
TOTAL_READS=$(grep "Sequenced Read Pairs" "${STATS_FILE}" | awk '{print $NF}')
UNMAPPED_LINE=$(grep "unmapped" "${STATS_FILE}" | head -n 1)
UNMAPPED_PERCENT=$(echo ${UNMAPPED_LINE} | awk -F'[()]' '{print $2}' | awk -F'%' '{print $1}')
TOTAL_ALIGNED=$(awk -v unmapped=$UNMAPPED_PERCENT 'BEGIN {printf "%.2f", 100 - unmapped}')
LIGATION_MOTIF_LINE=$(grep "Ligation Motif Present" "${STATS_FILE}")
LIGATION_MOTIF_PERCENT=$(echo ${LIGATION_MOTIF_LINE} | awk -F'[()]' '{print $2}' | awk -F'%' '{print $1}')
UNIQUE_READS_LINE=$(grep "Total Unique" "${STATS_FILE}")
TOTAL_UNIQUE=$(echo ${UNIQUE_READS_LINE} | awk '{print $3}')
UNIQUE_PERCENT_OF_TOTAL=$(echo ${UNIQUE_READS_LINE} | awk -F'[,)]' '{print $2}' | sed 's/%//')
DUPS_READS_LINE=$(grep "Total Duplicates" "${STATS_FILE}")
TOTAL_DUPS=$(echo ${DUPS_READS_LINE} | awk '{print $3}')
DUPS_PERCENT_OF_TOTAL=$(echo ${DUPS_READS_LINE} | awk -F'[,)]' '{print $2}' | sed 's/%//')
{
    echo "# Juicer Hi-C 数据质量评估报告"; echo ""; echo "**报告生成时间:** $(date)"; echo "**分析目录:** \`${TOP_DIR}\`"; echo ""; echo "---"; echo "";
    echo "## 一、核心结论"; echo "";
    if (( $(echo "$UNIQUE_PERCENT_OF_TOTAL >= 10" | bc -l) )); then
        echo "> **数据质量“可接受”或“优秀”。** 核心指标——最终有效交互率 (基于总 reads) 达到了 **${UNIQUE_PERCENT_OF_TOTAL}%**，可以安排后续测序。";
    else
        echo "> **数据质量“较差”，需谨慎评估。** 核心指标——最终有效交互率 (基于总 reads) 仅为 **${UNIQUE_PERCENT_OF_TOTAL}%**，不建议立即进行大规模测序。";
    fi;
    echo ""; echo "---"; echo ""; echo "## 二、关键产出文件"; echo "";
    echo "- **核心统计文件 (MAPQ>=30):** \`aligned/inter_30.txt\`";
    echo "- **Juicebox 可视化文件 (MAPQ>=30):** \`aligned/inter_30.hic\`";
    echo "  - *您可以下载此文件，并使用 Juicebox 桌面软件打开，以交互式地浏览互作热图。*"; echo ""; echo "---"; echo ""; echo "## 三、详细统计数据解读 (基于 inter_30.txt)"; echo "";
    echo "| 指标 (Metric) | 数量 / 比例 | 解读 |"; echo "|:---|:---|:---|"
    echo "| **总测序 Reads 对数 (Sequenced Read Pairs)** | \`${TOTAL_READS}\` | 原始数据总量。 |"
    echo "| **比对率 (Alignment Rate)** | \`${TOTAL_ALIGNED}%\` | Reads 能够成功比对到参考基因组的比例，越高越好。 |"
    echo "| **含连接基序比例 (Ligation Motif Present)** | \`${LIGATION_MOTIF_PERCENT}%\` | 证明 Hi-C 邻近连接实验成功的直接证据，越高越好。 |"
    echo "| **PCR 重复 Reads 比例 (Total Duplicates)** | \`${DUPS_PERCENT_OF_TOTAL}%\` (共 \`${TOTAL_DUPS}\` 对) | 文库的 PCR 扩增冗余，越低说明文库复杂度越高。 |"
    echo "| **<span style='color:red'>最终有效交互率 (Total Unique)</span>** | **\`${UNIQUE_PERCENT_OF_TOTAL}%\`** (共 \`${TOTAL_UNIQUE}\` 对) | **【黄金指标】** 最终用于分析的“精华数据”占总数据的比例。**这是评估文库质量最核心的依据。** |"; echo ""
} > "${TOP_DIR}/aligned/QC_SUMMARY_REPORT.md"

echo "INFO: --- 质量评估报告生成完毕 ---"
echo "INFO: 报告文件路径: ${TOP_DIR}/aligned/QC_SUMMARY_REPORT.md"
echo ""
echo "INFO: 作业结束于: $(date)"