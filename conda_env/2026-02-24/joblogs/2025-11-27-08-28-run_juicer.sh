#!/bin/bash
# =============================================================================
#       Juicer "单机版" 分析流程 - 超算提交脚本 (CSUB) - 最终完整修正版 v10
# =============================================================================
#
# 最终修正 v10:
#   - 提供了完整的、无省略的脚本。
#   - 解决了酶切位点文件生成在错误位置并导致目标文件为空的问题。
#   - 脚本现在会先在当前目录生成文件，然后将其排序并移动到正确的 restriction_sites/ 目录。
#
# =============================================================================

# =============================================================================
# --- 计算节点环境与用户配置 ---
# =============================================================================
echo "INFO: 作业开始于: $(date)"
echo "INFO: 运行于计算节点: $(hostname)"

set -e # 若有任何命令执行失败，则立即退出脚本
set -o pipefail # 管道中任何一个命令失败，都视为整个管道失败

# 合并数据
cat F25A040009292_ARAjyxceD/upload/Est-1/*_1.fq.gz > 01.data/raw/Est1_1.fq.gz
cat F25A040009292_ARAjyxceD/upload/Est-1/*_2.fq.gz > 01.data/raw/Est1_2.fq.gz

# 数据过滤
biopytools fastp -i 01.data/raw -o 01.data/clean

# 移动数据
cp 01.data/clean/Est1_1.clean.fq.gz fastq/Est1_R1.fastq.gz
cp 01.data/clean/Est1_2.clean.fq.gz fastq/Est1_R2.fastq.gz

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

echo "INFO: 正在设置精确的 Java 环境..."
JAVA_BIN_PATH="/share/org/YZWL/yzwl_lixg/miniforge3/envs/juicer_v.1.6/bin"
export PATH="${JAVA_BIN_PATH}:${PATH}"

# 2. 用户配置
JUICER_DIR="/share/org/YZWL/yzwl_lixg/software/juicer"
TOP_DIR="/share/org/YZWL/yzwl_lixg/project/22.Est1_hic/03.hic"
GENOME_ID="Est1"
RESTRICTION_ENZYME="MboI"
N_THREADS=88

# =============================================================================
# --- (C) 自动化文件准备 ---
# =============================================================================
echo "INFO: --- 开始自动化文件准备 ---"
cd "${TOP_DIR}"

# --- 定义所有必需文件的绝对路径 ---
GENOME_FA_FILE="${TOP_DIR}/genome.fa"
CHROM_SIZES_FILE="${TOP_DIR}/references/${GENOME_ID}.chrom.sizes"
SITE_FILE="${TOP_DIR}/restriction_sites/${GENOME_ID}_${RESTRICTION_ENZYME}.txt"
# 临时文件的名字，与 generate_site_positions.py 的默认输出匹配
SITE_FILE_TEMP="./${GENOME_ID}_${RESTRICTION_ENZYME}.txt" 

# --- 1. 创建 Juicer 需要的目录 ---
mkdir -p references restriction_sites
echo "INFO: 输入目录结构检查完毕。"

# --- 2. 准备 BWA 索引 ---
if [ ! -f "${GENOME_FA_FILE}.bwt" ]; then
    echo "INFO: BWA 索引未找到, 正在构建..."
    bwa index "${GENOME_FA_FILE}"
    echo "INFO: BWA 索引构建完成。"
else
    echo "INFO: BWA 索引已存在。"
fi

# --- 3. 准备染色体大小文件 ---
if [ ! -f "${CHROM_SIZES_FILE}" ]; then
    echo "INFO: 正在生成染色体大小文件..."
    samtools faidx "${GENOME_FA_FILE}"
    awk -v OFS='\t' '{print $1, $2}' "${GENOME_FA_FILE}.fai" > "${CHROM_SIZES_FILE}"
    echo "INFO: 染色体大小文件已生成。"
else
    echo "INFO: 染色体大小文件已存在。"
fi

# --- 4. 准备酶切位点文件 (关键修正) ---
if [ -f "${SITE_FILE}" ] && [ -s "${SITE_FILE}" ]; then # 检查文件是否存在且不为空
    echo "INFO: 酶切位点文件已存在且不为空, 跳过。"
else
    echo "INFO: 正在生成酶切位点文件..."
    
    # !! 关键修正: 移除重定向，让脚本在当前目录生成文件 !!
    python2 "${JUICER_DIR}/misc/generate_site_positions.py" ${RESTRICTION_ENZYME} ${GENOME_ID} "${GENOME_FA_FILE}"
    
    # 检查临时文件是否生成且不为空
    if [ ! -s "${SITE_FILE_TEMP}" ]; then
        echo "错误: generate_site_positions.py 未能正确生成酶切位点文件！"
        exit 1
    fi

    # !! 关键修正: 将生成的文件排序并移动到正确的位置 !!
    echo "INFO: 将生成的酶切位点文件移动到目标位置..."
    sort -k1,1V -k2,2n "${SITE_FILE_TEMP}" > "${SITE_FILE}"
    rm "${SITE_FILE_TEMP}" # 删除当前目录下的临时文件
    
    echo "INFO: 酶切位点文件已成功生成并移动到: ${SITE_FILE}"
fi
echo "INFO: --- 所有必需的输入文件准备工作完成 ---"

# =============================================================================
# --- (D) 启动 Juicer 主流程 ---
# =============================================================================
echo "INFO: --- 准备启动 Juicer ---"
JUICER_SH="${JUICER_DIR}/scripts/juicer.sh"

# 在启动前，强制清理上一次运行可能留下的临时和输出目录
echo "INFO: 正在清理旧的 'aligned' 和 'splits' 目录..."
rm -rf aligned splits
echo "INFO: 清理完成。"

echo "INFO: 即将以 ${N_THREADS} 个线程运行 Juicer..."

# 运行 Juicer 主脚本
# 这次，我们提供了所有它需要的输入文件
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
# --- (E) 自动生成质量评估报告 ---
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
TOTAL_READS=$(grep "Sequenced Read Pairs" "${STATS_FILE}" | awk '{print $NF}')
UNMAPPED_LINE=$(grep "unmapped" "${STATS_FILE}" | head -n 1)
UNMAPPED_PERCENT=$(echo ${UNMAPPED_LINE} | awk -F'[()]' '{print $2}' | awk -F'%' '{print $1}')
TOTAL_ALIGNED=$(awk -v unmapped=$UNMAPPED_PERCENT 'BEGIN {printf "%.2f", 100 - unmapped}')

LIGATION_MOTIF_LINE=$(grep "Ligation Motif Present" "${STATS_FILE}")
LIGATION_MOTIF_PERCENT=$(echo ${LIGATION_MOTIF_LINE} | awk -F'[()]' '{print $2}' | awk -F'%' '{print $1}')

UNIQUE_READS_LINE=$(grep "Total Unique" "${STATS_FILE}")
TOTAL_UNIQUE=$(echo ${UNIQUE_READS_LINE} | awk '{print $3}')
UNIQUE_PERCENT_OF_TOTAL=$(echo ${UNIQUE_READS_LINE} | awk -F'[,)]' '{print $2}' | sed 's/%//') # 黄金指标

DUPS_READS_LINE=$(grep "Total Duplicates" "${STATS_FILE}")
TOTAL_DUPS=$(echo ${DUPS_READS_LINE} | awk '{print $3}')
DUPS_PERCENT_OF_TOTAL=$(echo ${DUPS_READS_LINE} | awk -F'[,)]' '{print $2}' | sed 's/%//')

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
        echo "> **数据质量“可接受”或“优秀”。** 核心指标——最终有效交互率 (基于总 reads) 达到了 **${UNIQUE_PERCENT_OF_TOTAL}%**，超过了进行大规模测序的推荐阈值 (10%)。可以安排后续测序。"
    else
        echo "> **数据质量“较差”，需谨慎评估。** 核心指标——最终有效交互率 (基于总 reads) 仅为 **${UNIQUE_PERCENT_OF_TOTAL}%**，低于推荐阈值 (10%)。不建议立即进行大规模测序，请与测序公司讨论实验细节。"
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
    echo "| **含连接基序比例 (Ligation Motif Present)** | \`${LIGATION_MOTIF_PERCENT}%\` | 证明 Hi-C 邻近连接实验成功的直接证据，越高越好。 |"
    echo "| **PCR 重复 Reads 比例 (Total Duplicates)** | \`${DUPS_PERCENT_OF_TOTAL}%\` (共 \`${TOTAL_DUPS}\` 对) | 文库的 PCR 扩增冗余，越低说明文库复杂度越高。 |"
    echo "| **<span style='color:red'>最终有效交互率 (Total Unique)</span>** | **\`${UNIQUE_PERCENT_OF_TOTAL}%\`** (共 \`${TOTAL_UNIQUE}\` 对) | **【黄金指标】** 最终用于分析的“精华数据”占总数据的比例。**这是评估文库质量最核心的依据。** |"
    echo ""
} > "${TOP_DIR}/aligned/QC_SUMMARY_REPORT.md"

echo "INFO: --- 质量评估报告生成完毕 ---"
echo "INFO: 报告文件路径: ${TOP_DIR}/aligned/QC_SUMMARY_REPORT.md"
echo "INFO: 您可以查看此文件以获取详细的评估结果。"
echo ""
echo "INFO: 作业结束于: $(date)"