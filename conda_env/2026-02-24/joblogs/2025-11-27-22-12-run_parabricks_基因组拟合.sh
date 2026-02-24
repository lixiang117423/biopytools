#!/bin/bash
set -e # 脚本中任何命令执行失败，则立即退出

# =================================================================
#               📝 日志记录函数 (Logging Function) 📝
# =================================================================
log_info() {
    # 格式: [YYYY-MM-DD HH:MM:SS] INFO - Message
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] INFO - $1"
}

#!/bin/bash

# ================= 配置区域 =================
# 输入文件夹路径
IN_DIR="/share/org/YZWL/yzwl_lixg/database/soybean/genome/clean"

# 输出文件夹路径
OUT_DIR="/share/org/YZWL/yzwl_lixg/project/21.野生大豆群体/01.data/clean/genome_cell"

# 测序读长 (Read Length)，通常为150
READ_LEN=150

# 插入片段大小 (Insert size)，通常为500
INSERT_SIZE=500

# 生成的数据量计算
# 目标 30G 数据: 30 * 10^9 / (2 * 150) = 100,000,000 pairs
NUM_READS=100000000

# 错误率设置 (默认: 0.02，这里设为0，如果需要模拟真实测序错误请去掉-e 0)
# -e 0 表示生成完全纯净的数据(clean)，不含测序错误
ERROR_RATE=0 

# ===========================================

# 检查并创建输出目录
if [ ! -d "$OUT_DIR" ]; then
    mkdir -p "$OUT_DIR"
    echo "创建输出目录: $OUT_DIR"
fi

# 循环处理 S1 到 S36
for i in {1..36}; do
    SAMPLE="S${i}"
    REF_FILE="${IN_DIR}/${SAMPLE}.fa"
    
    OUT_R1="${OUT_DIR}/${SAMPLE}_1.clean.fq"
    OUT_R2="${OUT_DIR}/${SAMPLE}_2.clean.fq"

    # 检查输入参考基因组文件是否存在
    if [ -f "$REF_FILE" ]; then
        echo "========================================"
        echo "正在处理样品: ${SAMPLE}"
        echo "参考基因组: ${REF_FILE}"
        echo "目标数据量: 30G (Reads数: ${NUM_READS})"
        
        # 运行 wgsim
        # -N: reads对数
        # -1: Read1 长度
        # -2: Read2 长度
        # -d: 插入片段平均长度
        # -e: 测序错误率 (这里设为0以获得 clean data)
        # -r: 突变率 (默认0.001)
        # -R: Indel比例 (默认0.15)
        # -X: Indel扩展概率 (默认0.30)
        wgsim -N ${NUM_READS} \
              -1 ${READ_LEN} \
              -2 ${READ_LEN} \
              -d ${INSERT_SIZE} \
              -e ${ERROR_RATE} \
              ${REF_FILE} \
              ${OUT_R1} \
              ${OUT_R2}

        # 检查 wgsim 是否成功生成文件
        if [ -f "$OUT_R1" ]; then
            echo "wgsim 模拟完成，正在压缩文件..."
            
            # 使用 gzip 压缩 (如果有 pigz 建议替换 gzip 以加速)
            gzip -f ${OUT_R1} &
            gzip -f ${OUT_R2} &
            
            # 等待压缩任务完成再处理下一个样品，防止瞬间占用过多CPU/IO
            wait 
            
            echo "样品 ${SAMPLE} 处理完毕."
        else
            echo "错误: wgsim 生成文件失败 - ${SAMPLE}"
        fi
        
    else
        echo "跳过: 未找到文件 ${REF_FILE}"
    fi
done

echo "wgsm所有任务已完成。"

# =================================================================
#               🌟 配置部分 (Configuration Section) 🌟
#       只需要在这里修改路径和参数，下游代码会自动更新
# =================================================================

# --- 1. 项目基础路径 ---
PROJECT_BASE="/share/org/YZWL/yzwl_lixg/project/"
ANALYSIS_BASE="${PROJECT_BASE}/21.野生大豆群体"

# --- 2. 输入/输出目录 ---
# 原始 FastQ 数据所在的目录
FASTQ_INPUT_DIR="${ANALYSIS_BASE}/01.data/clean/genome_cell" # 假设FastQ原始数据在这里

# 本次分析的各个步骤目录
REF_GENOME_DIR="${ANALYSIS_BASE}/01.data/genome" # 参考基因组所在的目录
MAPPING_DIR="${ANALYSIS_BASE}/02.mapping"
GATK_DIR="${ANALYSIS_BASE}/03.gatk_joint"
FILTER_DIR="${ANALYSIS_BASE}/04.filtered_snp_indel"

# --- 3. 关键文件 ---
# 参考基因组 FastA 文件路径 (请确保此文件已存在于 REF_GENOME_DIR 中)
REF_GENOME_FASTA="${REF_GENOME_DIR}/genome.fa"

# --- 4. 工具参数 ---
THREADS_MAPPING=16 # 测序比对使用的线程数,默认只有4个GPU，每个GPU搭配16线程，所以这里只用了64个线程
THREADS_FILTER=16  # 变异过滤使用的线程数，GPU所在的节点只有64个线程，最大只能设置为64


# =================================================================
#                   🚀 分析流程 (Analysis Pipeline) 🚀
#         这部分代码通常不需要修改，除非流程本身发生变化
# =================================================================

log_info "✨ Parabricks + GATK joint Calling的 SNP/INDEL 鉴定流程 ✨"
log_info "============================================================="

# --- 步骤 1: 测序数据比对 (Mapping) ---
log_info "📊 步骤 1: 运行测序数据比对 (biopytools parabricks)..."
# 注意: biopytools parabricks 工具本身就会输出带时间戳的日志，所以这里我们只打印开始和结束信息
biopytools parabricks \
  -i "${FASTQ_INPUT_DIR}" \
  -o "${MAPPING_DIR}" \
  -r "${REF_GENOME_FASTA}" \
  -t "${THREADS_MAPPING}" \
  --no-joint-calling # 根据您的需求，这里保留了不进行联合变异检测的选项
log_info "✅ 测序数据比对完成。"

log_info "====================================================="
log_info "🎉 流程执行成功！所有步骤均已完成！"
log_info "📂 最终过滤后的VCF文件可以在这里找到: ${FILTER_DIR}"
log_info "====================================================="