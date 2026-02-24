#!/bin/bash
set -e # 脚本中任何命令执行失败，则立即退出

# =================================================================
#               📝 日志记录函数 (Logging Function) 📝
# =================================================================
log_info() {
    # 格式: [YYYY-MM-DD HH:MM:SS] INFO - Message
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] INFO - $1"
}

# =================================================================
#               🌟 配置部分 (Configuration Section) 🌟
#       只需要在这里修改路径和参数，下游代码会自动更新
# =================================================================

# --- 1. 项目基础路径 ---
PROJECT_BASE="/share/org/YZWL/yzwl_lixg/project"
ANALYSIS_BASE="${PROJECT_BASE}/94.rice_gas"

# --- 2. 输入/输出目录 ---
# 原始 FastQ 数据所在的目录
FASTQ_INPUT_DIR="${ANALYSIS_BASE}/01.data/clean/resequence " # 假设FastQ原始数据在这里

# 本次分析的各个步骤目录
REF_GENOME_DIR="${ANALYSIS_BASE}/01.data/genome" # 参考基因组所在的目录
MAPPING_DIR="${ANALYSIS_BASE}/02.mapping"

# --- 3. 关键文件 ---
# 参考基因组 FastA 文件路径 (请确保此文件已存在于 REF_GENOME_DIR 中)
REF_GENOME_FASTA="${REF_GENOME_DIR}/MSU.fa"

# --- 4. 工具参数 ---
THREADS_MAPPING=64 # 测序比对使用的线程数,默认只有4个GPU，每个GPU搭配16线程，所以这里只用了64个线程
THREADS_FILTER=64  # 变异过滤使用的线程数，GPU所在的节点只有64个线程，最大只能设置为64


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