#!/bin/bash
# =================================================================
#   重测序全基因组变异检测全流程分析脚本
#   Version: 2.0
#   支持: 自动化质控 -> 比对 -> 变异检测 -> 过滤
# =================================================================

# set -euo pipefail  # 严格错误处理：错误退出、未定义变量报错、管道错误

# =================================================================
#               📝 用户配置区域 (User Configuration)
# =================================================================

# 1. 核心输入路径 (必须修改)
PROJECT_BASE="${PROJECT_BASE:-/share/org/YZWL/yzwl_lixg/project/19.大豆疫霉菌/99.测试全自动流程}"
RAW_FASTQ_DIR="${PROJECT_BASE}/01.data/raw"
REF_GENOME_FA="${PROJECT_BASE}/01.data/genome/Phytophthora_sojae_JS2.cds.fa"

# 2. 工具路径
GTX_BIN="${GTX_BIN:-/share/org/YZWL/yzwl_lixg/software/gtx/bin/gtx}"
GTX_CMD_GEN_SCRIPT="${GTX_CMD_GEN_SCRIPT:-${HOME}/software/scripts/51.生成GTX按染色体合并gVCF的脚本.sh}"

# 3. 线程资源配置
THREADS_MAPPING="${THREADS_MAPPING:-88}"
THREADS_GTX="${THREADS_GTX:-88}"
THREADS_FILTER="${THREADS_FILTER:-88}"

# 4. 样本阈值
GATK_THRESHOLD=10      # < 10 使用 GATK
GTX_SINGLE_THRESHOLD=50 # < 50 使用 GTX 单机模式
GTX_WINDOW_SIZE=20000000 # GTX 分块窗口大小 (20Mb)

# 5. 过滤参数
SNP_MIN_DP=5
INDEL_MIN_DP=5

# =================================================================
#               ⚙️ 系统路径规划
# =================================================================
CLEAN_FASTQ_DIR="${PROJECT_BASE}/01.data/clean"
MAPPING_DIR="${PROJECT_BASE}/02.mapping"
GVCF_DIR="${MAPPING_DIR}/vcf"
JOINT_DIR="${PROJECT_BASE}/03.joint_calling"
FILTER_DIR="${PROJECT_BASE}/04.filtered_snp_indel"
SCRIPT_DIR="${PROJECT_BASE}/00.scripts"
LOG_DIR="${PROJECT_BASE}/99.logs"

# 日志配置
LOG_FILE="${LOG_DIR}/pipeline_$(date '+%Y%m%d_%H%M%S').log"
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m'

# =================================================================
#               🛠️ 工具函数库
# =================================================================

# 日志函数（同时输出到终端和文件）
log_info() {
    local msg="[INFO] $(date '+%Y-%m-%d %H:%M:%S') - $1"
    echo -e "${GREEN}${msg}${NC}" | tee -a "${LOG_FILE}"
}

log_warn() {
    local msg="[WARN] $(date '+%Y-%m-%d %H:%M:%S') - $1"
    echo -e "${YELLOW}${msg}${NC}" | tee -a "${LOG_FILE}"
}

log_error() {
    local msg="[ERROR] $(date '+%Y-%m-%d %H:%M:%S') - $1"
    echo -e "${RED}${msg}${NC}" | tee -a "${LOG_FILE}"
}

log_step() {
    local msg="$1"
    echo -e "\n${BLUE}========================================${NC}" | tee -a "${LOG_FILE}"
    echo -e "${BLUE}${msg}${NC}" | tee -a "${LOG_FILE}"
    echo -e "${BLUE}========================================${NC}" | tee -a "${LOG_FILE}"
}

# 检查命令是否存在
check_command() {
    if ! command -v "$1" &> /dev/null; then
        log_error "必需的命令未找到: $1"
        exit 1
    fi
}

# 检查文件是否存在
check_file() {
    if [ ! -f "$1" ]; then
        log_error "文件不存在: $1"
        exit 1
    fi
}

# 检查目录是否为空
check_dir_not_empty() {
    if [ -z "$(ls -A "$1" 2>/dev/null)" ]; then
        log_error "目录为空或不存在: $1"
        exit 1
    fi
}

# 创建目录
safe_mkdir() {
    mkdir -p "$1" || {
        log_error "无法创建目录: $1"
        exit 1
    }
}

# 计算样本数量
count_samples() {
    local dir="$1"
    local pattern="$2"
    find "${dir}" -name "${pattern}" 2>/dev/null | wc -l
}

# =================================================================
#               ✅ 预检查模块
# =================================================================

pre_flight_checks() {
    log_step "🔍 执行预检查 (Pre-flight Checks)"
    
    # 检查必需命令
    log_info "检查必需工具..."
    for cmd in bwa samtools gatk biopytools bcftools tabix python3; do
        check_command "$cmd"
    done
    
    # 检查参考基因组
    log_info "检查参考基因组..."
    check_file "${REF_GENOME_FA}"
    
    # 检查原始数据
    log_info "检查原始数据目录..."
    check_dir_not_empty "${RAW_FASTQ_DIR}"
    
    # 检查 GTX（如果需要）
    if [ -n "${GTX_BIN}" ] && [ "${GTX_BIN}" != "skip" ]; then
        check_file "${GTX_BIN}"
    fi
    
    # 创建所有必需目录
    log_info "创建工作目录..."
    for dir in "${CLEAN_FASTQ_DIR}" "${MAPPING_DIR}" "${JOINT_DIR}" \
               "${FILTER_DIR}" "${SCRIPT_DIR}" "${LOG_DIR}"; do
        safe_mkdir "${dir}"
    done
    
    log_info "✅ 预检查通过"
}

# =================================================================
#               📊 基因组索引模块
# =================================================================

build_genome_index() {
    log_step "📊 Step 0: 构建基因组索引"
    
    # BWA 索引
    if [ ! -f "${REF_GENOME_FA}.bwt" ]; then
        log_info "构建 BWA 索引..."
        bwa index "${REF_GENOME_FA}" 2>&1 | tee -a "${LOG_FILE}"
    else
        log_info "BWA 索引已存在，跳过"
    fi
    
    # SAMtools 索引
    if [ ! -f "${REF_GENOME_FA}.fai" ]; then
        log_info "构建 SAMtools 索引..."
        samtools faidx "${REF_GENOME_FA}" 2>&1 | tee -a "${LOG_FILE}"
    else
        log_info "SAMtools 索引已存在，跳过"
    fi
    
    # GATK 字典
    local ref_dict="${REF_GENOME_FA%.fa}.dict"
    if [ ! -f "${ref_dict}" ]; then
        log_info "构建 GATK 字典..."
        gatk CreateSequenceDictionary \
            -R "${REF_GENOME_FA}" \
            -O "${ref_dict}" 2>&1 | tee -a "${LOG_FILE}"
    else
        log_info "GATK 字典已存在，跳过"
    fi
    
    log_info "✅ 基因组索引准备完成"
}

# =================================================================
#               🧹 质控模块
# =================================================================

run_quality_control() {
    log_step "🧹 Step 1: Fastp 质量控制"
    
    local raw_count=$(count_samples "${RAW_FASTQ_DIR}" "*.fq.gz")
    log_info "检测到 ${raw_count} 个原始 FASTQ 文件"
    
    if [ "${raw_count}" -eq 0 ]; then
        log_error "未找到原始 FASTQ 文件 (*.fq.gz)"
        exit 1
    fi
    
    log_info "开始质控处理..."
    biopytools fastp \
        -i "${RAW_FASTQ_DIR}" \
        -o "${CLEAN_FASTQ_DIR}" \
        --read1-suffix "_1.clean.fq.gz" \
        --read2-suffix "_2.clean.fq.gz" 2>&1 | tee -a "${LOG_FILE}"
    
    local clean_count=$(count_samples "${CLEAN_FASTQ_DIR}" "*.fq.gz")
    log_info "✅ 质控完成，生成 ${clean_count} 个清洁文件"
}

# =================================================================
#               🗺️ 比对模块
# =================================================================

run_mapping() {
    log_step "🗺️ Step 2: 序列比对 (Mapping)"
    
    log_info "使用 ${THREADS_MAPPING} 线程进行比对..."
    biopytools parabricks \
        -i "${CLEAN_FASTQ_DIR}" \
        -o "${MAPPING_DIR}" \
        -r "${REF_GENOME_FA}" \
        -t "${THREADS_MAPPING}" \
        --read1-pattern "*_1.clean.fq.gz" \
        --read2-pattern "*_2.clean.fq.gz" \
        --no-joint-calling 2>&1 | tee -a "${LOG_FILE}"
    
    local gvcf_count=$(count_samples "${GVCF_DIR}" "*.g.vcf.gz")
    log_info "✅ 比对完成，生成 ${gvcf_count} 个 gVCF 文件"
}

# =================================================================
#               🧬 变异检测模块 - GATK 模式
# =================================================================

run_gatk_joint_calling() {
    log_info "👉 使用 GATK GenotypeGVCFs 模式"
    
    biopytools gatk-joint \
        -i "${GVCF_DIR}" \
        -o "${JOINT_DIR}" \
        -r "${REF_GENOME_FA}" 2>&1 | tee -a "${LOG_FILE}"
    
    # 自动识别输出文件
    local raw_vcf="${JOINT_DIR}/joint_genotyping_raw.vcf.gz"
    if [ -f "${JOINT_DIR}/joint_genotyping_merged_filtered.vcf.gz" ]; then
        raw_vcf="${JOINT_DIR}/joint_genotyping_merged_filtered.vcf.gz"
    fi
    
    if [ ! -f "${raw_vcf}" ]; then
        log_error "GATK 未生成预期的 VCF 文件"
        exit 1
    fi
    
    echo "${raw_vcf}"
}

# =================================================================
#               🧬 变异检测模块 - GTX 单机模式
# =================================================================

run_gtx_single_machine() {
    log_info "👉 使用 GTX 单机模式 (样本数: 10-49)"
    
    local output_vcf="${JOINT_DIR}/gtx_joint_raw.vcf.gz"
    local tmp_dir="${JOINT_DIR}/tmp_gtx"
    safe_mkdir "${tmp_dir}"
    
    # 动态构建参数数组
    local gtx_args=()
    gtx_args+=("-r" "${REF_GENOME_FA}")
    gtx_args+=("-o" "${output_vcf}")
    gtx_args+=("-t" "${THREADS_GTX}")
    gtx_args+=("--tmp-dir" "${tmp_dir}")
    
    # 安全读取 gVCF 文件列表
    log_info "收集 gVCF 文件列表..."
    while IFS= read -r gvcf_file; do
        gtx_args+=("-v" "${gvcf_file}")
    done < <(find "${GVCF_DIR}" -name "*.g.vcf.gz" -type f)
    
    log_info "开始 GTX 联合检测..."
    faketime '2020-10-20 00:00:00' "${GTX_BIN}" joint "${gtx_args[@]}" 2>&1 | tee -a "${LOG_FILE}"
    
    if [ ! -f "${output_vcf}" ]; then
        log_error "GTX 未生成预期的 VCF 文件"
        exit 1
    fi
    
    log_info "清理临时文件..."
    rm -rf "${tmp_dir}"
    
    echo "${output_vcf}"
}

# =================================================================
#               🧬 变异检测模块 - GTX 集群模式（生成脚本）
# =================================================================

generate_gtx_cluster_scripts() {
    log_warn "👉 大规模样本模式 (>= 50)，生成集群任务脚本"
    
    local chunks_dir="${JOINT_DIR}/chunks"
    local gtx_job_script="${JOINT_DIR}/01.run_gtx_jobs.sh"
    local merge_py_script="${SCRIPT_DIR}/02.merge_vcf.py"
    local final_vcf="${JOINT_DIR}/merged_all.vcf.gz"
    
    safe_mkdir "${chunks_dir}"
    
    # 1. 检查 GTX 命令生成脚本
    check_file "${GTX_CMD_GEN_SCRIPT}"
    
    # 2. 调用外部脚本生成分块命令
    log_info "⚙️ 生成分块变异检测命令 (窗口: ${GTX_WINDOW_SIZE} bp)..."
    bash "${GTX_CMD_GEN_SCRIPT}" \
        -g "${GTX_BIN}" \
        -r "${REF_GENOME_FA}" \
        -i "${GVCF_DIR}" \
        -o "${chunks_dir}" \
        -w "${GTX_WINDOW_SIZE}" \
        -s "${gtx_job_script}" \
        -t "${THREADS_GTX}" 2>&1 | tee -a "${LOG_FILE}"
    
    chmod +x "${gtx_job_script}"
    
    # 3. 生成 Python 合并脚本
    log_info "📝 生成 VCF 合并脚本..."
    cat << 'PYTHON_SCRIPT' > "${merge_py_script}"
#!/usr/bin/env python3
"""
VCF 文件自然排序合并脚本
用途: 合并按染色体区间分块的 GTX joint calling 结果
"""
import os
import sys
import glob
import re
import subprocess
import tempfile

def natural_sort_key(filename):
    """自然排序：正确处理数字序列（如 chr1, chr2, chr10）"""
    basename = os.path.basename(filename)
    return [int(text) if text.isdigit() else text.lower() 
            for text in re.split(r'([0-9]+)', basename)]

def main():
    if len(sys.argv) != 3:
        print("用法: python3 merge_vcf.py <input_dir> <output_vcf>", file=sys.stderr)
        sys.exit(1)
    
    input_dir = sys.argv[1]
    output_file = sys.argv[2]
    
    # 查找所有 VCF 文件并自然排序
    vcf_files = sorted(
        glob.glob(os.path.join(input_dir, "*.joint.vcf.gz")),
        key=natural_sort_key
    )
    
    if not vcf_files:
        print(f"错误: 在 {input_dir} 中未找到 *.joint.vcf.gz 文件", file=sys.stderr)
        sys.exit(1)
    
    print(f"找到 {len(vcf_files)} 个 VCF 文件，开始合并...")
    
    # 创建文件列表
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as tmp:
        for vcf in vcf_files:
            tmp.write(f"{vcf}\n")
        list_path = tmp.name
    
    try:
        # 使用 bcftools concat 合并
        print("执行合并命令...")
        subprocess.check_call(
            f"bcftools concat -f {list_path} -a -O z -o {output_file} --threads 48",
            shell=True
        )
        
        # 建立索引
        print("建立索引...")
        subprocess.check_call(f"tabix -p vcf {output_file}", shell=True)
        
        print(f"✅ 合并成功: {output_file}")
        
    except subprocess.CalledProcessError as e:
        print(f"❌ 合并失败: {e}", file=sys.stderr)
        sys.exit(1)
    finally:
        os.remove(list_path)

if __name__ == "__main__":
    main()
PYTHON_SCRIPT
    
    chmod +x "${merge_py_script}"
    
    # 4. 打印手动操作指南
    cat << MANUAL_GUIDE

========================================================================
🛑 自动化流程已暂停 (样本数 >= ${GTX_SINGLE_THRESHOLD})
========================================================================

由于样本量较大，为避免长时间占用管理节点，请按以下步骤手动操作：

📋 第一步：投递变异检测任务
-----------------------------------
任务脚本: ${gtx_job_script}

批量投递gVCF文件合并任务:
  batch_sub -i ${gtx_job_script} -j gtx_joint -s 5 -m 800

📋 第二步：等待任务完成并合并结果
-----------------------------------
  # 检查所有分块是否完成
  expected_chunks=\$(wc -l < ${gtx_job_script})
  actual_chunks=\$(ls ${chunks_dir}/*.joint.vcf.gz 2>/dev/null | wc -l)
  echo "预期: \${expected_chunks}, 实际: \${actual_chunks}"
  
  # 执行合并
  python3 ${merge_py_script} \\
      ${chunks_dir} \\
      ${final_vcf}

📋 第三步：运行变异过滤
-----------------------------------
  biopytools filter-snp-indel \\
      -i ${final_vcf} \\
      -o ${FILTER_DIR} \\
      -t ${THREADS_FILTER} \\
      --snp-dp ${SNP_MIN_DP} \\
      --indel-dp ${INDEL_MIN_DP}

========================================================================
💡 提示: 所有操作日志保存在 ${LOG_FILE}
========================================================================

MANUAL_GUIDE
}

# =================================================================
#               🧬 变异检测主控模块
# =================================================================

run_joint_calling() {
    log_step "🧬 Step 3: 联合变异检测 (Joint Calling)"
    
    local sample_count=$(count_samples "${GVCF_DIR}" "*.g.vcf.gz")
    log_info "检测到 ${sample_count} 个 gVCF 样本"
    
    if [ "${sample_count}" -eq 0 ]; then
        log_error "未找到任何 gVCF 文件"
        exit 1
    fi
    
    local raw_vcf=""
    
    # 根据样本数选择策略
    if [ "${sample_count}" -ge "${GTX_SINGLE_THRESHOLD}" ]; then
        # >= 50 样本：生成脚本后退出
        generate_gtx_cluster_scripts
        exit 0  # 正常退出，等待用户手动操作
        
    elif [ "${sample_count}" -ge "${GATK_THRESHOLD}" ]; then
        # 10-49 样本：GTX 单机模式
        raw_vcf=$(run_gtx_single_machine)
        
    else
        # < 10 样本：GATK 模式
        raw_vcf=$(run_gatk_joint_calling)
    fi
    
    log_info "✅ 变异检测完成"
    echo "${raw_vcf}"
}

# =================================================================
#               🧹 变异过滤模块
# =================================================================

run_variant_filtering() {
    local input_vcf="$1"
    
    log_step "🧹 Step 4: 变异过滤 (Filtering)"
    
    check_file "${input_vcf}"
    
    log_info "输入 VCF: ${input_vcf}"
    log_info "过滤参数: SNP DP >= ${SNP_MIN_DP}, InDel DP >= ${INDEL_MIN_DP}"
    
    biopytools filter-snp-indel \
        -i "${input_vcf}" \
        -o "${FILTER_DIR}" \
        -t "${THREADS_FILTER}" \
        --snp-dp "${SNP_MIN_DP}" \
        --indel-dp "${INDEL_MIN_DP}" 2>&1 | tee -a "${LOG_FILE}"
    
    log_info "✅ 过滤完成"
}

# =================================================================
#               🎯 主流程入口
# =================================================================

main() {
    log_step "✨ 全基因组重测序自动化分析流程启动 ✨"
    log_info "项目路径: ${PROJECT_BASE}"
    log_info "日志文件: ${LOG_FILE}"
    
    # 执行流程
    pre_flight_checks
    build_genome_index
    run_quality_control
    run_mapping
    
    # 变异检测（可能在此退出）
    local final_vcf
    final_vcf=$(run_joint_calling)
    
    # 如果没有退出，继续过滤
    run_variant_filtering "${final_vcf}"
    
    # 最终报告
    log_step "🎉 全流程执行成功！"
    log_info "📂 最终结果目录: ${FILTER_DIR}"
    log_info "📊 日志文件: ${LOG_FILE}"
    log_info "====================================================="
}

# =================================================================
#               🚀 脚本执行
# =================================================================

# 捕获中断信号
trap 'log_error "脚本被用户中断"; exit 130' INT TERM

# 启动主流程
main "$@"