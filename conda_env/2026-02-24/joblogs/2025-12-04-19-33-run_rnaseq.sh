#!/bin/bash

# ================= 路径与参数配置 =================

# 1. 输入文件路径
GENOME_FA="/share/org/YZWL/yzwl_lixg/project/94.rice_gas/18.根部转录组/01.data/genome/MSU.fa"
GTF_FILE="/share/org/YZWL/yzwl_lixg/project/94.rice_gas/18.根部转录组/01.data/genome/MSU.gtf"
FASTQ_DIR="/share/org/YZWL/yzwl_lixg/project/94.rice_gas/18.根部转录组/01.data/clean"

# 2. 输出文件夹 (指定路径)
# 所有结果将保存在此目录下
BASE_OUT_DIR="/share/org/YZWL/yzwl_lixg/project/94.rice_gas/18.根部转录组/02.mapping"

# 3. 线程设置 (总计 88)
TOTAL_THREADS=88
# 为了防止管道阻塞，分配给比对和排序不同的线程
# HISAT2 负责计算密集型任务，分配 80 线程
HISAT2_THREADS=80
# Samtools Sort 负责IO密集型任务，分配 8 线程
SAMTOOLS_THREADS=8

# ================= 初始化环境 =================

# 创建输出目录结构
ALIGN_DIR="${BASE_OUT_DIR}/alignment"
EXPR_DIR="${BASE_OUT_DIR}/expression"
LOG_DIR="${BASE_OUT_DIR}/logs"

mkdir -p "${ALIGN_DIR}"
mkdir -p "${EXPR_DIR}"
mkdir -p "${LOG_DIR}"

echo "========================================================"
echo "项目启动时间: $(date)"
echo "输出目录: ${BASE_OUT_DIR}"
echo "输入 FASTQ: ${FASTQ_DIR}"
echo "线程数: ${TOTAL_THREADS}"
echo "========================================================"

# ================= 1. HISAT2 索引检查与构建 =================

# 假设索引文件前缀与 fasta 文件名一致
INDEX_BASE="${GENOME_FA%.*}"

# 检查 .1.ht2 文件是否存在 (HISAT2 索引的标志)
if [ ! -f "${INDEX_BASE}.1.ht2" ]; then
    echo "[$(date)] 未检测到 HISAT2 索引，开始在基因组目录构建索引..."
    echo "注意：请确保对基因组目录有写入权限，否则请手动更改索引输出路径。"
    
    hisat2-build -p ${TOTAL_THREADS} "${GENOME_FA}" "${INDEX_BASE}" > "${LOG_DIR}/index_build.log" 2>&1
    
    if [ $? -ne 0 ]; then
        echo "Error: 索引构建失败！请检查日志: ${LOG_DIR}/index_build.log"
        exit 1
    fi
    echo "[$(date)] 索引构建完成。"
else
    echo "[$(date)] 检测到已有索引: ${INDEX_BASE}，跳过构建。"
fi

# ================= 2. 批量处理样品 =================

# 遍历 R1 文件
for R1 in ${FASTQ_DIR}/*_1.clean.fq.gz
do
    # 如果没有找到文件（目录为空），则跳出
    [ -e "$R1" ] || break

    # 获取样品名 (去除路径和后缀)
    # 例如: /path/to/SRR12079365_1.clean.fq.gz -> SRR12079365
    SAMPLE=$(basename "$R1" _1.clean.fq.gz)
    
    # 构建 R2 文件名
    R2="${FASTQ_DIR}/${SAMPLE}_2.clean.fq.gz"
    
    echo ""
    echo ">>> 开始处理样品: ${SAMPLE} <<<"
    echo "    R1: $R1"
    
    # 检查 R2 是否存在
    if [ ! -f "$R2" ]; then
        echo "    Warning: 未找到对应的 R2 文件 ($R2)，跳过该样品。"
        continue
    fi

    # 定义输出文件路径
    BAM_FILE="${ALIGN_DIR}/${SAMPLE}.sorted.bam"
    GTF_OUT="${EXPR_DIR}/${SAMPLE}/${SAMPLE}.gtf"
    ABUND_OUT="${EXPR_DIR}/${SAMPLE}/${SAMPLE}.gene_abund.tab"
    SAMPLE_LOG_DIR="${LOG_DIR}/${SAMPLE}"
    
    mkdir -p "${SAMPLE_LOG_DIR}"
    mkdir -p "${EXPR_DIR}/${SAMPLE}"

    # --- Step 1: Hisat2 比对 & Samtools 排序 ---
    echo "[$(date)] [${SAMPLE}] Step 1: HISAT2 Alignment..."
    
    # --dta: 下游对接 StringTie 必须参数
    # --summary-file: 输出比对统计信息
    (hisat2 -p ${HISAT2_THREADS} --dta -x "${INDEX_BASE}" \
        -1 "${R1}" -2 "${R2}" \
        --summary-file "${SAMPLE_LOG_DIR}/hisat2_summary.txt" \
        2> "${SAMPLE_LOG_DIR}/hisat2_stderr.log" \
        | samtools sort -@ ${SAMTOOLS_THREADS} -o "${BAM_FILE}" - ) 2> "${SAMPLE_LOG_DIR}/samtools_sort.log"

    # 检查比对步骤是否成功
    if [ $? -ne 0 ]; then
        echo "    Error: [${SAMPLE}] 比对或排序失败！详情请查看 ${SAMPLE_LOG_DIR} 下的日志。"
        echo "    -> 跳过此样品，继续下一个..."
        # 清理可能产生的损坏文件
        rm -f "${BAM_FILE}"
        continue
    fi
    
    # 建立索引（方便后续IGV查看，可选但推荐）
    samtools index -@ ${SAMTOOLS_THREADS} "${BAM_FILE}"

    # --- Step 2: StringTie 定量 ---
    echo "[$(date)] [${SAMPLE}] Step 2: StringTie Quantification..."
    
    # -e: 仅限制在参考注释 GTF 中定量 (不预测新转录本)
    # -B: 生成 Ballgown 输入文件 (可选)
    stringtie -p ${TOTAL_THREADS} -e -B \
        -G "${GTF_FILE}" \
        -o "${GTF_OUT}" \
        -A "${ABUND_OUT}" \
        "${BAM_FILE}" > "${SAMPLE_LOG_DIR}/stringtie.log" 2>&1

    if [ $? -ne 0 ]; then
        echo "    Error: [${SAMPLE}] StringTie 定量失败！详情请查看 ${SAMPLE_LOG_DIR}/stringtie.log"
        echo "    -> 跳过此样品，继续下一个..."
        continue
    fi

    echo "[$(date)] [${SAMPLE}] 处理成功完成。"

done

echo "========================================================"
echo "[$(date)] 所有任务流程结束。"
echo "结果目录: ${BASE_OUT_DIR}"
echo "========================================================"
