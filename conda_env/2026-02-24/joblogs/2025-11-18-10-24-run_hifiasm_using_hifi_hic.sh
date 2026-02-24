#!/bin/bash
#########################################################################
# 🧬 HiFi + Hi-C 基因组组装脚本
# 软件: hifiasm
# 样本: OV53_1
# 预估基因组大小: 1.4-1.5G
#########################################################################

set -e  # 遇到错误立即退出
set -u  # 使用未定义变量时报错

# 合并文件
zcat /share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/01.data/hic/raw/ov53-1-HIC1/*_1.fq.gz >> /share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/01.data/hic/OV53_1-hic_R1.fastq
zcat /share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/01.data/hic/raw/ov53-1-HIC1/*_2.fq.gz >> /share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/01.data/hic/OV53_1-hic_R2.fastq

gzip /share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/01.data/hic/OV53_1-hic_R1.fastq
gzip /share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/01.data/hic/OV53_1-hic_R2.fastq


# 📁 定义路径变量
HIFI_DATA="/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/01.data/hifi/OV53_1-hifi.fq"
HIC_R1="/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/01.data/hic/OV53_1-hic_R1.fastq.gz"
HIC_R2="/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/01.data/hic/OV53_1-hic_R2.fastq.gz"
BASE_DIR="/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/02.assembly"

# ⚙️ 定义参数
PREFIX="OV53_1"
THREADS=88
GENOME_SIZE="1.45g"  # 取中间值1.45G
N_HAP=2              # 二倍体

# 📂 定义目录结构
WORK_DIR="${BASE_DIR}/${PREFIX}"
RAW_DIR="${WORK_DIR}/01.raw_output"      # 原始输出文件
FASTA_DIR="${WORK_DIR}/02.fasta"         # FASTA格式文件
LOG_DIR="${WORK_DIR}/03.logs"            # 日志文件
STAT_DIR="${WORK_DIR}/04.statistics"     # 统计信息

# 📋 打印运行信息
echo "=========================================="
echo "🧬 开始基因组组装"
echo "=========================================="
echo "📌 样本名称: ${PREFIX}"
echo "📌 基因组大小: ${GENOME_SIZE}"
echo "📌 线程数: ${THREADS}"
echo "📌 倍性: ${N_HAP}"
echo "📌 HiFi数据: ${HIFI_DATA}"
echo "📌 Hi-C R1: ${HIC_R1}"
echo "📌 Hi-C R2: ${HIC_R2}"
echo "📌 工作目录: ${WORK_DIR}"
echo "=========================================="
echo ""

# 🔍 检查输入文件
echo "🔍 检查输入文件..."
if [ ! -f "${HIFI_DATA}" ]; then
    echo "❌ 错误: HiFi文件不存在: ${HIFI_DATA}"
    exit 1
fi

if [ ! -f "${HIC_R1}" ]; then
    echo "❌ 错误: Hi-C R1文件不存在: ${HIC_R1}"
    exit 1
fi

if [ ! -f "${HIC_R2}" ]; then
    echo "❌ 错误: Hi-C R2文件不存在: ${HIC_R2}"
    exit 1
fi
echo "✅ 所有输入文件检查通过"
echo ""

# 📁 创建目录结构
echo "📁 创建目录结构..."
mkdir -p ${RAW_DIR}
mkdir -p ${FASTA_DIR}
mkdir -p ${LOG_DIR}
mkdir -p ${STAT_DIR}
echo "✅ 目录结构已创建:"
echo "   📂 ${WORK_DIR}/"
echo "      ├── 📂 01.raw_output/     (原始GFA等文件)"
echo "      ├── 📂 02.fasta/          (FASTA格式文件)"
echo "      ├── 📂 03.logs/           (日志文件)"
echo "      └── 📂 04.statistics/     (统计信息)"
echo ""

# ⏰ 记录开始时间
START_TIME=$(date +%s)
echo "⏰ 开始时间: $(date '+%Y-%m-%d %H:%M:%S')"
echo ""

# 🚀 运行hifiasm组装
echo "=========================================="
echo "🚀 运行hifiasm组装 (HiFi + Hi-C模式)"
echo "=========================================="

cd ${RAW_DIR}

hifiasm \
    -o ${PREFIX} \
    -t ${THREADS} \
    --hg-size ${GENOME_SIZE} \
    --n-hap ${N_HAP} \
    --h1 ${HIC_R1} \
    --h2 ${HIC_R2} \
    --primary \
    -l 3 \
    ${HIFI_DATA} 2>&1 | tee ${LOG_DIR}/hifiasm_assembly.log

# ✅ 检查组装是否成功
if [ $? -eq 0 ]; then
    echo ""
    echo "=========================================="
    echo "✅ 组装完成！"
    echo "=========================================="
else
    echo ""
    echo "=========================================="
    echo "❌ 组装失败，请检查错误信息"
    echo "=========================================="
    exit 1
fi

# 🔄 转换GFA为FASTA格式
echo ""
echo "=========================================="
echo "🔄 转换GFA为FASTA格式"
echo "=========================================="

cd ${RAW_DIR}

# 定义需要转换的GFA文件和对应的FASTA文件名
declare -A GFA_FILES=(
    ["${PREFIX}.bp.hap1.p_ctg.gfa"]="${FASTA_DIR}/${PREFIX}.hap1.primary.fa"
    ["${PREFIX}.bp.hap2.p_ctg.gfa"]="${FASTA_DIR}/${PREFIX}.hap2.primary.fa"
    ["${PREFIX}.bp.p_ctg.gfa"]="${FASTA_DIR}/${PREFIX}.primary.fa"
    ["${PREFIX}.bp.a_ctg.gfa"]="${FASTA_DIR}/${PREFIX}.alternate.fa"
)

# 转换函数
convert_gfa_to_fasta() {
    local gfa_file=$1
    local fasta_file=$2
    
    if [ -f "${gfa_file}" ]; then
        echo "🔹 转换: $(basename ${gfa_file}) -> $(basename ${fasta_file})"
        awk '/^S/{print ">"$2; print $3}' ${gfa_file} > ${fasta_file}
        
        # 统计序列数和总长度
        local num_seqs=$(grep -c "^>" ${fasta_file})
        local total_len=$(awk '/^>/{next}{sum+=length($0)}END{print sum}' ${fasta_file})
        echo "   ✓ 序列数: ${num_seqs}"
        echo "   ✓ 总长度: ${total_len} bp"
        echo ""
    else
        echo "⚠️  警告: 文件不存在 - ${gfa_file}"
        echo ""
    fi
}

# 执行转换
for gfa_file in "${!GFA_FILES[@]}"; do
    convert_gfa_to_fasta "${gfa_file}" "${GFA_FILES[$gfa_file]}"
done

echo "✅ FASTA文件转换完成"
echo ""

# 📊 生成统计信息
echo "=========================================="
echo "📊 生成组装统计信息"
echo "=========================================="

STAT_FILE="${STAT_DIR}/${PREFIX}_assembly_statistics.txt"

{
    echo "=========================================="
    echo "基因组组装统计报告"
    echo "样本: ${PREFIX}"
    echo "日期: $(date '+%Y-%m-%d %H:%M:%S')"
    echo "=========================================="
    echo ""
    
    for fasta_file in "${GFA_FILES[@]}"; do
        if [ -f "${fasta_file}" ]; then
            echo "📄 文件: $(basename ${fasta_file})"
            echo "----------------------------------------"
            
            # 序列数量
            num_seqs=$(grep -c "^>" ${fasta_file})
            echo "序列数量: ${num_seqs}"
            
            # 总长度
            total_len=$(awk '/^>/{next}{sum+=length($0)}END{print sum}' ${fasta_file})
            echo "总长度: ${total_len} bp ($(echo "scale=2; ${total_len}/1000000000" | bc) Gb)"
            
            # 最长序列
            max_len=$(awk '/^>/{if(seq){print length(seq)}; seq=""; next}{seq=seq$0}END{if(seq){print length(seq)}}' ${fasta_file} | sort -rn | head -1)
            echo "最长序列: ${max_len} bp"
            
            # 计算N50
            awk '/^>/{if(seq){print length(seq)}; seq=""; next}{seq=seq$0}END{if(seq){print length(seq)}}' ${fasta_file} | \
            sort -rn | \
            awk '{len[NR]=$1; sum+=$1} END {
                target=sum/2; 
                cumsum=0; 
                for(i=1; i<=NR; i++){
                    cumsum+=len[i]; 
                    if(cumsum>=target){
                        print "N50: "len[i]" bp"; 
                        break
                    }
                }
            }'
            
            echo ""
        fi
    done
    
    echo "=========================================="
    echo "文件位置"
    echo "=========================================="
    echo "原始输出: ${RAW_DIR}"
    echo "FASTA文件: ${FASTA_DIR}"
    echo "日志文件: ${LOG_DIR}"
    echo "统计文件: ${STAT_DIR}"
    echo ""
    
} | tee ${STAT_FILE}

echo "✅ 统计信息已保存至: ${STAT_FILE}"
echo ""

# ⏰ 记录结束时间
END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))
HOURS=$((ELAPSED / 3600))
MINUTES=$(((ELAPSED % 3600) / 60))
SECONDS=$((ELAPSED % 60))

echo ""
echo "=========================================="
echo "⏰ 运行时间统计"
echo "=========================================="
echo "开始时间: $(date -d @${START_TIME} '+%Y-%m-%d %H:%M:%S')"
echo "结束时间: $(date -d @${END_TIME} '+%Y-%m-%d %H:%M:%S')"
echo "总耗时: ${HOURS}小时 ${MINUTES}分钟 ${SECONDS}秒"
echo ""

# 📊 最终输出文件清单
echo "=========================================="
echo "📊 输出文件清单"
echo "=========================================="
echo ""
echo "📂 ${WORK_DIR}/"
echo "   │"
echo "   ├── 📂 01.raw_output/     (原始输出文件)"
echo "   │   ├── 🔸 ${PREFIX}.bp.hap1.p_ctg.gfa"
echo "   │   ├── 🔸 ${PREFIX}.bp.hap2.p_ctg.gfa"
echo "   │   ├── 🔸 ${PREFIX}.bp.p_ctg.gfa"
echo "   │   ├── 🔸 ${PREFIX}.bp.a_ctg.gfa"
echo "   │   └── 🔸 其他中间文件 (.bin, .ovlp等)"
echo "   │"
echo "   ├── 📂 02.fasta/          (FASTA格式)"
echo "   │   ├── 🧬 ${PREFIX}.hap1.primary.fa      (单倍型1主要序列)"
echo "   │   ├── 🧬 ${PREFIX}.hap2.primary.fa      (单倍型2主要序列)"
echo "   │   ├── 🧬 ${PREFIX}.primary.fa           (主要组装)"
echo "   │   └── 🧬 ${PREFIX}.alternate.fa         (备选组装)"
echo "   │"
echo "   ├── 📂 03.logs/           (日志文件)"
echo "   │   └── 📄 hifiasm_assembly.log"
echo "   │"
echo "   └── 📂 04.statistics/     (统计信息)"
echo "       └── 📄 ${PREFIX}_assembly_statistics.txt"
echo ""
echo "=========================================="
echo "🎉 所有任务完成！"
echo "=========================================="
echo ""
echo "💡 下一步建议:"
echo "   1. 查看统计报告: cat ${STAT_FILE}"
echo "   2. 检查组装质量: 使用QUAST或BUSCO"
echo "   3. 可视化组装: 使用Bandage查看GFA文件"
echo ""