#!/bin/bash
# =============================================================================
#  🧬 YaHS 高速染色体挂载流程 - 优化增强版 (v3.0)
#  作者: 基于用户v2.0版本改进
#  日期: 2025-11-21
#  主要改进:
#    - 增强错误处理和日志系统
#    - 添加断点续跑功能
#    - 优化内存和磁盘使用
#    - 自动资源检测和参数调整
#    - 改进统计信息输出
# =============================================================================

# =============================================================================
# --- 🔧 基础配置 ---
# =============================================================================
set -e 
set -o pipefail

# 颜色输出定义
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# 日志函数
log_info() { echo -e "${BLUE}ℹ️  INFO:${NC} $1" | tee -a pipeline.log; }
log_success() { echo -e "${GREEN}✅ SUCCESS:${NC} $1" | tee -a pipeline.log; }
log_warning() { echo -e "${YELLOW}⚠️  WARNING:${NC} $1" | tee -a pipeline.log; }
log_error() { echo -e "${RED}❌ ERROR:${NC} $1" | tee -a pipeline.log; }

# 错误处理
trap 'log_error "脚本在第 $LINENO 行失败，退出码: $?"; exit 1' ERR

# =============================================================================
# --- 💻 环境与参数配置 (用户修改区) ---
# =============================================================================
log_info "作业开始于: $(date '+%Y-%m-%d %H:%M:%S')"
log_info "运行于计算节点: $(hostname)"
log_info "当前用户: $(whoami)"

# 1. 📂 路径设置
WORK_DIR="/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/04.YaHS"
REF_FA="/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/04.YaHS/OV53_1.primary.fa"
R1_FQ="/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/04.YaHS/fastq/OV53_1-hic_R1.fastq.gz"
R2_FQ="/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/04.YaHS/fastq/OV53_1-hic_R2.fastq.gz"

# 2. ⚙️ 软件路径 & 参数
JUICER_JAR="/share/org/YZWL/yzwl_lixg/software/juicer/scripts/common/juicer_tools.jar"
# 【新增】指定 YaHS 自带的 juicer 工具绝对路径
YAHS_JUICER_TOOL="/share/org/YZWL/yzwl_lixg/miniforge3/envs/yahs_v.1.2.2/bin/juicer" 
ENZYME_SEQ="GATC"  # MboI / DpnII / Arima (根据实验设计修改)
THREADS=88
MIN_LEN=10000
MIN_MAPQ=30        # 最小比对质量

# 3. 🔧 资源配置 (自动检测)
TOTAL_MEM=$(free -g | awk '/^Mem:/{print $2}')
SORT_MEM="${TOTAL_MEM}G"
JAVA_MEM="${TOTAL_MEM}G"

log_info "检测到系统总内存: ${TOTAL_MEM}G"
log_info "排序内存配置: ${SORT_MEM}"
log_info "Java堆内存配置: ${JAVA_MEM}"

# 4. 环境激活
export PATH="/share/org/YZWL/yzwl_lixg/miniforge3/envs/yahs_v.1.2.2/bin:$PATH"

# =============================================================================
# --- 🔍 环境检查 ---
# =============================================================================
log_info "执行环境检查..."

# 检查必需软件
REQUIRED_TOOLS=("bwa" "samtools" "yahs" "java" "awk" "sort")
for tool in "${REQUIRED_TOOLS[@]}"; do
    if ! command -v $tool &> /dev/null; then
        log_error "未找到必需工具: $tool"
        exit 1
    fi
    # 检查 YaHS 转换工具
    if [ ! -x "${YAHS_JUICER_TOOL}" ]; then
        log_error "YaHS自带juicer工具未找到或不可执行: ${YAHS_JUICER_TOOL}"
        exit 1
    fi
done
log_success "所有必需软件检查通过"

# 检查输入文件
log_info "检查输入文件完整性..."
for file in "${REF_FA}" "${R1_FQ}" "${R2_FQ}"; do
    if [ ! -f "${file}" ]; then
        log_error "文件不存在: ${file}"
        exit 1
    fi
    # 检查文件是否可读且非空
    if [ ! -r "${file}" ] || [ ! -s "${file}" ]; then
        log_error "文件不可读或为空: ${file}"
        exit 1
    fi
done
log_success "输入文件检查通过"

# 检查Juicer工具
if [ ! -f "${JUICER_JAR}" ]; then
    log_error "Juicer JAR文件未找到: ${JUICER_JAR}"
    exit 1
fi

# 检查磁盘空间 (至少需要100GB可用空间)
AVAILABLE_SPACE=$(df -BG "${WORK_DIR}" | awk 'NR==2 {print $4}' | sed 's/G//')
if [ "${AVAILABLE_SPACE}" -lt 100 ]; then
    log_warning "可用磁盘空间不足100GB (当前: ${AVAILABLE_SPACE}GB)，可能导致空间不足"
fi

# =============================================================================
# --- 📁 工作目录设置 ---
# =============================================================================
mkdir -p "${WORK_DIR}"
cd "${WORK_DIR}" || exit 1

# 创建子目录
mkdir -p logs tmp_files results

# 重定向所有输出到日志
exec > >(tee -a logs/pipeline_$(date +%Y%m%d_%H%M%S).log)
exec 2>&1

log_info "工作目录: ${WORK_DIR}"

# =============================================================================
# --- 步骤 1: 建立索引 (Indexing) ---
# =============================================================================
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
log_info "步骤 1: 检查/构建基因组索引"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

# BWA索引
if [ ! -f "${REF_FA}.bwt" ]; then
    log_info "构建BWA索引..."
    bwa index "${REF_FA}" 2>&1 | tee logs/bwa_index.log
    log_success "BWA索引构建完成"
else
    log_info "发现已有BWA索引，跳过"
fi

# SAMtools索引
if [ ! -f "${REF_FA}.fai" ]; then
    log_info "构建SAMtools索引..."
    samtools faidx "${REF_FA}" 2>&1 | tee logs/samtools_faidx.log
    log_success "SAMtools索引构建完成"
else
    log_info "发现已有SAMtools索引，跳过"
fi

# 验证索引完整性
for ext in amb ann bwt pac sa fai; do
    if [ ! -f "${REF_FA}.${ext}" ]; then
        log_error "索引文件缺失: ${REF_FA}.${ext}"
        exit 1
    fi
done
log_success "所有索引文件完整"

# =============================================================================
# --- 步骤 2: Hi-C 比对与处理 (Mapping & Processing) ---
# =============================================================================
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
log_info "步骤 2: Hi-C测序数据比对与预处理"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

FINAL_BAM="results/aligned_sorted_dedup.bam"

if [ -f "${FINAL_BAM}" ] && [ -f "${FINAL_BAM}.bai" ]; then
    log_info "发现已有BAM文件及索引，跳过比对步骤"
    
    # 显示现有BAM统计信息
    log_info "现有BAM文件统计:"
    samtools flagstat -@ ${THREADS} "${FINAL_BAM}" | tee logs/existing_bam_stats.txt
else
    log_info "开始BWA比对流程..."
    
    # 创建临时目录
    mkdir -p tmp_files/nsort tmp_files/sort
    
    # 完整的比对流程
    log_info "执行BWA MEM比对 (使用 -5SP 参数处理Hi-C数据)..."
    
    bwa mem -5SP -t ${THREADS} "${REF_FA}" "${R1_FQ}" "${R2_FQ}" 2> logs/bwa_mem.log | \
    samtools view -@ ${THREADS} -bS - > tmp_files/aligned.bam
    
    log_info "按read name排序..."
    samtools sort -n -@ ${THREADS} -m 4G \
        -T tmp_files/nsort/split \
        -o tmp_files/aligned_nsorted.bam \
        tmp_files/aligned.bam 2> logs/samtools_nsort.log
    
    log_info "修复mate pair信息..."
    samtools fixmate -m -@ ${THREADS} \
        tmp_files/aligned_nsorted.bam \
        tmp_files/aligned_fixmate.bam 2> logs/samtools_fixmate.log
    
    log_info "按坐标排序..."
    samtools sort -@ ${THREADS} -m 4G \
        -T tmp_files/sort/split \
        -o tmp_files/aligned_sorted.bam \
        tmp_files/aligned_fixmate.bam 2> logs/samtools_sort.log
    
    log_info "标记并移除PCR重复..."
    samtools markdup -r -@ ${THREADS} \
        tmp_files/aligned_sorted.bam \
        "${FINAL_BAM}" 2> logs/samtools_markdup.log
    
    log_info "构建BAM索引..."
    samtools index -@ ${THREADS} "${FINAL_BAM}"
    
    # 清理中间文件
    log_info "清理临时文件..."
    rm -rf tmp_files/nsort tmp_files/sort
    rm -f tmp_files/aligned.bam
    rm -f tmp_files/aligned_nsorted.bam tmp_files/aligned_fixmate.bam
    rm -f tmp_files/aligned_sorted.bam
    
    log_success "比对流程完成"
fi

# 生成详细的比对统计
log_info "生成比对统计信息..."
samtools flagstat -@ ${THREADS} "${FINAL_BAM}" > results/alignment_stats.txt
samtools stats -@ ${THREADS} "${FINAL_BAM}" > results/alignment_detailed_stats.txt

# 输出关键统计
echo ""
log_info "比对质量总结:"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
cat results/alignment_stats.txt
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

# 提取关键指标
TOTAL_READS=$(grep "in total" results/alignment_stats.txt | awk '{print $1}')
MAPPED_READS=$(grep "mapped (" results/alignment_stats.txt | head -1 | awk '{print $1}')
MAPPING_RATE=$(grep "mapped (" results/alignment_stats.txt | head -1 | awk '{print $5}')
PROPERLY_PAIRED=$(grep "properly paired" results/alignment_stats.txt | awk '{print $1}')

log_info "总reads数: ${TOTAL_READS}"
log_info "比对成功reads: ${MAPPED_READS} (${MAPPING_RATE})"
log_info "正确配对reads: ${PROPERLY_PAIRED}"

# 警告检查 - 使用简单的字符串比较
MAPPING_PCT=$(echo "$MAPPING_RATE" | grep -oE '[0-9]+\.[0-9]+' | head -1)
if [ -n "$MAPPING_PCT" ]; then
    # 使用bc进行浮点数比较
    IS_LOW=$(echo "$MAPPING_PCT < 70" | bc -l 2>/dev/null || echo "0")
    if [ "$IS_LOW" = "1" ]; then
        log_warning "比对率偏低 (${MAPPING_RATE})，建议检查数据质量"
    fi
fi

# =============================================================================
# --- 步骤 3: 运行 YaHS 组装 (Scaffolding) ---
# =============================================================================
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
log_info "步骤 3: 执行YaHS染色体挂载"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

OUT_PREFIX="results/yahs_out"

if [ -f "${OUT_PREFIX}_scaffolds_final.fa" ]; then
    log_info "发现已有scaffold文件，跳过YaHS步骤"
else
    log_info "运行YaHS (酶切位点: ${ENZYME_SEQ}, MAPQ阈值: ${MIN_MAPQ})..."
    
    yahs -e "${ENZYME_SEQ}" \
         -q ${MIN_MAPQ} \
         -o "${OUT_PREFIX}" \
         -l ${MIN_LEN} \
         "${REF_FA}" \
         "${FINAL_BAM}" --no-contig-ec  2>&1 | tee logs/yahs.log
    
    # 验证输出
    if [ ! -f "${OUT_PREFIX}_scaffolds_final.fa" ]; then
        log_error "YaHS未能生成scaffolds文件，请检查日志: logs/yahs.log"
        exit 1
    fi
    
    if [ ! -f "${OUT_PREFIX}_scaffolds_final.agp" ]; then
        log_error "YaHS未能生成AGP文件"
        exit 1
    fi
    
    log_success "YaHS挂载完成"
fi

# =============================================================================
# --- 步骤 4: 生成 Hi-C 可视化文件 ---
# =============================================================================
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
log_info "步骤 4: 生成Hi-C热图文件 (.hic)"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

HIC_FILE="results/yahs_out_final.hic"

if [ -f "${HIC_FILE}" ] && [ $(stat -c%s "${HIC_FILE}") -gt 100000 ]; then
    log_info "发现已有.hic文件 ($(du -h ${HIC_FILE} | cut -f1))，跳过生成步骤"
else
    if [ ! -f "${OUT_PREFIX}.bin" ]; then
        log_error "YaHS .bin文件未找到: ${OUT_PREFIX}.bin"
        exit 1
    fi
    
    log_info "转换YaHS输出为Juicer格式..."
    
    # # 生成比对文件
    # juicer pre "${OUT_PREFIX}.bin" \
    #            "${OUT_PREFIX}_scaffolds_final.agp" \
    #            "${REF_FA}.fai" 2> logs/juicer_pre.log | \
    # sort -k2,2d -k6,6d -T tmp_files/ --parallel=${THREADS} -S32G | \
    # awk 'NF' > tmp_files/alignments_sorted.txt

    # 【修改后】：
    "${YAHS_JUICER_TOOL}" pre "${OUT_PREFIX}.bin" \
               "${OUT_PREFIX}_scaffolds_final.agp" \
               "${REF_FA}.fai" 2> logs/juicer_pre.log | \
    sort -k2,2d -k6,6d -T tmp_files/ --parallel=${THREADS} -S32G | \
    awk 'NF' > tmp_files/alignments_sorted.txt
    
    # 验证中间文件
    if [ ! -s tmp_files/alignments_sorted.txt ]; then
        log_error "alignments_sorted.txt为空，请检查: logs/juicer_pre.log"
        cat logs/juicer_pre.log
        exit 1
    fi
    
    ALIGN_LINES=$(wc -l < tmp_files/alignments_sorted.txt)
    log_info "生成了 ${ALIGN_LINES} 行比对记录"
    
    # 生成染色体大小文件 (基于最终scaffold)
    log_info "生成染色体大小文件..."
    samtools faidx "${OUT_PREFIX}_scaffolds_final.fa"
    cut -f1,2 "${OUT_PREFIX}_scaffolds_final.fa.fai" > results/chrom.sizes.final
    
    # 验证染色体文件
    CHR_COUNT=$(wc -l < results/chrom.sizes.final)
    log_info "检测到 ${CHR_COUNT} 条染色体/scaffold"
    
    if [ ${CHR_COUNT} -eq 0 ]; then
        log_error "染色体大小文件为空"
        exit 1
    fi
    
    # 生成.hic文件
    log_info "生成.hic文件 (可能需要较长时间)..."
    java -Xmx120G -Xms8G -jar "${JUICER_JAR}" pre \
        tmp_files/alignments_sorted.txt \
        "${HIC_FILE}" \
        results/chrom.sizes.final 2>&1 | tee logs/juicer_tools.log
    
    # 验证.hic文件
    if [ -f "${HIC_FILE}" ] && [ $(stat -c%s "${HIC_FILE}") -gt 100000 ]; then
        HIC_SIZE=$(du -h "${HIC_FILE}" | cut -f1)
        log_success ".hic文件生成成功 (大小: ${HIC_SIZE})"
        
        # 清理大型中间文件
        rm -f tmp_files/alignments_sorted.txt
    else
        log_error ".hic文件过小或生成失败 (检查logs/juicer_tools.log)"
        log_warning "保留中间文件用于排查: tmp_files/alignments_sorted.txt"
        exit 1
    fi
fi

# =============================================================================
# --- 步骤 5: 生成 JBAT 文件 (手动纠错用) ---
# =============================================================================
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
log_info "步骤 5: 生成JBAT文件 (用于Juicebox Assembly Tools手动纠错)"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

JBAT_HIC="results/out_JBAT.hic"
JBAT_ASSEMBLY="results/out_JBAT.assembly"

if [ -f "${JBAT_HIC}" ] && [ -s "${JBAT_HIC}" ]; then
    log_info "发现已有JBAT文件，跳过生成步骤"
else
    log_info "生成JBAT格式文件..."
    
    # juicer pre -a -o results/out_JBAT \
    #            "${OUT_PREFIX}.bin" \
    #            "${OUT_PREFIX}_scaffolds_final.agp" \
    #            "${REF_FA}.fai" > logs/out_JBAT.log 2>&1

    # 【修改后】：
    "${YAHS_JUICER_TOOL}" pre -a -o results/out_JBAT \
               "${OUT_PREFIX}.bin" \
               "${OUT_PREFIX}_scaffolds_final.agp" \
               "${REF_FA}.fai" > logs/out_JBAT.log 2>&1
    
    if [ -f "results/out_JBAT.txt" ] && [ -s "results/out_JBAT.txt" ]; then
        # 提取assembly大小信息
        if grep -q "PRE_C_SIZE" logs/out_JBAT.log; then
            grep "PRE_C_SIZE" logs/out_JBAT.log | \
                awk '{print $2" "$3}' > results/jbat_chrom_sizes.txt
        else
            # 备用方案：计算总长度
            TOTAL_BP=$(grep -v '>' "${REF_FA}" | tr -d '\n' | wc -c)
            echo "assembly ${TOTAL_BP}" > results/jbat_chrom_sizes.txt
        fi

        # 先对大文件进行排序 (这是解决42G文件OOM的唯一方法)
        log_info "正在对JBAT文件进行排序 (解决内存溢出关键步骤)..."
        # 使用多线程和60G内存缓存进行排序
        sort --parallel=${THREADS} -S 60G -T tmp_files/ -k2,2d -k6,6d results/out_JBAT.txt > results/out_JBAT_sorted.txt
        
        # 生成JBAT .hic文件
        log_info "生成JBAT .hic文件..."
        java -Xmx120G -Xms8G -jar "${JUICER_JAR}" pre \
            results/out_JBAT_sorted.txt \
            results/out_JBAT.hic.part \
            results/jbat_chrom_sizes.txt 2>&1 | tee logs/juicer_jbat.log
        
        if [ -s "results/out_JBAT.hic.part" ]; then
            mv results/out_JBAT.hic.part "${JBAT_HIC}"
            log_success "JBAT文件生成成功"
            log_info "可使用Juicebox打开以下文件进行手动纠错:"
            log_info "  - ${JBAT_HIC}"
            log_info "  - ${JBAT_ASSEMBLY}"
        else
            log_warning "JBAT .hic文件生成失败，请检查: logs/juicer_jbat.log"
        fi
    else
        log_warning "JBAT文本文件生成失败"
    fi
fi

# =============================================================================
# --- 步骤 6: 组装质量评估 ---
# =============================================================================
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
log_info "步骤 6: 组装质量统计与评估"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

SCAFFOLD_FA="${OUT_PREFIX}_scaffolds_final.fa"

if [ ! -f "${SCAFFOLD_FA}" ]; then
    log_error "Scaffold文件不存在: ${SCAFFOLD_FA}"
    exit 1
fi

# 计算基础统计
log_info "计算基础统计指标..."

# Scaffold数量和总长度
N_SCAFFOLDS=$(grep -c "^>" "${SCAFFOLD_FA}")
TOTAL_LENGTH=$(awk '/^>/ {next} {sum += length($0)} END {print sum}' "${SCAFFOLD_FA}")

# 提取所有scaffold长度并排序
awk '/^>/ {if (seq) print length(seq); seq=""; next} {seq=seq$0} END {if (seq) print length(seq)}' \
    "${SCAFFOLD_FA}" | sort -rn > tmp_files/scaffold_lengths.txt

# 计算N50, N90, L50, L90
awk -v total=${TOTAL_LENGTH} '
BEGIN {
    n50=0; n90=0; l50=0; l90=0
    cum=0; count=0
}
{
    len[NR]=$1
    cum += $1
    count++
    
    if (cum >= total*0.5 && n50==0) {
        n50 = $1
        l50 = count
    }
    if (cum >= total*0.9 && n90==0) {
        n90 = $1
        l90 = count
        exit
    }
}
END {
    print "N50\t" n50
    print "N90\t" n90
    print "L50\t" l50
    print "L90\t" l90
}
' tmp_files/scaffold_lengths.txt > results/assembly_metrics.txt

# 最长scaffold
MAX_SCAFFOLD=$(head -1 tmp_files/scaffold_lengths.txt)

# GC含量
GC_CONTENT=$(awk '/^>/ {next} {
    seq = seq $0
} END {
    gsub(/[^GCgc]/, "", seq)
    gc_count = length(seq)
    gsub(/[^ATGCatgc]/, "", seq)
    total = length(seq)
    if (total > 0) printf "%.2f", (gc_count/total)*100
    else print 0
}' "${SCAFFOLD_FA}")

# 输出统计报告
echo ""
echo "╔════════════════════════════════════════════════════════════╗"
echo "║           📊 组装质量统计报告                              ║"
echo "╠════════════════════════════════════════════════════════════╣"
printf "║ %-30s %28s ║\n" "Scaffold总数:" "${N_SCAFFOLDS}"
printf "║ %-30s %28s ║\n" "总长度:" "$(numfmt --grouping ${TOTAL_LENGTH} 2>/dev/null || echo ${TOTAL_LENGTH}) bp"
printf "║ %-30s %28s ║\n" "最长Scaffold:" "$(numfmt --grouping ${MAX_SCAFFOLD} 2>/dev/null || echo ${MAX_SCAFFOLD}) bp"
printf "║ %-30s %28s ║\n" "GC含量:" "${GC_CONTENT}%"

# 读取并显示N50等指标
while IFS=$'\t' read -r metric value; do
    if [[ $metric == N* ]]; then
        printf "║ %-30s %28s ║\n" "${metric}:" "$(numfmt --grouping ${value} 2>/dev/null || echo ${value}) bp"
    else
        printf "║ %-30s %28s ║\n" "${metric}:" "${value}"
    fi
done < results/assembly_metrics.txt

echo "╚════════════════════════════════════════════════════════════╝"

# 生成长度分布统计
log_info "生成scaffold长度分布..."
awk '
BEGIN {
    bins[0]="<1Kb"; bins[1]="1-10Kb"; bins[2]="10-100Kb"; 
    bins[3]="100Kb-1Mb"; bins[4]="1-10Mb"; bins[5]=">10Mb"
}
{
    len=$1
    if (len<1000) count[0]++
    else if (len<10000) count[1]++
    else if (len<100000) count[2]++
    else if (len<1000000) count[3]++
    else if (len<10000000) count[4]++
    else count[5]++
}
END {
    print "\nScaffold长度分布:"
    print "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    for (i=0; i<6; i++) {
        printf "  %-12s : %6d\n", bins[i], count[i]+0
    }
}
' tmp_files/scaffold_lengths.txt | tee -a results/assembly_metrics.txt

# =============================================================================
# --- 步骤 7: 生成最终报告 ---
# =============================================================================
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
log_info "生成最终分析报告..."
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
