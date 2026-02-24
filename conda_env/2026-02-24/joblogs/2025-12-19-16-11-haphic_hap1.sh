#!/bin/bash

################################################################################
# HapHiC染色体挂载流程脚本
# 用途: 使用Hi-C数据将基因组组装挂载至染色体水平
# 作者: Auto-generated
# 日期: 2025-12-19
################################################################################

set -e  # 遇到错误立即退出
set -u  # 使用未定义变量时报错
set -o pipefail  # 管道命令中任何一个失败都会导致整个管道失败

# ============================================================================
# 参数设置
# ============================================================================

# 显示使用说明
usage() {
    cat << EOF
用法: $0 -f <assembly.fa> -1 <hic_R1.fq.gz> -2 <hic_R2.fq.gz> -n <nchrs> [选项]

必需参数:
    -f    组装的FASTA文件 (可以是contig或scaffold)
    -1    Hi-C测序数据R1 (fastq.gz格式)
    -2    Hi-C测序数据R2 (fastq.gz格式)
    -n    染色体数量

可选参数:
    -o    输出目录 [默认: haphic_output]
    -t    线程数 [默认: 28]
    -p    进程数(用于排序和定向) [默认: 8]
    -r    限制性酶切位点 [默认: GATC]
    -c    组装纠错轮次 [默认: 0, 不纠错]
    -a    移除等位contig间的Hi-C信号(倍性) [默认: 0, 不移除]
    -g    hifiasm的GFA文件(多个文件用逗号分隔)
    -q    使用快速查看模式 [默认: 否]
    -h    显示此帮助信息

示例:
    # 基础用法
    $0 -f asm.fa -1 hic_R1.fq.gz -2 hic_R2.fq.gz -n 24
    
    # 使用组装纠错和多线程
    $0 -f asm.fa -1 hic_R1.fq.gz -2 hic_R2.fq.gz -n 24 -c 2 -t 32 -p 16
    
    # 四倍体单倍型分型组装
    $0 -f asm.fa -1 hic_R1.fq.gz -2 hic_R2.fq.gz -n 24 -a 4
    
    # 使用Arima双酶方案
    $0 -f asm.fa -1 hic_R1.fq.gz -2 hic_R2.fq.gz -n 24 -r "GATC,GANTC"

EOF
    exit 1
}

# ============================================================================
# 解析命令行参数
# ============================================================================

# 默认参数
ASSEMBLY="OV53_1.hap1.primary.fa"
HIC_R1="OV53_1-hic_R1.fastq.gz"
HIC_R2="OV53_1-hic_R1.fastq.gz"
NCHRS="12"
OUTDIR="haphic_output"
THREADS=88
PROCESSES=88
RE_SITE="GATC"
CORRECT_ROUNDS=2
PLOIDY=2
GFA_FILES=""
QUICK_VIEW=false

# 解析参数
while getopts "f:1:2:n:o:t:p:r:c:a:g:qh" opt; do
    case $opt in
        f) ASSEMBLY="$OPTARG" ;;
        1) HIC_R1="$OPTARG" ;;
        2) HIC_R2="$OPTARG" ;;
        n) NCHRS="$OPTARG" ;;
        o) OUTDIR="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        p) PROCESSES="$OPTARG" ;;
        r) RE_SITE="$OPTARG" ;;
        c) CORRECT_ROUNDS="$OPTARG" ;;
        a) PLOIDY="$OPTARG" ;;
        g) GFA_FILES="$OPTARG" ;;
        q) QUICK_VIEW=true ;;
        h) usage ;;
        *) usage ;;
    esac
done

# 检查必需参数
if [[ -z "$ASSEMBLY" ]] || [[ -z "$HIC_R1" ]] || [[ -z "$HIC_R2" ]] || [[ -z "$NCHRS" ]]; then
    echo "错误: 缺少必需参数!"
    usage
fi

# 检查输入文件是否存在
if [[ ! -f "$ASSEMBLY" ]]; then
    echo "错误: 文件不存在: $ASSEMBLY"
    exit 1
fi

if [[ ! -f "$HIC_R1" ]]; then
    echo "错误: 文件不存在: $HIC_R1"
    exit 1
fi

if [[ ! -f "$HIC_R2" ]]; then
    echo "错误: 文件不存在: $HIC_R2"
    exit 1
fi

# ============================================================================
# 软件路径设置
# ============================================================================

HAPHIC="/share/org/YZWL/yzwl_lixg/miniforge3/envs/haphic/bin/haphic"
BWA="bwa"
SAMTOOLS="samtools"
SAMBLASTER="samblaster"
FILTER_BAM="/share/org/YZWL/yzwl_lixg/miniforge3/envs/haphic/bin/../utils/filter_bam"

# 检查软件是否可用
for cmd in "$BWA" "$SAMTOOLS" "$SAMBLASTER"; do
    if ! command -v $cmd &> /dev/null; then
        echo "错误: 未找到命令 $cmd, 请确保已安装并添加到PATH"
        exit 1
    fi
done

if [[ ! -f "$HAPHIC" ]]; then
    echo "错误: 未找到HapHiC程序: $HAPHIC"
    exit 1
fi

# ============================================================================
# 创建输出目录
# ============================================================================

# 获取输入文件的绝对路径(在切换目录之前)
ASSEMBLY_ABS=$(readlink -f "$ASSEMBLY")
HIC_R1_ABS=$(readlink -f "$HIC_R1")
HIC_R2_ABS=$(readlink -f "$HIC_R2")

mkdir -p "$OUTDIR"
cd "$OUTDIR"

# 记录日志
LOG_FILE="haphic_pipeline.log"
exec > >(tee -a "$LOG_FILE")
exec 2>&1

echo "============================================================================"
echo "HapHiC染色体挂载流程开始"
echo "开始时间: $(date)"
echo "============================================================================"
echo "输入参数:"
echo "  组装文件: $ASSEMBLY"
echo "  Hi-C R1: $HIC_R1"
echo "  Hi-C R2: $HIC_R2"
echo "  染色体数: $NCHRS"
echo "  输出目录: $OUTDIR"
echo "  线程数: $THREADS"
echo "  进程数: $PROCESSES"
echo "  限制性酶切位点: $RE_SITE"
echo "  纠错轮次: $CORRECT_ROUNDS"
echo "  倍性(移除等位信号): $PLOIDY"
echo "  GFA文件: ${GFA_FILES:-未指定}"
echo "  快速查看模式: $QUICK_VIEW"
echo "============================================================================"
echo ""

# ============================================================================
# 步骤1: 创建BWA索引
# ============================================================================

echo "[$(date)] 步骤1: 为组装文件创建BWA索引..."

# 复制组装文件到工作目录
cp "$ASSEMBLY_ABS" asm.fa

if [[ ! -f "asm.fa.bwt" ]]; then
    $BWA index asm.fa
    echo "[$(date)] BWA索引创建完成"
else
    echo "[$(date)] BWA索引已存在,跳过"
fi

echo ""

# ============================================================================
# 步骤2: Hi-C数据比对
# ============================================================================

echo "[$(date)] 步骤2: 将Hi-C数据比对到组装序列..."

if [[ ! -f "HiC.bam" ]]; then
    $BWA mem -5SP -t $THREADS asm.fa "$HIC_R1_ABS" "$HIC_R2_ABS" | \
        $SAMBLASTER | \
        $SAMTOOLS view - -@ $((THREADS/2)) -S -h -b -F 3340 -o HiC.bam
    echo "[$(date)] Hi-C数据比对完成"
else
    echo "[$(date)] HiC.bam已存在,跳过比对"
fi

echo ""

# ============================================================================
# 步骤3: 过滤BAM文件
# ============================================================================

echo "[$(date)] 步骤3: 过滤BAM文件 (MAPQ>=1, NM<3)..."

if [[ ! -f "HiC.filtered.bam" ]]; then
    $FILTER_BAM HiC.bam 1 --nm 3 --threads $((THREADS/2)) | \
        $SAMTOOLS view - -b -@ $((THREADS/2)) -o HiC.filtered.bam
    echo "[$(date)] BAM文件过滤完成"
else
    echo "[$(date)] HiC.filtered.bam已存在,跳过过滤"
fi

echo ""

# ============================================================================
# 步骤4: 运行HapHiC pipeline
# ============================================================================

echo "[$(date)] 步骤4: 运行HapHiC挂载流程..."

# 构建HapHiC命令
HAPHIC_CMD="$HAPHIC pipeline asm.fa HiC.filtered.bam $NCHRS"
HAPHIC_CMD="$HAPHIC_CMD --threads $THREADS"
HAPHIC_CMD="$HAPHIC_CMD --processes $PROCESSES"
HAPHIC_CMD="$HAPHIC_CMD --RE \"$RE_SITE\""

# 添加可选参数
if [[ $CORRECT_ROUNDS -gt 0 ]]; then
    HAPHIC_CMD="$HAPHIC_CMD --correct_nrounds $CORRECT_ROUNDS"
    echo "  启用组装纠错: $CORRECT_ROUNDS 轮"
fi

if [[ $PLOIDY -gt 0 ]]; then
    HAPHIC_CMD="$HAPHIC_CMD --remove_allelic_links $PLOIDY"
    echo "  移除等位contig间的Hi-C信号 (倍性: $PLOIDY)"
fi

if [[ -n "$GFA_FILES" ]]; then
    HAPHIC_CMD="$HAPHIC_CMD --gfa \"$GFA_FILES\""
    echo "  使用GFA文件: $GFA_FILES"
fi

if [[ "$QUICK_VIEW" = true ]]; then
    HAPHIC_CMD="$HAPHIC_CMD --quick_view"
    echo "  使用快速查看模式"
fi

echo "执行命令: $HAPHIC_CMD"
echo ""

eval $HAPHIC_CMD

echo "[$(date)] HapHiC挂载流程完成"
echo ""

# ============================================================================
# 步骤5: 整理输出结果
# ============================================================================

echo "[$(date)] 步骤5: 整理输出结果..."

# 创建结果目录
mkdir -p final_results

# 复制主要输出文件
if [[ -f "04.build/scaffolds.fa" ]]; then
    cp 04.build/scaffolds.fa final_results/
    cp 04.build/scaffolds.agp final_results/
    cp 04.build/scaffolds.raw.agp final_results/
    echo "  已复制 scaffolds.fa, scaffolds.agp, scaffolds.raw.agp"
fi

# 如果存在纠错后的contig
if [[ -f "01.cluster/corrected_asm.fa" ]]; then
    cp 01.cluster/corrected_asm.fa final_results/
    echo "  已复制 corrected_asm.fa"
fi

# 复制Juicebox脚本
if [[ -f "04.build/juicebox.sh" ]]; then
    cp 04.build/juicebox.sh final_results/
    echo "  已复制 juicebox.sh"
fi

echo ""

# ============================================================================
# 统计结果
# ============================================================================

echo "[$(date)] 生成结果统计..."

if [[ -f "final_results/scaffolds.fa" ]]; then
    echo ""
    echo "============================================================================"
    echo "最终Scaffold统计信息:"
    echo "============================================================================"
    
    # 统计scaffold数量和总长度
    awk '/^>/ {if (seqlen) print seqlen; seq=""; seqlen=0; next} 
         {seqlen += length($0)} 
         END {if (seqlen) print seqlen}' final_results/scaffolds.fa | \
        awk '{sum+=$1; count++; lens[count]=$1} 
             END {
                 asort(lens);
                 print "Scaffold数量:", count;
                 print "总长度:", sum, "bp";
                 print "最长Scaffold:", lens[count], "bp";
                 print "最短Scaffold:", lens[1], "bp";
                 
                 # 计算N50
                 half=sum/2;
                 cumsum=0;
                 for(i=count; i>=1; i--) {
                     cumsum+=lens[i];
                     if(cumsum>=half) {
                         print "Scaffold N50:", lens[i], "bp";
                         break;
                     }
                 }
             }'
    echo "============================================================================"
fi

echo ""

# ============================================================================
# 完成
# ============================================================================

echo "============================================================================"
echo "HapHiC染色体挂载流程全部完成!"
echo "结束时间: $(date)"
echo "============================================================================"
echo ""
echo "主要输出文件位于: $OUTDIR/final_results/"
echo "  - scaffolds.fa: 最终scaffold序列"
echo "  - scaffolds.agp: SALSA风格的AGP文件"
echo "  - scaffolds.raw.agp: YaHS风格的AGP文件"
if [[ $CORRECT_ROUNDS -gt 0 ]]; then
    echo "  - corrected_asm.fa: 纠错后的contig序列"
fi
echo "  - juicebox.sh: Juicebox可视化脚本"
echo ""
echo "后续步骤:"
echo "1. 使用Juicebox进行可视化和手动调整:"
echo "   cd final_results && bash juicebox.sh"
echo ""
echo "2. 使用HapHiC生成Hi-C互作热图:"
echo "   $HAPHIC plot final_results/scaffolds.raw.agp HiC.filtered.bam"
echo ""
echo "完整日志保存在: $OUTDIR/$LOG_FILE"
echo "============================================================================"