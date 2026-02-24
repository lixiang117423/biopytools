#!/bin/bash

# ==============================================================================
# 🧬 ALLHiC Pipeline Automation Script (Bio-protocol version)
# 📅 Date: $(date +%F)
# 🎯 Goal: Chromosome Scaffolding for Diploid Genome
# ==============================================================================

# ---------------------- 🛠️ 参数配置区 (Configuration) ----------------------

# 1. 📂 工作目录
WORK_DIR="/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/05.allhic"

# 2. 🧬 输入文件 (原始路径)
REF_FA_RAW="/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/05.allhic/OV53_1.primary.fa"
R1_FQ_RAW="/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/04.YaHS/fastq/OV53_1-hic_R1.fastq.gz"
R2_FQ_RAW="/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/04.YaHS/fastq/OV53_1-hic_R2.fastq.gz"

# 3. 🔢 基因组参数
CHROMOSOME_K=12                 # 预期染色体数 (k)
RE_NAME="MBOI"                  # 酶名称 (用于 PreprocessSAMs.pl)
RE_MOTIF="GATC"                 # 酶切序列 (MboI = GATC, 用于 Partition)

# 4. ⚙️ 系统资源
THREADS=64                      # CPU 线程数 (根据服务器情况调整)

# ==============================================================================

# 🛑 错误处理函数
set -e
handle_error() {
    echo -e "\n❌ [ERROR] 脚本在第 $1 行发生错误! 请检查日志。"
    exit 1
}
trap 'handle_error $LINENO' ERR

# 🚀 开始流程
echo -e "\n🚀 ================= STARTING ALLHiC PIPELINE ================="
echo "📂 Work Dir : $WORK_DIR"
echo "🧬 Ref Fasta: $REF_FA_RAW"
echo "✂️  Enzyme   : $RE_NAME ($RE_MOTIF)"
echo "🔢 Target K : $CHROMOSOME_K"
echo -e "=============================================================\n"

# 创建并进入工作目录
mkdir -p "$WORK_DIR"
cd "$WORK_DIR"

# 🔗 Step 0: 建立软链接 (保持目录整洁)
echo -e "🔗 [Step 0] Linking input files..."
ln -sf "$REF_FA_RAW" draft.asm.fasta
ln -sf "$R1_FQ_RAW" Lib_R1.fastq.gz
ln -sf "$R2_FQ_RAW" Lib_R2.fastq.gz

# ------------------------------------------------------------------------------
# 🔧 Step B: Correction of the draft contigs (纠错)
# ------------------------------------------------------------------------------
echo -e "\n🔧 [Step B] Correcting draft contigs (ALLHiC_corrector)..."

if [ ! -f "seq.HiCcorrected.fasta" ]; then
    # B.1.a 建索引
    echo "   Building index for draft assembly..."
    bwa index draft.asm.fasta
    samtools faidx draft.asm.fasta

    # B.1.b 比对 Hi-C reads
    echo "   Mapping reads to draft assembly..."
    bwa mem -5SPM -t "$THREADS" draft.asm.fasta Lib_R1.fastq.gz Lib_R2.fastq.gz \
        | samtools view -hF 256 - \
        | samtools sort -@ "$THREADS" -o sorted.bam -T tmp.ali

    samtools index -@ "$THREADS" sorted.bam

    # B.2 运行纠错
    echo "   Running ALLHiC_corrector..."
    ALLHiC_corrector -m sorted.bam -r draft.asm.fasta -o seq.HiCcorrected.fasta -t "$THREADS"
    echo "✅ Step B Done. Corrected assembly: seq.HiCcorrected.fasta"
else
    echo "⚠️  seq.HiCcorrected.fasta already exists. Skipping Step B."
fi

# ------------------------------------------------------------------------------
# 🔍 Step C: Map Hi-C reads to corrected assembly (重新比对与过滤)
# ------------------------------------------------------------------------------
echo -e "\n🔍 [Step C] Remapping and filtering Hi-C signals..."

CORRECTED_FA="seq.HiCcorrected.fasta"
PRUNED_BAM="sample.unique.REduced.paired_only.bam"

if [ ! -f "$PRUNED_BAM" ]; then
    # C.1 对纠错后的基因组建索引
    echo "   Indexing corrected assembly..."
    bwa index "$CORRECTED_FA"

    # 比对
    echo "   Mapping reads to corrected assembly..."
    bwa mem -5SPM -t "$THREADS" "$CORRECTED_FA" Lib_R1.fastq.gz Lib_R2.fastq.gz \
        | samtools view -hF 256 - \
        | samtools sort -@ "$THREADS" -o sample.bwa_mem.bam -T tmp.ali2

    # C.2 过滤 MAPQ < 30
    echo "   Filtering alignments (MAPQ >= 30)..."
    samtools view -bq 30 sample.bwa_mem.bam > sample.unique.bam

    # C.3 预处理 BAM (Pruning)
    # 注意: 这里使用 MBOI 作为参数
    echo "   Pruning BAM file using PreprocessSAMs.pl (Enzyme: $RE_NAME)..."
    PreprocessSAMs.pl sample.unique.bam "$CORRECTED_FA" "$RE_NAME"
    
    echo "✅ Step C Done. Pruned BAM: $PRUNED_BAM"
else
    echo "⚠️  $PRUNED_BAM already exists. Skipping Step C."
fi

# ------------------------------------------------------------------------------
# 📦 Step D: Partition (聚类)
# ------------------------------------------------------------------------------
echo -e "\n📦 [Step D] Partitioning contigs into $CHROMOSOME_K groups..."

# 注意: 这里 ALLHiC_partition 通常使用 Motif (GATC)
# 文档中示例是 HINDIII (AAGCTT)，此处用 -e GATC 适配 MboI
ALLHiC_partition -r "$CORRECTED_FA" -e "$RE_MOTIF" -k "$CHROMOSOME_K" -b "$PRUNED_BAM"

echo "✅ Step D Done. Check *.counts_*.txt files."

# ------------------------------------------------------------------------------
# ⚙️ Step E: Optimization (排序与定向)
# ------------------------------------------------------------------------------
echo -e "\n⚙️  [Step E] Optimizing ordering and orientation..."

# 生成 optimize 命令列表
# 这一步 ALLHiC 会生成 .clm 文件，我们需要找到它
CLM_FILE=$(ls sample.unique.REduced.paired_only.clm 2>/dev/null || echo "")

if [ -z "$CLM_FILE" ]; then
    echo "❌ Error: .clm file not found! Partition step might have failed."
    exit 1
fi

echo "   Generating command list..."
> cmd.list
for i in $(seq 1 "$CHROMOSOME_K"); do
    # 查找对应的 count 文件
    COUNT_FILE=$(ls sample.unique.REduced.paired_only.counts_${RE_MOTIF}.${CHROMOSOME_K}g${i}.txt 2>/dev/null)
    
    if [ -n "$COUNT_FILE" ]; then
        echo "allhic optimize $COUNT_FILE $CLM_FILE" >> cmd.list
    else
        echo "⚠️  Warning: Count file for group $i not found."
    fi
done

# 使用 ParaFly 并行运行
echo "   Running Optimization with ParaFly (Threads: $THREADS)..."
ParaFly -c cmd.list -CPU "$THREADS"

echo "✅ Step E Done."

# ------------------------------------------------------------------------------
# 🏗️ Step F: Building (构建最终序列)
# ------------------------------------------------------------------------------
echo -e "\n🏗️  [Step F] Building chromosome-scale scaffolds..."

ALLHiC_build "$CORRECTED_FA"

echo -e "\n🎉 ================= ALLHiC PIPELINE FINISHED ================= 🎉"
echo "📊 最终结果文件:"
echo "   1. groups.agp  (Contig 位置信息)"
echo "   2. groups.asm.fasta (最终染色体序列)"
echo -e "==============================================================="
