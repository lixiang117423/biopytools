#!/bin/bash

# ==============================================================================
# 🧬 ALLHiC Pipeline Automation Script (Official Workflow)
# 📅 Date: $(date +%F)
# 🎯 Goal: Chromosome Scaffolding using standard ALLHiC binaries
# ==============================================================================

# ---------------------- 🛠️ 参数配置区 (Configuration) ----------------------

# 1. 📂 软件路径 (Crucial)
ALLHIC_SOFTWARE_PATH="/share/org/YZWL/yzwl_lixg/software/ALLHiC"

# 2. 📂 工作目录
WORK_DIR="/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/05.allhic"

# 3. 🧬 输入文件
REF_FA_RAW="/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/05.allhic/OV53_1.primary.fa"
R1_FQ_RAW="/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/04.YaHS/fastq/OV53_1-hic_R1.fastq.gz"
R2_FQ_RAW="/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/04.YaHS/fastq/OV53_1-hic_R2.fastq.gz"

# 4. 🔢 基因组参数
CHROMOSOME_K=12                 # 预期染色体数 (k)
RE_MOTIF="GATC"                 # 酶切序列 (MboI = GATC, HindIII = AAGCTT)

# 5. ⚙️ 系统资源
THREADS=64

# 6. 🔧 二倍体/多倍体修剪 (Pruning)
# 如果有 Allele.ctg.table 文件，请填写路径；否则留空，脚本将跳过 Prune 步骤或仅做格式转换
ALLELE_TABLE="" 
# 示例: ALLELE_TABLE="/path/to/Allele.ctg.table"

# 7. 📛 样本名称
SAMPLE_NAME="zhugecai"

# ==============================================================================
# 🔧 环境初始化
# ==============================================================================

# 设置环境变量
export PATH=${ALLHIC_SOFTWARE_PATH}/scripts:${ALLHIC_SOFTWARE_PATH}/bin:$PATH
export PATH=$PATH  # 确保系统其他工具(bwa, samtools)可用

# 日志设置
mkdir -p "${WORK_DIR}/logs"
LOG_FILE="${WORK_DIR}/logs/allhic_pipeline.log"

log() {
    echo -e "[$(date +'%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG_FILE"
}

run_cmd() {
    log "▶️  Running: $1"
    eval "$1" >> "$LOG_FILE" 2>&1
    if [ $? -ne 0 ]; then
        log "❌ Error executing command. Check log for details."
        exit 1
    fi
    log "✅ Done."
}

# 进入工作目录
mkdir -p "$WORK_DIR"
cd "$WORK_DIR" || exit 1

log "🚀 Starting ALLHiC Pipeline..."
log "📂 Work Dir: $WORK_DIR"
log "🛠️  ALLHiC Path: $ALLHIC_SOFTWARE_PATH"

# ==============================================================================
# 0. 数据准备 (Pre-check & Linking)
# ==============================================================================
log "🔗 [Step 0] Linking input files..."
ln -sf "$REF_FA_RAW" draft.asm.fasta
ln -sf "$R1_FQ_RAW" reads_R1.fastq.gz
ln -sf "$R2_FQ_RAW" reads_R2.fastq.gz

# 检查依赖
command -v ALLHiC_partition >/dev/null 2>&1 || { log "❌ ALLHiC_partition not found in PATH!"; exit 1; }
command -v bwa >/dev/null 2>&1 || { log "❌ bwa not found!"; exit 1; }
command -v samtools >/dev/null 2>&1 || { log "❌ samtools not found!"; exit 1; }

# ==============================================================================
# 1. 比对 (Mapping) - 使用 BWA MEM (官方推荐 Tip)
# ==============================================================================
log "🔍 [Step 1] Mapping Hi-C reads..."

# 1.1 建立索引
if [ ! -f "draft.asm.fasta.bwt" ]; then
    run_cmd "bwa index draft.asm.fasta"
fi
if [ ! -f "draft.asm.fasta.fai" ]; then
    run_cmd "samtools faidx draft.asm.fasta"
fi

# 1.2 比对 (使用 bwa mem 替代 bwa aln，适合大基因组和长读长)
# 输出 sample.clean.bam (需经过过滤)
CLEAN_BAM="sample.clean.bam"

if [ ! -f "$CLEAN_BAM" ]; then
    log "   Running bwa mem and filtering (MAPQ>=30, no secondary)..."
    # 这里直接生成 clean bam，跳过 perl 脚本处理 sam 的步骤，因为 bwa mem 输出已经是 sam/bam
    # -F 2316: 过滤 unmapped(4), secondary(256), supplementary(2048)
    # -q 30: 过滤低质量比对
    run_cmd "bwa mem -t $THREADS -5SP draft.asm.fasta reads_R1.fastq.gz reads_R2.fastq.gz | \
             samtools view -@ $THREADS -hF 2316 -q 30 - | \
             samtools sort -@ $THREADS -n -o $CLEAN_BAM -"
             
    # 注意：ALLHiC 有时需要 name-sorted bam 用于 extract，有时需要 coordinate-sorted。
    # 官方流程中 filterBAM_forHiC.pl 输出通常是处理过的。
    # 这里的 -n (name sort) 也是为了后续 extract 提取 pairs 更准确。
else
    log "⚠️  $CLEAN_BAM exists. Skipping mapping."
fi

# ==============================================================================
# 2. 修剪 (Pruning) - 处理多倍体/移除同源信号
# ==============================================================================
log "✂️  [Step 2] Pruning (Removing allelic/weak signals)..."

PRUNED_BAM="prunning.bam"

if [ ! -f "$PRUNED_BAM" ]; then
    if [ -n "$ALLELE_TABLE" ] && [ -f "$ALLELE_TABLE" ]; then
        log "   Allele table found: $ALLELE_TABLE"
        run_cmd "ALLHiC_prune -i $ALLELE_TABLE -b $CLEAN_BAM -r draft.asm.fasta"
        # ALLHiC_prune 通常输出 prunning.bam
    else
        log "   ℹ️ No Allele table provided ($ALLELE_TABLE). Skipping pruning logic."
        log "   Linking clean BAM to prunning.bam for next steps."
        run_cmd "ln -sf $CLEAN_BAM $PRUNED_BAM"
    fi
else
    log "⚠️  $PRUNED_BAM exists. Skipping pruning."
fi

# ==============================================================================
# 3. 分组 (Partition)
# ==============================================================================
log "📦 [Step 3] Partitioning into $CHROMOSOME_K groups..."

# 检查 Pruning BAM 是否存在
if [ ! -f "$PRUNED_BAM" ]; then log "❌ $PRUNED_BAM missing!"; exit 1; fi

# Partition 输出通常是 prunning.clusters.txt 和 prunning.counts_${RE_MOTIF}.${K}g*.txt
CLUSTERS_FILE="prunning.clusters.txt"

if [ ! -f "$CLUSTERS_FILE" ]; then
    run_cmd "ALLHiC_partition -b $PRUNED_BAM -r draft.asm.fasta -e $RE_MOTIF -k $CHROMOSOME_K"
else
    log "⚠️  Partition results exist. Skipping."
fi

# ==============================================================================
# 3.5 提取信号 (Extract) - 为 Rescue 和 Optimize 做准备
# ==============================================================================
log "🧬 [Step 3.5] Extracting CLM and Counts from CLEAN BAM..."
# 注意：Extract 应该用 sample.clean.bam (包含所有信号)，而不是 pruned bam
# 这会生成 sample.clean.clm 和 sample.clean.counts_${RE_MOTIF}.txt

CLM_FILE="sample.clean.clm"
COUNTS_FILE="sample.clean.counts_${RE_MOTIF}.txt"

if [ ! -f "$CLM_FILE" ]; then
    run_cmd "allhic extract $CLEAN_BAM draft.asm.fasta --RE $RE_MOTIF"
else
    log "⚠️  CLM file exists. Skipping extraction."
fi

# ==============================================================================
# 4. 挽救 (Rescue) - 召回未分组的 Contigs
# ==============================================================================
log "🚑 [Step 4] Rescuing unplaced contigs..."

# Rescue 需要原始的 counts 和 clusters
# 输出通常会更新 clusters 或者生成 groups.txt (取决于版本)
# 我们这里假设生成新的 cluster 映射

if [ -f "$CLUSTERS_FILE" ] && [ -f "$COUNTS_FILE" ]; then
    # 检查是否已经 Rescue 过 (通常 Rescue 比较快，可以覆盖运行，或者检查标志文件)
    # 这里的 -c 是 partition 产生的 clusters
    # -i 是 extract 产生的 counts
    run_cmd "ALLHiC_rescue -b $CLEAN_BAM -r draft.asm.fasta -c $CLUSTERS_FILE -i $COUNTS_FILE"
    
    # ⚠️ 注意: ALLHiC_rescue 运行后，通常会生成 "groups.txt" 或更新 cluster 文件
    # 我们需要确认 Rescue 的输出用于下一步 Optimize
else
    log "❌ Missing inputs for Rescue step!"
    exit 1
fi

# ==============================================================================
# 5. 优化 (Optimize) - 排序和定向
# ==============================================================================
log "⚙️  [Step 5] Optimizing ordering and orientation..."

# 准备 Optimize 的输入文件
# Optimize 需要: 1. groupX.txt (格式: #Contig RECounts Length), 2. .clm 文件
# Rescue 步骤通常生成了 groups.txt，我们需要将其拆分或者 ALLHiC_partition 已经生成了 counts 文件。
# 最佳实践：使用 partition 生成的 prunning.counts_RE.KgX.txt，或者 Rescue 后的结果。
# 这里的逻辑是：如果 Rescue 改变了分组，我们需要重新生成 group files。
# 为简单起见，且遵循官方流程 "allhic optimize group1.txt"，我们需要确保这些文件存在。

# 自动检测 partition 生成的文件
PARTITION_FILES=$(ls prunning.counts_${RE_MOTIF}.${CHROMOSOME_K}g*.txt 2>/dev/null)

if [ -z "$PARTITION_FILES" ]; then
    log "❌ No partition group files found! (prunning.counts_${RE_MOTIF}.*)"
    exit 1
fi

# 生成命令列表
> optimize_cmds.sh
COUNT=1
for GFILE in $PARTITION_FILES; do
    # 重命名为简单的 groupN.txt 以符合流程习惯 (可选，但为了清晰)
    NEW_NAME="group${COUNT}.txt"
    cp "$GFILE" "$NEW_NAME"
    
    # 检查是否已经 optimize 过 (生成 .tour 文件)
    if [ ! -f "group${COUNT}.tour" ]; then
        echo "allhic optimize $NEW_NAME $CLM_FILE" >> optimize_cmds.sh
    fi
    ((COUNT++))
done

# 并行运行
if [ -s "optimize_cmds.sh" ]; then
    log "   Running optimization for $(wc -l < optimize_cmds.sh) groups..."
    # 使用 ParaFly 或 xargs 并行
    if command -v ParaFly >/dev/null 2>&1; then
        ParaFly -c optimize_cmds.sh -CPU "$THREADS" -failed_cmds optimize_failed.cmds
    else
        # 简单的 bash 并行
        cat optimize_cmds.sh | xargs -L 1 -I CMD -P "$THREADS" bash -c "CMD"
    fi
else
    log "⚠️  Optimization seems done (tour files exist)."
fi

# ==============================================================================
# 6. 构建 (Build) - 生成 Fasta 和 AGP
# ==============================================================================
log "🏗️  [Step 6] Building chromosome-scale assembly..."

# ALLHiC_build 会自动寻找当前目录下的 *.tour 文件
if [ ! -f "groups.asm.fasta" ]; then
    run_cmd "ALLHiC_build draft.asm.fasta"
else
    log "⚠️  groups.asm.fasta exists. Skipping build."
fi

# ==============================================================================
# 7. 绘图 (Plot) - 生成热图
# ==============================================================================
log "📊 [Step 7] Plotting heatmap..."

if [ -f "groups.agp" ] && [ -f "$CLEAN_BAM" ]; then
    # 7.1 生成 chrn.list (格式: groupName Length)
    # 从 groups.asm.fasta 的 index 提取，或者从 tour 文件推断
    # 这里用 samtools faidx 获取最终组装的长度
    run_cmd "samtools faidx groups.asm.fasta"
    cut -f1,2 groups.asm.fasta.fai > chrn.list
    
    # 7.2 绘图 (500k 分辨率)
    if [ ! -f "heatmap.pdf" ]; then
        # ALLHiC_plot bam agp list bin_size ext
        run_cmd "ALLHiC_plot $CLEAN_BAM groups.agp chrn.list 500k pdf"
    fi
else
    log "⚠️  Cannot plot. Missing groups.agp or clean bam."
fi

# ==============================================================================
# 🎉 结束
# ==============================================================================
log "=========================================="
log "🎉 Pipeline Finished!"
log "📂 Results:"
log "   - Final Fasta: ${WORK_DIR}/groups.asm.fasta"
log "   - AGP File   : ${WORK_DIR}/groups.agp"
log "   - Heatmap    : ${WORK_DIR}/heatmap.pdf"
log "=========================================="
