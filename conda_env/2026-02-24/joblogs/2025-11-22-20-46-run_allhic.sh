# #!/bin/bash

# # ==============================================================================
# # 🧬 ALLHiC Pipeline Automation Script (Updated for New Version)
# # 📅 Date: $(date +%F)
# # 🎯 Goal: Chromosome Scaffolding for Diploid Genome
# # ==============================================================================

# # ---------------------- 🛠️ 参数配置区 (Configuration) ----------------------

# # 1. 📂 工作目录
# WORK_DIR="/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/05.allhic"

# # 2. 🧬 输入文件 (原始路径)
# REF_FA_RAW="/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/05.allhic/OV53_1.primary.fa"
# R1_FQ_RAW="/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/04.YaHS/fastq/OV53_1-hic_R1.fastq.gz"
# R2_FQ_RAW="/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/04.YaHS/fastq/OV53_1-hic_R2.fastq.gz"

# # 3. 🔢 基因组参数
# CHROMOSOME_K=12                 # 预期染色体数 (k)
# RE_MOTIF="GATC"                 # 酶切序列 (MboI = GATC)

# # 4. ⚙️ 系统资源
# THREADS=64                      # CPU 线程数

# # 5. 📛 样本名称
# SAMPLE_NAME="sample"            # 输出文件前缀

# # ==============================================================================

# # 🛑 错误处理函数
# set -e
# handle_error() {
#     echo -e "\n❌ [ERROR] 脚本在第 $1 行发生错误! 请检查日志。"
#     exit 1
# }
# trap 'handle_error $LINENO' ERR

# # 🚀 开始流程
# echo -e "\n🚀 ================= STARTING ALLHiC PIPELINE ================="
# echo "📂 Work Dir : $WORK_DIR"
# echo "🧬 Ref Fasta: $REF_FA_RAW"
# echo "✂️  Enzyme   : MboI ($RE_MOTIF)"
# echo "🔢 Target K : $CHROMOSOME_K"
# echo "🧵 Threads  : $THREADS"
# echo -e "=============================================================\n"

# # 创建并进入工作目录
# mkdir -p "$WORK_DIR"
# cd "$WORK_DIR"

# # 🔗 Step 0: 建立软链接
# echo -e "🔗 [Step 0] Linking input files..."
# ln -sf "$REF_FA_RAW" draft.asm.fasta
# ln -sf "$R1_FQ_RAW" Lib_R1.fastq.gz
# ln -sf "$R2_FQ_RAW" Lib_R2.fastq.gz
# echo "✅ Step 0 Done."

# # ------------------------------------------------------------------------------
# # 🔍 Step 1: Build Index and Map Hi-C Reads
# # ------------------------------------------------------------------------------
# echo -e "\n🔍 [Step 1] Building index and mapping Hi-C reads..."

# if [ ! -f "${SAMPLE_NAME}.bam" ]; then
#     # 1.1 建索引
#     if [ ! -f "draft.asm.fasta.bwt" ]; then
#         echo "   Building BWA index..."
#         bwa index draft.asm.fasta
#     fi
    
#     if [ ! -f "draft.asm.fasta.fai" ]; then
#         echo "   Building samtools index..."
#         samtools faidx draft.asm.fasta
#     fi

#     # 1.2 比对 Hi-C reads
#     echo "   Mapping Hi-C reads (this may take a while)..."
#     bwa mem -5SPM -t "$THREADS" draft.asm.fasta \
#         Lib_R1.fastq.gz Lib_R2.fastq.gz \
#         | samtools view -@ "$THREADS" -bS - \
#         | samtools sort -@ "$THREADS" -o ${SAMPLE_NAME}.bam -

#     # 1.3 索引 BAM
#     echo "   Indexing BAM file..."
#     samtools index -@ "$THREADS" ${SAMPLE_NAME}.bam
    
#     echo "✅ Step 1 Done. BAM file: ${SAMPLE_NAME}.bam"
# else
#     echo "⚠️  ${SAMPLE_NAME}.bam already exists. Skipping Step 1."
# fi

# # ------------------------------------------------------------------------------
# # 🔪 Step 2: Prune BAM (Filter and Extract Hi-C Links)
# # ------------------------------------------------------------------------------
# echo -e "\n🔪 [Step 2] Pruning BAM file (allhic prune)..."

# PRUNED_BAM="${SAMPLE_NAME}.clean.bam"

# if [ ! -f "$PRUNED_BAM" ]; then
#     # 2.1 Filter MAPQ < 1 (保留高质量比对)
#     echo "   Filtering low-quality alignments..."
#     samtools view -@ "$THREADS" -bq 1 ${SAMPLE_NAME}.bam > ${SAMPLE_NAME}.filtered.bam
    
#     # 2.2 运行 allhic prune
#     # 注意: 如果是二倍体且有 allele table, 可以添加 --alleles alleles.table
#     echo "   Running allhic prune..."
#     allhic prune --bam ${SAMPLE_NAME}.filtered.bam \
#                  --fasta draft.asm.fasta \
#                  --output $PRUNED_BAM
    
#     echo "✅ Step 2 Done. Clean BAM: $PRUNED_BAM"
# else
#     echo "⚠️  $PRUNED_BAM already exists. Skipping Step 2."
# fi

# # ------------------------------------------------------------------------------
# # 🧬 Step 3: Extract Hi-C Link Information
# # ------------------------------------------------------------------------------
# echo -e "\n🧬 [Step 3] Extracting Hi-C link information..."

# CLM_FILE="${SAMPLE_NAME}.clm"

# if [ ! -f "$CLM_FILE" ]; then
#     echo "   Running allhic extract..."
#     allhic extract --bam $PRUNED_BAM \
#                    --fasta draft.asm.fasta \
#                    --RE "$RE_MOTIF" \
#                    --output $CLM_FILE
    
#     echo "✅ Step 3 Done. CLM file: $CLM_FILE"
# else
#     echo "⚠️  $CLM_FILE already exists. Skipping Step 3."
# fi

# # ------------------------------------------------------------------------------
# # 📦 Step 4: Partition Contigs into K Groups
# # ------------------------------------------------------------------------------
# echo -e "\n📦 [Step 4] Partitioning contigs into $CHROMOSOME_K groups..."

# # 检查是否已经有 partition 结果
# if [ ! -f "groups.txt" ]; then
#     echo "   Running allhic partition..."
#     allhic partition --clm $CLM_FILE \
#                      --fasta draft.asm.fasta \
#                      -k $CHROMOSOME_K
    
#     echo "✅ Step 4 Done. Check groups.txt and counts_*.txt files."
# else
#     echo "⚠️  groups.txt already exists. Skipping Step 4."
# fi

# # ------------------------------------------------------------------------------
# # ⚙️ Step 5: Optimize Each Group (Ordering & Orientation)
# # ------------------------------------------------------------------------------
# echo -e "\n⚙️  [Step 5] Optimizing ordering and orientation for each group..."

# # 检查是否需要运行 optimize
# NEED_OPTIMIZE=false
# for i in $(seq 1 "$CHROMOSOME_K"); do
#     TOUR_FILE="group${i}.tour"
#     if [ ! -f "$TOUR_FILE" ]; then
#         NEED_OPTIMIZE=true
#         break
#     fi
# done

# if [ "$NEED_OPTIMIZE" = true ]; then
#     # 生成优化命令列表
#     echo "   Generating optimization commands..."
#     > cmd.list
    
#     for i in $(seq 1 "$CHROMOSOME_K"); do
#         COUNT_FILE="counts_${RE_MOTIF}.${CHROMOSOME_K}g${i}.txt"
        
#         if [ -f "$COUNT_FILE" ]; then
#             echo "allhic optimize --clm $CLM_FILE --counts $COUNT_FILE --output group${i}" >> cmd.list
#         else
#             echo "⚠️  Warning: $COUNT_FILE not found for group $i"
#         fi
#     done
    
#     # 并行运行优化
#     echo "   Running optimization (parallel with $THREADS threads)..."
#     if command -v ParaFly &> /dev/null; then
#         ParaFly -c cmd.list -CPU "$THREADS" -failed_cmds cmd.list.failed
#     else
#         echo "   ParaFly not found, running sequentially..."
#         while read -r cmd; do
#             eval "$cmd"
#         done < cmd.list
#     fi
    
#     echo "✅ Step 5 Done. Check group*.tour files."
# else
#     echo "⚠️  All .tour files exist. Skipping Step 5."
# fi

# # ------------------------------------------------------------------------------
# # 🏗️ Step 6: Build Final Chromosome-Scale Assembly
# # ------------------------------------------------------------------------------
# echo -e "\n🏗️  [Step 6] Building chromosome-scale scaffolds..."

# if [ ! -f "groups.asm.fasta" ]; then
#     echo "   Running allhic build..."
    
#     # 收集所有 tour 文件
#     TOUR_FILES=""
#     for i in $(seq 1 "$CHROMOSOME_K"); do
#         if [ -f "group${i}.tour" ]; then
#             TOUR_FILES="$TOUR_FILES group${i}.tour"
#         fi
#     done
    
#     if [ -z "$TOUR_FILES" ]; then
#         echo "❌ Error: No .tour files found!"
#         exit 1
#     fi
    
#     allhic build --fasta draft.asm.fasta \
#                  --tours $TOUR_FILES
    
#     echo "✅ Step 6 Done."
# else
#     echo "⚠️  groups.asm.fasta already exists. Skipping Step 6."
# fi

# # ------------------------------------------------------------------------------
# # 📊 Step 7: Assessment (Optional)
# # ------------------------------------------------------------------------------
# echo -e "\n📊 [Step 7] Assessing assembly quality (optional)..."

# if [ -f "groups.asm.fasta" ]; then
#     echo "   Running allhic assess..."
#     allhic assess --fasta groups.asm.fasta \
#                   --bam $PRUNED_BAM \
#                   --output assessment
#     echo "✅ Step 7 Done. Check assessment files."
# fi

# # ==============================================================================
# # 🎉 完成
# # ==============================================================================
# echo -e "\n🎉 ================= ALLHiC PIPELINE FINISHED ================= 🎉"
# echo "📊 主要输出文件:"
# echo "   1. groups.txt          - Contig 分组信息"
# echo "   2. groups.agp          - AGP 格式 scaffold 信息"
# echo "   3. groups.asm.fasta    - 最终染色体序列"
# echo "   4. group*.tour         - 每个组的排序结果"
# echo "   5. assessment/         - 质量评估结果"
# echo -e "===============================================================\n"

# # 清理临时文件 (可选)
# read -p "是否删除中间文件以节省空间? (y/N): " -n 1 -r
# echo
# if [[ $REPLY =~ ^[Yy]$ ]]; then
#     echo "🧹 Cleaning up intermediate files..."
#     rm -f ${SAMPLE_NAME}.filtered.bam
#     rm -f cmd.list cmd.list.completed cmd.list.failed
#     echo "✅ Cleanup done."
# fi

# echo "🎊 All done! Check your results in $WORK_DIR"

#!/bin/bash

# ==============================================================================
# 🧬 ALLHiC Pipeline Automation Script (Corrected for Actual Commands)
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
RE_MOTIF="GATC"                 # 酶切序列 (MboI = GATC)

# 4. ⚙️ 系统资源
THREADS=64                      # CPU 线程数

# 5. 📛 样本名称
SAMPLE_NAME="sample"            # 输出文件前缀

# 6. 🔧 二倍体处理 (可选)
USE_ALLELE_PRUNING=false        # 如果是二倍体且有 allele table, 设为 true
ALLELE_TABLE=""                 # alleles.table 文件路径 (如果使用)

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
echo "✂️  Enzyme   : MboI ($RE_MOTIF)"
echo "🔢 Target K : $CHROMOSOME_K"
echo "🧵 Threads  : $THREADS"
echo -e "=============================================================\n"

# 创建并进入工作目录
mkdir -p "$WORK_DIR"
cd "$WORK_DIR"

# 🔗 Step 0: 建立软链接
echo -e "🔗 [Step 0] Linking input files..."
ln -sf "$REF_FA_RAW" draft.asm.fasta
ln -sf "$R1_FQ_RAW" Lib_R1.fastq.gz
ln -sf "$R2_FQ_RAW" Lib_R2.fastq.gz
echo "✅ Step 0 Done."

# ------------------------------------------------------------------------------
# 🔍 Step 1: Build Index and Map Hi-C Reads
# ------------------------------------------------------------------------------
echo -e "\n🔍 [Step 1] Building index and mapping Hi-C reads..."

if [ ! -f "${SAMPLE_NAME}.bam" ]; then
    # 1.1 建索引
    if [ ! -f "draft.asm.fasta.bwt" ]; then
        echo "   Building BWA index..."
        bwa index draft.asm.fasta
    fi
    
    if [ ! -f "draft.asm.fasta.fai" ]; then
        echo "   Building samtools index..."
        samtools faidx draft.asm.fasta
    fi

    # 1.2 比对 Hi-C reads
    echo "   Mapping Hi-C reads (this may take a while)..."
    bwa mem -5SPM -t "$THREADS" draft.asm.fasta \
        Lib_R1.fastq.gz Lib_R2.fastq.gz \
        | samtools view -@ "$THREADS" -bS - \
        | samtools sort -@ "$THREADS" -o ${SAMPLE_NAME}.bam -

    # 1.3 索引 BAM
    echo "   Indexing BAM file..."
    samtools index -@ "$THREADS" ${SAMPLE_NAME}.bam
    
    echo "✅ Step 1 Done. BAM file: ${SAMPLE_NAME}.bam"
else
    echo "⚠️  ${SAMPLE_NAME}.bam already exists. Skipping Step 1."
fi

# ------------------------------------------------------------------------------
# 🔪 Step 2: Filter BAM (MAPQ filtering)
# ------------------------------------------------------------------------------
echo -e "\n🔪 [Step 2] Filtering BAM file (MAPQ >= 1)..."

FILTERED_BAM="${SAMPLE_NAME}.filtered.bam"

if [ ! -f "$FILTERED_BAM" ]; then
    echo "   Filtering low-quality alignments..."
    samtools view -@ "$THREADS" -bq 1 ${SAMPLE_NAME}.bam > $FILTERED_BAM
    samtools index -@ "$THREADS" $FILTERED_BAM
    
    echo "✅ Step 2 Done. Filtered BAM: $FILTERED_BAM"
else
    echo "⚠️  $FILTERED_BAM already exists. Skipping Step 2."
fi

# ------------------------------------------------------------------------------
# 🧬 Step 3: Extract Hi-C Link Information
# ------------------------------------------------------------------------------
echo -e "\n🧬 [Step 3] Extracting Hi-C link information..."

# 根据实际命令格式: allhic extract bamfile fastafile [flags]
COUNTS_FILE="counts_${RE_MOTIF}.txt"
PAIRS_FILE="pairs.txt"

if [ ! -f "$COUNTS_FILE" ] || [ ! -f "$PAIRS_FILE" ]; then
    echo "   Running allhic extract..."
    allhic extract $FILTERED_BAM draft.asm.fasta --RE "$RE_MOTIF"
    
    echo "✅ Step 3 Done. Output: $COUNTS_FILE and $PAIRS_FILE"
else
    echo "⚠️  $COUNTS_FILE and $PAIRS_FILE already exist. Skipping Step 3."
fi

# ------------------------------------------------------------------------------
# 🔪 Step 3b: Prune Allelic Links (Optional, for Diploids)
# ------------------------------------------------------------------------------
if [ "$USE_ALLELE_PRUNING" = true ] && [ -n "$ALLELE_TABLE" ]; then
    echo -e "\n🔪 [Step 3b] Pruning allelic links..."
    
    PRUNED_PAIRS="pairs.pruned.txt"
    
    if [ ! -f "$PRUNED_PAIRS" ]; then
        echo "   Running allhic prune..."
        # 格式: allhic prune alleles.table pairs.txt
        allhic prune "$ALLELE_TABLE" "$PAIRS_FILE"
        
        # prune 会生成新的 pairs 文件，需要重命名
        if [ -f "pairs.pruned.txt" ]; then
            PAIRS_FILE="pairs.pruned.txt"
        fi
        
        echo "✅ Step 3b Done. Pruned pairs: $PAIRS_FILE"
    else
        echo "⚠️  $PRUNED_PAIRS already exists. Skipping Step 3b."
        PAIRS_FILE="pairs.pruned.txt"
    fi
fi

# ------------------------------------------------------------------------------
# 📦 Step 4: Partition Contigs into K Groups
# ------------------------------------------------------------------------------
echo -e "\n📦 [Step 4] Partitioning contigs into $CHROMOSOME_K groups..."

# 检查是否已经有 partition 结果
CLUSTERS_FILE="clusters.txt"

if [ ! -f "$CLUSTERS_FILE" ]; then
    echo "   Running allhic partition..."
    # 格式: allhic partition counts_RE.txt pairs.txt k
    allhic partition "$COUNTS_FILE" "$PAIRS_FILE" "$CHROMOSOME_K"
    
    echo "✅ Step 4 Done. Check $CLUSTERS_FILE and counts_${RE_MOTIF}.${CHROMOSOME_K}g*.txt files."
else
    echo "⚠️  $CLUSTERS_FILE already exists. Skipping Step 4."
fi

# ------------------------------------------------------------------------------
# ⚙️ Step 5: Optimize Each Group (Ordering & Orientation)
# ------------------------------------------------------------------------------
echo -e "\n⚙️  [Step 5] Optimizing ordering and orientation for each group..."

# 生成 CLM 文件名 (根据 extract 输出)
CLM_FILE="${SAMPLE_NAME}.filtered.clm"
if [ ! -f "$CLM_FILE" ]; then
    # 尝试其他可能的命名
    CLM_FILE=$(ls *.clm 2>/dev/null | head -1)
    if [ -z "$CLM_FILE" ]; then
        echo "❌ Error: .clm file not found! Check extract step output."
        exit 1
    fi
fi

# 检查是否需要运行 optimize
NEED_OPTIMIZE=false
for i in $(seq 1 "$CHROMOSOME_K"); do
    TOUR_FILE="group${i}.tour"
    if [ ! -f "$TOUR_FILE" ]; then
        NEED_OPTIMIZE=true
        break
    fi
done

if [ "$NEED_OPTIMIZE" = true ]; then
    # 生成优化命令列表
    echo "   Generating optimization commands..."
    > cmd.list
    
    for i in $(seq 1 "$CHROMOSOME_K"); do
        COUNT_FILE="counts_${RE_MOTIF}.${CHROMOSOME_K}g${i}.txt"
        
        if [ -f "$COUNT_FILE" ]; then
            # 格式: allhic optimize counts_RE.txt clmfile
            echo "allhic optimize $COUNT_FILE $CLM_FILE" >> cmd.list
        else
            echo "⚠️  Warning: $COUNT_FILE not found for group $i"
        fi
    done
    
    # 并行运行优化
    echo "   Running optimization (parallel with $THREADS threads)..."
    if command -v ParaFly &> /dev/null; then
        ParaFly -c cmd.list -CPU "$THREADS" -failed_cmds cmd.list.failed
    else
        echo "   ParaFly not found, running sequentially..."
        while read -r cmd; do
            echo "   Executing: $cmd"
            eval "$cmd"
        done < cmd.list
    fi
    
    echo "✅ Step 5 Done. Check group*.tour files."
else
    echo "⚠️  All .tour files exist. Skipping Step 5."
fi

# ------------------------------------------------------------------------------
# 🏗️ Step 6: Build Final Chromosome-Scale Assembly
# ------------------------------------------------------------------------------
echo -e "\n🏗️  [Step 6] Building chromosome-scale scaffolds..."

OUTPUT_FASTA="groups.asm.fasta"

if [ ! -f "$OUTPUT_FASTA" ]; then
    echo "   Running allhic build..."
    
    # 收集所有 tour 文件
    TOUR_FILES=""
    for i in $(seq 1 "$CHROMOSOME_K"); do
        if [ -f "group${i}.tour" ]; then
            TOUR_FILES="$TOUR_FILES group${i}.tour"
        fi
    done
    
    if [ -z "$TOUR_FILES" ]; then
        echo "❌ Error: No .tour files found!"
        exit 1
    fi
    
    # 格式: allhic build tourfile1 tourfile2 ... contigs.fasta asm.chr.fasta
    allhic build $TOUR_FILES draft.asm.fasta $OUTPUT_FASTA
    
    echo "✅ Step 6 Done. Output: $OUTPUT_FASTA and groups.agp"
else
    echo "⚠️  $OUTPUT_FASTA already exists. Skipping Step 6."
fi

# ------------------------------------------------------------------------------
# 📊 Step 7: Assessment (Optional)
# ------------------------------------------------------------------------------
echo -e "\n📊 [Step 7] Assessing assembly quality (optional)..."

if [ -f "$OUTPUT_FASTA" ] && [ -f "groups.agp" ]; then
    echo "   Running allhic assess for each chromosome..."
    
    # 从 AGP 文件提取染色体列表
    CHR_LIST=$(awk '{print $1}' groups.agp | sort -u)
    
    for chr in $CHR_LIST; do
        if [ ! -f "assess_${chr}.txt" ]; then
            echo "   Assessing $chr..."
            # 格式: allhic assess bamfile bedfile chr1
            # 注意: 需要先从 agp 生成 bed 文件
            awk -v chr="$chr" '$1==chr {print $1"\t"$2"\t"$3}' groups.agp > ${chr}.bed
            
            if [ -f "${chr}.bed" ]; then
                allhic assess $FILTERED_BAM ${chr}.bed $chr > assess_${chr}.txt 2>&1 || true
            fi
        fi
    done
    
    echo "✅ Step 7 Done. Check assess_*.txt files."
fi

# ==============================================================================
# 🎉 完成
# ==============================================================================
echo -e "\n🎉 ================= ALLHiC PIPELINE FINISHED ================= 🎉"
echo "📊 主要输出文件:"
echo "   1. clusters.txt        - Contig 分组信息"
echo "   2. groups.agp          - AGP 格式 scaffold 信息"
echo "   3. $OUTPUT_FASTA       - 最终染色体序列"
echo "   4. group*.tour         - 每个组的排序结果"
echo "   5. counts_${RE_MOTIF}.txt - Hi-C link counts"
echo "   6. pairs.txt           - Hi-C pairs"
echo "   7. assess_*.txt        - 质量评估结果 (如果运行)"
echo -e "===============================================================\n"

# 生成统计报告
if [ -f "$OUTPUT_FASTA" ]; then
    echo "📈 Assembly Statistics:"
    echo "   Total scaffolds: $(grep -c '^>' $OUTPUT_FASTA)"
    echo "   Total length: $(awk '/^>/ {next} {sum+=length($0)} END {print sum}' $OUTPUT_FASTA) bp"
fi

# 清理临时文件 (可选)
read -p "是否删除中间文件以节省空间? (y/N): " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    echo "🧹 Cleaning up intermediate files..."
    rm -f ${SAMPLE_NAME}.filtered.bam ${SAMPLE_NAME}.filtered.bam.bai
    rm -f cmd.list cmd.list.completed cmd.list.failed
    rm -f *.bed
    echo "✅ Cleanup done."
fi

echo "🎊 All done! Check your results in $WORK_DIR"