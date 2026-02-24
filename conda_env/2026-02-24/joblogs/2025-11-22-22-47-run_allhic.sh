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
    
    # 统计过滤后的reads数量
    TOTAL_READS=$(samtools view -c ${SAMPLE_NAME}.bam)
    FILTERED_READS=$(samtools view -c $FILTERED_BAM)
    echo "   Total reads: $TOTAL_READS"
    echo "   Filtered reads (MAPQ>=1): $FILTERED_READS"
    echo "   Retention rate: $(echo "scale=2; $FILTERED_READS * 100 / $TOTAL_READS" | bc)%"
    
    echo "✅ Step 2 Done. Filtered BAM: $FILTERED_BAM"
else
    echo "⚠️  $FILTERED_BAM already exists. Skipping Step 2."
    # 仍然显示统计信息
    FILTERED_READS=$(samtools view -c $FILTERED_BAM)
    echo "   Filtered reads: $FILTERED_READS"
fi

# ------------------------------------------------------------------------------
# 🧬 Step 3: Extract Hi-C Link Information
# ------------------------------------------------------------------------------
echo -e "\n🧬 [Step 3] Extracting Hi-C link information..."

# allhic extract 会生成带 BAM 前缀的文件名
BASE_NAME=$(basename $FILTERED_BAM .bam)
COUNTS_FILE="${BASE_NAME}.counts_${RE_MOTIF}.txt"
PAIRS_FILE="${BASE_NAME}.pairs.txt"
CLM_FILE="${BASE_NAME}.clm"

if [ ! -f "$COUNTS_FILE" ] || [ ! -f "$CLM_FILE" ]; then
    echo "   Running allhic extract..."
    # 添加 --minLinks 参数，降低阈值以获取更多 pairs
    allhic extract $FILTERED_BAM draft.asm.fasta --RE "$RE_MOTIF" --minLinks 3
    
    echo "✅ Step 3 Done. Output files:"
    echo "      - $COUNTS_FILE ($(wc -l < $COUNTS_FILE 2>/dev/null || echo 0) lines)"
    echo "      - $PAIRS_FILE ($(wc -l < $PAIRS_FILE 2>/dev/null || echo 0) lines)"
    echo "      - $CLM_FILE ($(wc -l < $CLM_FILE 2>/dev/null || echo 0) lines)"
    
    # 检查 pairs 文件是否为空（这是已知问题，可以从 CLM 恢复）
    if [ ! -s "$PAIRS_FILE" ]; then
        echo "   ⚠️  WARNING: $PAIRS_FILE is empty!"
        echo "   This is a known issue with some allhic versions."
        echo "   We will generate it from CLM file in the next step."
    fi
else
    echo "⚠️  Extract output files already exist. Skipping Step 3."
fi

# 验证 CLM 文件（这是关键文件）
if [ ! -s "$CLM_FILE" ]; then
    echo ""
    echo "❌ CRITICAL ERROR: $CLM_FILE is empty or missing!"
    echo "   This file is essential for partition and optimize steps."
    echo ""
    exit 1
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
# partition 可以直接使用 CLM 文件，不需要 pairs.txt
PARTITION_CHECK="${BASE_NAME}.counts_${RE_MOTIF}.${CHROMOSOME_K}g1.txt"

if [ ! -f "$PARTITION_CHECK" ]; then
    echo "   Running allhic partition..."
    
    # 检查 pairs.txt 是否为空或格式不对
    if [ ! -s "$PAIRS_FILE" ]; then
        echo "   ⚠️  WARNING: $PAIRS_FILE is empty!"
        echo "   Generating pairs.txt from CLM file..."
        
        if [ -f "$CLM_FILE" ] && [ -s "$CLM_FILE" ]; then
            # 使用 Python 生成正确格式的 pairs.txt
            python3 << 'PYEOF'
import sys

# 读取 counts 文件获取 RE 信息
re_counts = {}
counts_file = sys.argv[1] if len(sys.argv) > 1 else 'sample.filtered.counts_GATC.txt'
with open(counts_file, 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        if len(fields) >= 2:
            contig = fields[0]
            count = int(fields[1])
            re_counts[contig] = count

# 生成 pairs.txt
clm_file = sys.argv[2] if len(sys.argv) > 2 else 'sample.filtered.clm'
pairs_file = sys.argv[3] if len(sys.argv) > 3 else 'sample.filtered.pairs.txt'

with open(pairs_file, 'w') as out:
    out.write("#X\tY\tContig1\tContig2\tRE1\tRE2\tObservedLinks\tExpectedLinksIfAdjacent\tLabel\n")
    
    with open(clm_file, 'r') as f:
        for line in f:
            fields = line.strip().split()
            if len(fields) < 3:
                continue
            
            contig1 = fields[0].rstrip('+-')
            contig2 = fields[1].rstrip('+-')
            
            re1 = re_counts.get(contig1, 0)
            re2 = re_counts.get(contig2, 0)
            observed_links = int(fields[2])
            
            out.write(f"0\t0\t{contig1}\t{contig2}\t{re1}\t{re2}\t{observed_links}\t0\t1\n")

print(f"Generated {pairs_file} successfully!")
PYEOF
            python3 -c "import sys; sys.argv = ['', '$COUNTS_FILE', '$CLM_FILE', '$PAIRS_FILE']; exec(open('/dev/stdin').read())" << 'PYEOF'
import sys

# 读取 counts 文件获取 RE 信息
re_counts = {}
with open(sys.argv[1], 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        if len(fields) >= 2:
            contig = fields[0]
            count = int(fields[1])
            re_counts[contig] = count

# 生成 pairs.txt
with open(sys.argv[3], 'w') as out:
    out.write("#X\tY\tContig1\tContig2\tRE1\tRE2\tObservedLinks\tExpectedLinksIfAdjacent\tLabel\n")
    
    with open(sys.argv[2], 'r') as f:
        for line in f:
            fields = line.strip().split()
            if len(fields) < 3:
                continue
            
            contig1 = fields[0].rstrip('+-')
            contig2 = fields[1].rstrip('+-')
            
            re1 = re_counts.get(contig1, 0)
            re2 = re_counts.get(contig2, 0)
            observed_links = int(fields[2])
            
            out.write(f"0\t0\t{contig1}\t{contig2}\t{re1}\t{re2}\t{observed_links}\t0\t1\n")

print(f"Generated pairs.txt: {sum(1 for _ in open(sys.argv[3]))} lines")
PYEOF
            
            echo "   ✓ Generated $PAIRS_FILE from CLM"
        else
            echo "   ❌ ERROR: CLM file is also missing or empty!"
            exit 1
        fi
    fi
    
    # 尝试 partition，使用更宽松的参数
    echo "   Executing partition with relaxed parameters..."
    echo "   Command: allhic partition $COUNTS_FILE $PAIRS_FILE $CHROMOSOME_K --minREs 5 --maxLinkDensity 5 --nonInformativeRatio 5"
    
    allhic partition "$COUNTS_FILE" "$PAIRS_FILE" "$CHROMOSOME_K" \
        --minREs 5 \
        --maxLinkDensity 5 \
        --nonInformativeRatio 5
    
    # 检查是否成功
    if [ ! -f "${BASE_NAME}.clusters.txt" ]; then
        echo "   ❌ ERROR: clusters.txt not generated!"
        exit 1
    fi
    
    CLUSTER_COUNT=$(grep -v "^#" ${BASE_NAME}.clusters.txt | wc -l)
    echo "   Generated $CLUSTER_COUNT clusters"
    
    if [ "$CLUSTER_COUNT" -eq 0 ]; then
        echo "   ⚠️  WARNING: No clusters generated! Trying with even more relaxed parameters..."
        
        allhic partition "$COUNTS_FILE" "$PAIRS_FILE" "$CHROMOSOME_K" \
            --minREs 3 \
            --maxLinkDensity 10 \
            --nonInformativeRatio 10
        
        CLUSTER_COUNT=$(grep -v "^#" ${BASE_NAME}.clusters.txt | wc -l)
        echo "   Generated $CLUSTER_COUNT clusters (second attempt)"
        
        if [ "$CLUSTER_COUNT" -eq 0 ]; then
            echo "   ❌ ERROR: Still no clusters! Your Hi-C data may not have strong enough signal."
            echo "   Suggestions:"
            echo "   1. Check Hi-C library quality"
            echo "   2. Try different k values (lower k might work better)"
            echo "   3. Check if enzyme is correct"
            exit 1
        fi
    fi
    
    # 列出实际生成的文件
    echo "   Generated partition files:"
    ls -lh ${BASE_NAME}.counts_${RE_MOTIF}.${CHROMOSOME_K}g*.txt 2>/dev/null || echo "   Warning: No group-specific counts files found"
    
    echo "✅ Step 4 Done."
else
    echo "⚠️  Partition results already exist. Skipping Step 4."
fi

# 验证 partition 结果
echo "   Checking partition results..."
if [ -f "${BASE_NAME}.clusters.txt" ]; then
    CLUSTER_COUNT=$(grep -v "^#" ${BASE_NAME}.clusters.txt | wc -l)
    echo "   ✓ Found $CLUSTER_COUNT clusters in ${BASE_NAME}.clusters.txt"
    
    if [ "$CLUSTER_COUNT" -eq 0 ]; then
        echo "   ❌ No valid clusters found! Cannot proceed to optimize step."
        exit 1
    fi
else
    echo "   ❌ clusters.txt not found!"
    exit 1
fi

# ------------------------------------------------------------------------------
# ⚙️ Step 5: Optimize Each Group (Ordering & Orientation)
# ------------------------------------------------------------------------------
echo -e "\n⚙️  [Step 5] Optimizing ordering and orientation for each group..."

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
    # 智能查找 partition 生成的文件
    echo "   Detecting partition output files..."
    
    # 尝试多种可能的命名模式
    PATTERN1="${BASE_NAME}.counts_${RE_MOTIF}.${CHROMOSOME_K}g"
    PATTERN2="counts_${RE_MOTIF}.${CHROMOSOME_K}g"
    
    # 检查哪种模式的文件存在
    if ls ${PATTERN1}1.txt &>/dev/null; then
        FILE_PATTERN="$PATTERN1"
        echo "   Using pattern: ${FILE_PATTERN}*.txt"
    elif ls ${PATTERN2}1.txt &>/dev/null; then
        FILE_PATTERN="$PATTERN2"
        echo "   Using pattern: ${FILE_PATTERN}*.txt"
    else
        echo "   ❌ Error: Cannot find partition output files!"
        echo "   Searched for:"
        echo "     - ${PATTERN1}*.txt"
        echo "     - ${PATTERN2}*.txt"
        echo "   Available files in directory:"
        ls -lh *.txt 2>/dev/null | head -20
        exit 1
    fi
    
    # 生成优化命令列表
    echo "   Generating optimization commands..."
    > cmd.list
    
    for i in $(seq 1 "$CHROMOSOME_K"); do
        COUNT_FILE="${FILE_PATTERN}${i}.txt"
        
        if [ -f "$COUNT_FILE" ]; then
            # 格式: allhic optimize counts_RE.txt clmfile
            echo "allhic optimize $COUNT_FILE $CLM_FILE" >> cmd.list
        else
            echo "⚠️  Warning: $COUNT_FILE not found for group $i"
        fi
    done
    
    # 检查是否有有效命令
    if [ ! -s cmd.list ]; then
        echo "❌ Error: No valid optimization commands generated!"
        exit 1
    fi
    
    echo "   Generated $(wc -l < cmd.list) optimization commands."
    
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
echo "   1. ${BASE_NAME}.counts_${RE_MOTIF}.${CHROMOSOME_K}g*.txt - Partition results"
echo "   2. groups.agp          - AGP 格式 scaffold 信息"
echo "   3. $OUTPUT_FASTA       - 最终染色体序列"
echo "   4. group*.tour         - 每个组的排序结果"
echo "   5. $COUNTS_FILE        - Hi-C link counts"
echo "   6. $PAIRS_FILE         - Hi-C pairs"
echo "   7. $CLM_FILE           - CLM matrix"
echo "   8. assess_*.txt        - 质量评估结果 (如果运行)"
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