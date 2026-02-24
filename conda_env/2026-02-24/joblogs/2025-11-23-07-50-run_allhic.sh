#!/bin/bash

# ==============================================================================
# 🧬 ALLHiC Pipeline Automation Script (Debug Version)
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
SAMPLE_NAME="zhugecai"            # 输出文件前缀

# 6. 🔧 二倍体处理 (可选)
USE_ALLELE_PRUNING=false        # 如果是二倍体且有 allele table, 设为 true
ALLELE_TABLE=""                 # alleles.table 文件路径 (如果使用)

# 7. 📝 日志设置
LOG_DIR="logs"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
LOG_FILE="${LOG_DIR}/allhic_${TIMESTAMP}.log"
CMD_LOG="${LOG_DIR}/commands_${TIMESTAMP}.sh"

# ==============================================================================

# 🛑 错误处理函数
set -e
set -o pipefail

handle_error() {
    local line_num=$1
    echo -e "\n❌ [ERROR] 脚本在第 $line_num 行发生错误!" | tee -a "$LOG_FILE"
    echo "   查看日志: $LOG_FILE" | tee -a "$LOG_FILE"
    echo "   查看命令: $CMD_LOG" | tee -a "$LOG_FILE"
    exit 1
}
trap 'handle_error $LINENO' ERR

# 📝 日志函数
log_info() {
    echo -e "[$(date +'%Y-%m-%d %H:%M:%S')] [INFO] $*" | tee -a "$LOG_FILE"
}

log_warn() {
    echo -e "[$(date +'%Y-%m-%d %H:%M:%S')] [WARN] $*" | tee -a "$LOG_FILE"
}

log_error() {
    echo -e "[$(date +'%Y-%m-%d %H:%M:%S')] [ERROR] $*" | tee -a "$LOG_FILE"
}

log_cmd() {
    echo -e "\n# $(date +'%Y-%m-%d %H:%M:%S')" >> "$CMD_LOG"
    echo "# $1" >> "$CMD_LOG"
    echo "$2" >> "$CMD_LOG"
    echo "" >> "$CMD_LOG"
}

# 🎯 执行命令并记录
run_cmd() {
    local description=$1
    shift
    local cmd="$@"
    
    log_info "▶️  Executing: $description"
    log_info "   Command: $cmd"
    log_cmd "$description" "$cmd"
    
    echo "----------------------------------------" | tee -a "$LOG_FILE"
    echo "COMMAND: $cmd" | tee -a "$LOG_FILE"
    echo "----------------------------------------" | tee -a "$LOG_FILE"
    
    # 执行命令并捕获输出
    if eval "$cmd" 2>&1 | tee -a "$LOG_FILE"; then
        log_info "✅ Success: $description"
        return 0
    else
        log_error "❌ Failed: $description"
        return 1
    fi
}

# 🚀 开始流程
mkdir -p "$LOG_DIR"

log_info "=========================================="
log_info "🚀 STARTING ALLHiC PIPELINE"
log_info "=========================================="
log_info "📂 Work Dir : $WORK_DIR"
log_info "🧬 Ref Fasta: $REF_FA_RAW"
log_info "✂️  Enzyme   : MboI ($RE_MOTIF)"
log_info "🔢 Target K : $CHROMOSOME_K"
log_info "🧵 Threads  : $THREADS"
log_info "📝 Log File : $LOG_FILE"
log_info "📝 Cmd Log  : $CMD_LOG"
log_info "=========================================="

# 创建并进入工作目录
mkdir -p "$WORK_DIR"
cd "$WORK_DIR"

# 🔗 Step 0: 建立软链接
log_info "\n🔗 [Step 0] Linking input files..."

if [ ! -L "draft.asm.fasta" ]; then
    run_cmd "Link reference fasta" "ln -sf '$REF_FA_RAW' draft.asm.fasta"
fi

if [ ! -L "Lib_R1.fastq.gz" ]; then
    run_cmd "Link R1 reads" "ln -sf '$R1_FQ_RAW' Lib_R1.fastq.gz"
fi

if [ ! -L "Lib_R2.fastq.gz" ]; then
    run_cmd "Link R2 reads" "ln -sf '$R2_FQ_RAW' Lib_R2.fastq.gz"
fi

log_info "✅ Step 0 Done."

# ------------------------------------------------------------------------------
# 🔍 Step 1: Build Index and Map Hi-C Reads
# ------------------------------------------------------------------------------
log_info "\n🔍 [Step 1] Building index and mapping Hi-C reads..."

if [ ! -f "${SAMPLE_NAME}.bam" ]; then
    # 1.1 建索引
    if [ ! -f "draft.asm.fasta.bwt" ]; then
        run_cmd "Build BWA index" "bwa index draft.asm.fasta"
    else
        log_info "   BWA index already exists, skipping..."
    fi
    
    if [ ! -f "draft.asm.fasta.fai" ]; then
        run_cmd "Build samtools index" "samtools faidx draft.asm.fasta"
    else
        log_info "   Samtools fai index already exists, skipping..."
    fi

    # 1.2 比对 Hi-C reads
    log_info "   Mapping Hi-C reads (this may take a while)..."
    run_cmd "Map Hi-C reads with BWA-MEM" \
        "bwa mem -5SPM -t $THREADS draft.asm.fasta Lib_R1.fastq.gz Lib_R2.fastq.gz | \
         samtools view -@ $THREADS -bS - | \
         samtools sort -@ $THREADS -o ${SAMPLE_NAME}.bam -"

    # 1.3 索引 BAM
    run_cmd "Index BAM file" "samtools index -@ $THREADS ${SAMPLE_NAME}.bam"
    
    log_info "✅ Step 1 Done. BAM file: ${SAMPLE_NAME}.bam"
else
    log_info "⚠️  ${SAMPLE_NAME}.bam already exists. Skipping Step 1."
fi

# 统计 BAM 信息
log_info "\n📊 BAM Statistics:"
run_cmd "BAM flagstat" "samtools flagstat ${SAMPLE_NAME}.bam | head -10"

# ------------------------------------------------------------------------------
# 🔪 Step 2: Filter BAM (MAPQ filtering)
# ------------------------------------------------------------------------------
log_info "\n🔪 [Step 2] Filtering BAM file (MAPQ >= 1)..."

FILTERED_BAM="${SAMPLE_NAME}.filtered.bam"

if [ ! -f "$FILTERED_BAM" ]; then
    run_cmd "Filter BAM by MAPQ" "samtools view -@ $THREADS -bq 1 ${SAMPLE_NAME}.bam > $FILTERED_BAM"
    run_cmd "Index filtered BAM" "samtools index -@ $THREADS $FILTERED_BAM"
    
    log_info "✅ Step 2 Done. Filtered BAM: $FILTERED_BAM"
else
    log_info "⚠️  $FILTERED_BAM already exists. Skipping Step 2."
fi

# 统计过滤结果
log_info "\n📊 Filter Statistics:"
TOTAL_READS=$(samtools view -c ${SAMPLE_NAME}.bam)
FILTERED_READS=$(samtools view -c $FILTERED_BAM)
RETENTION=$(echo "scale=2; $FILTERED_READS * 100 / $TOTAL_READS" | bc)
log_info "   Total reads: $TOTAL_READS"
log_info "   Filtered reads (MAPQ>=1): $FILTERED_READS"
log_info "   Retention rate: ${RETENTION}%"

# ------------------------------------------------------------------------------
# 🧬 Step 3: Extract Hi-C Link Information
# ------------------------------------------------------------------------------
log_info "\n🧬 [Step 3] Extracting Hi-C link information..."

BASE_NAME=$(basename $FILTERED_BAM .bam)
COUNTS_FILE="${BASE_NAME}.counts_${RE_MOTIF}.txt"
PAIRS_FILE="${BASE_NAME}.pairs.txt"
CLM_FILE="${BASE_NAME}.clm"

if [ ! -f "$COUNTS_FILE" ] || [ ! -f "$CLM_FILE" ]; then
    run_cmd "Extract Hi-C links" \
        "allhic extract $FILTERED_BAM draft.asm.fasta --RE '$RE_MOTIF' --minLinks 3"
    
    log_info "\n✅ Step 3 Done. Output files:"
    log_info "      - $COUNTS_FILE ($(wc -l < $COUNTS_FILE 2>/dev/null || echo 0) lines)"
    log_info "      - $PAIRS_FILE ($(wc -l < $PAIRS_FILE 2>/dev/null || echo 0) lines)"
    log_info "      - $CLM_FILE ($(wc -l < $CLM_FILE 2>/dev/null || echo 0) lines)"
    
    # 检查 pairs 文件
    if [ ! -s "$PAIRS_FILE" ]; then
        log_warn "⚠️  WARNING: $PAIRS_FILE is empty!"
        log_warn "   This is a known issue. Will generate it from CLM in next step."
    fi
else
    log_info "⚠️  Extract output files already exist. Skipping Step 3."
fi

# 验证关键文件
if [ ! -s "$CLM_FILE" ]; then
    log_error "❌ CRITICAL: $CLM_FILE is empty or missing!"
    exit 1
fi

# 显示 CLM 统计
log_info "\n📊 CLM File Statistics:"
CLM_SIZE=$(stat -f%z "$CLM_FILE" 2>/dev/null || stat -c%s "$CLM_FILE" 2>/dev/null)
CLM_LINES=$(wc -l < "$CLM_FILE")
log_info "   Size: $(numfmt --to=iec-i --suffix=B $CLM_SIZE 2>/dev/null || echo ${CLM_SIZE} bytes)"
log_info "   Lines: $CLM_LINES"
run_cmd "Show CLM header" "head -3 $CLM_FILE"

# ------------------------------------------------------------------------------
# 🔪 Step 3b: Generate proper pairs.txt if needed
# ------------------------------------------------------------------------------
if [ ! -s "$PAIRS_FILE" ]; then
    log_info "\n🔧 [Step 3b] Generating pairs.txt from CLM file..."
    
    cat > generate_pairs.py << 'PYEOF'
import sys

# 读取 counts 文件获取 RE 信息
re_counts = {}
counts_file = sys.argv[1]
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
clm_file = sys.argv[2]
pairs_file = sys.argv[3]

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

print(f"Generated {pairs_file}: {sum(1 for _ in open(pairs_file))} lines")
PYEOF

    run_cmd "Generate pairs.txt from CLM" \
        "python3 generate_pairs.py '$COUNTS_FILE' '$CLM_FILE' '$PAIRS_FILE'"
    
    log_info "✅ Generated $PAIRS_FILE"
fi

# ------------------------------------------------------------------------------
# 📦 Step 4: Partition Contigs into K Groups
# ------------------------------------------------------------------------------
log_info "\n📦 [Step 4] Partitioning contigs into $CHROMOSOME_K groups..."

PARTITION_CHECK="${BASE_NAME}.counts_${RE_MOTIF}.${CHROMOSOME_K}g1.txt"
CLUSTERS_FILE="${BASE_NAME}.clusters.txt"

if [ ! -f "$PARTITION_CHECK" ]; then
    log_info "   Attempting partition with standard parameters..."
    
    run_cmd "Partition contigs (attempt 1)" \
        "allhic partition '$COUNTS_FILE' '$PAIRS_FILE' $CHROMOSOME_K"
    
    # 检查结果
    if [ -f "$CLUSTERS_FILE" ]; then
        CLUSTER_COUNT=$(grep -v "^#" "$CLUSTERS_FILE" | wc -l)
        log_info "   Generated $CLUSTER_COUNT clusters"
        
        if [ "$CLUSTER_COUNT" -eq 0 ]; then
            log_warn "   ⚠️  No clusters generated! Trying relaxed parameters..."
            
            run_cmd "Partition contigs (attempt 2 - relaxed)" \
                "allhic partition '$COUNTS_FILE' '$PAIRS_FILE' $CHROMOSOME_K --minREs 5 --maxLinkDensity 5 --nonInformativeRatio 5"
            
            CLUSTER_COUNT=$(grep -v "^#" "$CLUSTERS_FILE" | wc -l)
            log_info "   Generated $CLUSTER_COUNT clusters (second attempt)"
            
            if [ "$CLUSTER_COUNT" -eq 0 ]; then
                log_warn "   ⚠️  Still no clusters! Trying very relaxed parameters..."
                
                run_cmd "Partition contigs (attempt 3 - very relaxed)" \
                    "allhic partition '$COUNTS_FILE' '$PAIRS_FILE' $CHROMOSOME_K --minREs 3 --maxLinkDensity 10 --nonInformativeRatio 10"
                
                CLUSTER_COUNT=$(grep -v "^#" "$CLUSTERS_FILE" | wc -l)
                log_info "   Generated $CLUSTER_COUNT clusters (third attempt)"
                
                if [ "$CLUSTER_COUNT" -eq 0 ]; then
                    log_error "   ❌ Still no clusters after 3 attempts!"
                    log_error "   Possible issues:"
                    log_error "   1. Hi-C signal too weak"
                    log_error "   2. K value too high (try smaller k, e.g., 6 or 8)"
                    log_error "   3. Wrong enzyme motif"
                    log_error "   4. Poor Hi-C library quality"
                    
                    log_info "\n   Suggestions:"
                    log_info "   - Try smaller k values: 6, 8, 10"
                    log_info "   - Check Hi-C library QC"
                    log_info "   - Verify enzyme is correct ($RE_MOTIF)"
                    log_info "   - Try with unfiltered BAM (${SAMPLE_NAME}.bam)"
                    
                    exit 1
                fi
            fi
        fi
    fi
    
    # 列出生成的文件
    log_info "\n   Partition output files:"
    run_cmd "List partition files" "ls -lh ${BASE_NAME}.counts_${RE_MOTIF}.${CHROMOSOME_K}g*.txt 2>/dev/null || echo 'No group files found'"
    
    log_info "✅ Step 4 Done."
else
    log_info "⚠️  Partition results already exist. Skipping Step 4."
fi

# 显示 clusters 内容
if [ -f "$CLUSTERS_FILE" ]; then
    log_info "\n📊 Clusters Summary:"
    run_cmd "Show clusters file" "cat '$CLUSTERS_FILE'"
fi

# ------------------------------------------------------------------------------
# ⚙️ Step 5: Optimize Each Group (Ordering & Orientation)
# ------------------------------------------------------------------------------
log_info "\n⚙️  [Step 5] Optimizing ordering and orientation for each group..."

# 检查是否有有效的 clusters
if [ ! -f "$CLUSTERS_FILE" ]; then
    log_error "❌ clusters.txt not found! Cannot proceed."
    exit 1
fi

CLUSTER_COUNT=$(grep -v "^#" "$CLUSTERS_FILE" | wc -l)
if [ "$CLUSTER_COUNT" -eq 0 ]; then
    log_error "❌ No valid clusters found! Cannot proceed to optimize."
    exit 1
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
    log_info "   Detecting partition output file pattern..."
    
    # 智能查找文件模式
    PATTERN1="${BASE_NAME}.counts_${RE_MOTIF}.${CHROMOSOME_K}g"
    PATTERN2="counts_${RE_MOTIF}.${CHROMOSOME_K}g"
    
    if ls ${PATTERN1}1.txt &>/dev/null; then
        FILE_PATTERN="$PATTERN1"
        log_info "   Using pattern: ${FILE_PATTERN}*.txt"
    elif ls ${PATTERN2}1.txt &>/dev/null; then
        FILE_PATTERN="$PATTERN2"
        log_info "   Using pattern: ${FILE_PATTERN}*.txt"
    else
        log_error "   ❌ Cannot find partition output files!"
        log_error "   Searched for:"
        log_error "     - ${PATTERN1}*.txt"
        log_error "     - ${PATTERN2}*.txt"
        run_cmd "List available txt files" "ls -lh *.txt"
        exit 1
    fi
    
    # 生成优化命令列表
    log_info "   Generating optimization commands..."
    > cmd.list
    
    VALID_CMDS=0
    for i in $(seq 1 "$CHROMOSOME_K"); do
        COUNT_FILE="${FILE_PATTERN}${i}.txt"
        
        if [ -f "$COUNT_FILE" ]; then
            echo "allhic optimize $COUNT_FILE $CLM_FILE" >> cmd.list
            VALID_CMDS=$((VALID_CMDS + 1))
            log_info "     Group $i: $COUNT_FILE"
        else
            log_warn "     ⚠️  Group $i: $COUNT_FILE not found"
        fi
    done
    
    if [ "$VALID_CMDS" -eq 0 ]; then
        log_error "❌ No valid optimization commands generated!"
        exit 1
    fi
    
    log_info "   Generated $VALID_CMDS optimization commands."
    log_info "   Commands saved to: cmd.list"
    
    # 显示命令列表
    run_cmd "Show optimization commands" "cat cmd.list"
    
    # 并行运行优化
    log_info "\n   Running optimization (parallel with $THREADS threads)..."
    if command -v ParaFly &> /dev/null; then
        run_cmd "Run optimization with ParaFly" \
            "ParaFly -c cmd.list -CPU $THREADS -failed_cmds cmd.list.failed"
    else
        log_info "   ParaFly not found, running sequentially..."
        while read -r cmd; do
            run_cmd "Optimize group" "$cmd"
        done < cmd.list
    fi
    
    log_info "✅ Step 5 Done. Check group*.tour files."
else
    log_info "⚠️  All .tour files exist. Skipping Step 5."
fi

# 列出生成的 tour 文件
log_info "\n📊 Generated tour files:"
run_cmd "List tour files" "ls -lh group*.tour 2>/dev/null || echo 'No tour files found'"

# ------------------------------------------------------------------------------
# 🏗️ Step 6: Build Final Chromosome-Scale Assembly
# ------------------------------------------------------------------------------
log_info "\n🏗️  [Step 6] Building chromosome-scale scaffolds..."

OUTPUT_FASTA="groups.asm.fasta"

if [ ! -f "$OUTPUT_FASTA" ]; then
    # 收集所有 tour 文件
    TOUR_FILES=""
    TOUR_COUNT=0
    for i in $(seq 1 "$CHROMOSOME_K"); do
        if [ -f "group${i}.tour" ]; then
            TOUR_FILES="$TOUR_FILES group${i}.tour"
            TOUR_COUNT=$((TOUR_COUNT + 1))
        fi
    done
    
    if [ -z "$TOUR_FILES" ]; then
        log_error "❌ No .tour files found! Cannot build assembly."
        exit 1
    fi
    
    log_info "   Found $TOUR_COUNT tour files"
    log_info "   Tour files: $TOUR_FILES"
    
    run_cmd "Build chromosome-scale assembly" \
        "allhic build $TOUR_FILES draft.asm.fasta $OUTPUT_FASTA"
    
    log_info "✅ Step 6 Done. Output: $OUTPUT_FASTA and groups.agp"
else
    log_info "⚠️  $OUTPUT_FASTA already exists. Skipping Step 6."
fi

# 统计最终结果
if [ -f "$OUTPUT_FASTA" ]; then
    log_info "\n📊 Final Assembly Statistics:"
    SCAFFOLD_COUNT=$(grep -c "^>" "$OUTPUT_FASTA")
    TOTAL_LENGTH=$(awk '/^>/ {next} {sum+=length($0)} END {print sum}' "$OUTPUT_FASTA")
    log_info "   Total scaffolds: $SCAFFOLD_COUNT"
    log_info "   Total length: $TOTAL_LENGTH bp"
fi

# ------------------------------------------------------------------------------
# 📊 Step 7: Assessment (Optional)
# ------------------------------------------------------------------------------
log_info "\n📊 [Step 7] Assessing assembly quality (optional)..."

if [ -f "$OUTPUT_FASTA" ] && [ -f "groups.agp" ]; then
    log_info "   Extracting chromosome list from AGP..."
    CHR_LIST=$(awk '{print $1}' groups.agp | sort -u)
    
    for chr in $CHR_LIST; do
        if [ ! -f "assess_${chr}.txt" ]; then
            log_info "   Assessing $chr..."
            
            # 从 agp 生成 bed 文件
            awk -v chr="$chr" '$1==chr {print $1"\t"$2"\t"$3}' groups.agp > ${chr}.bed
            
            if [ -f "${chr}.bed" ]; then
                run_cmd "Assess $chr" \
                    "allhic assess $FILTERED_BAM ${chr}.bed $chr > assess_${chr}.txt 2>&1 || true"
            fi
        fi
    done
    
    log_info "✅ Step 7 Done."
fi

# ==============================================================================
# 🎉 完成
# ==============================================================================
log_info "\n=========================================="
log_info "🎉 ALLHiC PIPELINE FINISHED"
log_info "=========================================="
log_info "📊 主要输出文件:"
log_info "   1. ${BASE_NAME}.clusters.txt - Contig 分组信息"
log_info "   2. groups.agp          - AGP 格式 scaffold 信息"
log_info "   3. $OUTPUT_FASTA       - 最终染色体序列"
log_info "   4. group*.tour         - 每个组的排序结果"
log_info "   5. $COUNTS_FILE        - Hi-C link counts"
log_info "   6. $PAIRS_FILE         - Hi-C pairs"
log_info "   7. $CLM_FILE           - CLM matrix"
log_info "   8. assess_*.txt        - 质量评估结果"
log_info "=========================================="
log_info "📝 完整日志: $LOG_FILE"
log_info "📝 命令记录: $CMD_LOG"
log_info "=========================================="

# 清理临时文件 (可选)
read -p "是否删除中间文件以节省空间? (y/N): " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    log_info "🧹 Cleaning up intermediate files..."
    run_cmd "Remove filtered BAM" "rm -f ${SAMPLE_NAME}.filtered.bam ${SAMPLE_NAME}.filtered.bam.bai"
    run_cmd "Remove temp files" "rm -f cmd.list cmd.list.completed cmd.list.failed"
    run_cmd "Remove BED files" "rm -f *.bed"
    run_cmd "Remove Python script" "rm -f generate_pairs.py"
    log_info "✅ Cleanup done."
fi

log_info "\n🎊 All done! Check your results in $WORK_DIR"
log_info "📧 If you encounter issues, please share:"
log_info "   - Log file: $LOG_FILE"
log_info "   - Command log: $CMD_LOG"