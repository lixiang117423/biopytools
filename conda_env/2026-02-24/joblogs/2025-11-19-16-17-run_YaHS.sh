#!/bin/bash
# =============================================================================
#  🧬 YaHS 高速染色体挂载流程 - 超算提交脚本
#  适配: Hifiasm Contigs + Hi-C FASTQ 数据
# =============================================================================

# --- (A) 作业调度系统指令 ---
#CSUB -L /bin/bash
#CSUB -J YaHS_Scaffolding
#CSUB -q c01
#CSUB -o Outlog/%J.out
#CSUB -e Outlog/%J.error
#CSUB -n 64
#CSUB -R "span[hosts=1]"
#CSUB -R "rusage[mem=640000]" 

# =============================================================================
# --- (B) 💻 环境与参数配置 (请修改此处) ---
# =============================================================================
set -e 
set -o pipefail

echo "ℹ️  INFO: 作业开始于: $(date)"
echo "🖥️  INFO: 运行于计算节点: $(hostname)"

# 1. 📂 路径设置
WORK_DIR="/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/04.YaHS"
REF_FA="/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/04.YaHS/OV53_1.primary.fa"
R1_FQ="/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/04.YaHS/fastq/OV53_1-hic_R1.fastq.gz"
R2_FQ="/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/04.YaHS/fastq/OV53_1-hic_R2.fastq.gz"

# 2. ⚙️ 软件路径 & 参数
JUICER_JAR="/share/org/YZWL/yzwl_lixg/software/juicer/scripts/common/juicer_tools.jar"
ENZYME_SEQ="GATC"
THREADS=88

# 确保软件在路径中
export PATH="/share/org/YZWL/yzwl_lixg/miniforge3/envs/yahs_v.1.2.2/bin:$PATH"

# =============================================================================
# --- (C) 🚀 流程开始 ---
# =============================================================================
mkdir -p ${WORK_DIR}
cd ${WORK_DIR}

# 检查输入文件
echo "检查输入文件..."
for file in "${REF_FA}" "${R1_FQ}" "${R2_FQ}"; do
    if [ ! -f "${file}" ]; then
        echo "❌ 错误: 文件不存在 - ${file}"
        exit 1
    fi
done
echo "✅ 输入文件检查通过"

# --- 步骤 1: 建立索引 (Indexing) ---
echo "Step 1: Checking/Building Indexes..."

if [ ! -f "${REF_FA}.bwt" ]; then
    echo "Building BWA index..."
    bwa index ${REF_FA}
else
    echo "BWA index found."
fi

if [ ! -f "${REF_FA}.fai" ]; then
    echo "Building FAI index..."
    samtools faidx ${REF_FA}
else
    echo "FAI index found."
fi

# --- 步骤 2: Hi-C 比对与处理 (Mapping & Processing) ---
echo "Step 2: Mapping Hi-C reads to Contigs..."

FINAL_BAM="aligned_sorted_dedup.bam"

if [ -f "${FINAL_BAM}" ] && [ -f "${FINAL_BAM}.bai" ]; then
    echo "✅ Found existing BAM file with index, skipping mapping."
else
    echo "🚀 Running BWA mapping pipeline..."
    
    bwa mem -5SP -t ${THREADS} ${REF_FA} ${R1_FQ} ${R2_FQ} 2> bwa_mem.log | \
    samtools sort -n -@ ${THREADS} -m 4G -T tmp_nsort - 2> samtools_nsort.log | \
    samtools fixmate -m -@ ${THREADS} - - 2> samtools_fixmate.log | \
    samtools sort -@ ${THREADS} -m 4G -T tmp_sort - 2> samtools_sort.log | \
    samtools markdup -r -@ ${THREADS} - ${FINAL_BAM} 2> samtools_markdup.log
    
    samtools index -@ ${THREADS} ${FINAL_BAM}
    echo "✅ Mapping and deduplication finished."
    
    # 输出比对统计
    echo "📊 Alignment Statistics:"
    samtools flagstat -@ ${THREADS} ${FINAL_BAM} > alignment_stats.txt
    cat alignment_stats.txt
fi

# --- 步骤 3: 运行 YaHS 组装 (Scaffolding) ---
echo "Step 3: Running YaHS scaffolding..."

OUT_PREFIX="yahs_out"

yahs -e ${ENZYME_SEQ} \
     -q 10 \
     -o ${OUT_PREFIX} \
     ${REF_FA} \
     ${FINAL_BAM} 2>&1 | tee yahs.log

if [ ! -f "${OUT_PREFIX}_scaffolds_final.fa" ]; then
    echo "❌ 错误: YaHS 未能生成 scaffolds 文件!"
    exit 1
fi

echo "✅ YaHS finished. Check ${OUT_PREFIX}_scaffolds_final.agp"

# =============================================================================
# --- (D) 📊 结果可视化与后处理 ---
# =============================================================================

# --- 步骤 4: 生成 Hi-C 热图 (用于 Juicebox 浏览) ---
echo "Step 4: Generating .hic file for visualization..."

if [ -f "${OUT_PREFIX}.bin" ]; then
    echo "转换 YaHS 输出为 Juicer 格式..."
    
    juicer pre ${OUT_PREFIX}.bin \
               ${OUT_PREFIX}_scaffolds_final.agp \
               ${REF_FA}.fai 2> juicer_pre.log | \
    sort -k2,2d -k6,6d -T ./ --parallel=${THREADS} -S32G | \
    awk 'NF' > alignments_sorted.txt
    
    if [ ! -s alignments_sorted.txt ]; then
        echo "❌ 错误: alignments_sorted.txt 为空!"
        exit 1
    fi
    
    # 生成 chrom.sizes
    cut -f1,2 ${REF_FA}.fai > chrom.sizes
    
    echo "生成 .hic 文件..."
    java -Xmx120G -Xms32G -jar ${JUICER_JAR} pre \
        alignments_sorted.txt \
        ${OUT_PREFIX}_final.hic \
        chrom.sizes 2>&1 | tee juicer_tools.log
        
    if [ -f "${OUT_PREFIX}_final.hic" ]; then
        echo "✅ Generated: ${OUT_PREFIX}_final.hic"
        rm alignments_sorted.txt
    else
        echo "❌ 错误: .hic 文件生成失败!"
        exit 1
    fi
else
    echo "❌ Error: YaHS .bin file not found!"
    exit 1
fi

# --- 步骤 5: 生成 JBAT 文件 (用于手动纠错) ---
echo "Step 5: Generating JBAT files for manual curation..."

juicer pre -a -o out_JBAT \
           ${OUT_PREFIX}.bin \
           ${OUT_PREFIX}_scaffolds_final.agp \
           ${REF_FA}.fai > out_JBAT.log 2>&1

if [ -f "out_JBAT.txt" ] && [ -s "out_JBAT.txt" ]; then
    # 提取 assembly size
    if grep -q "PRE_C_SIZE" out_JBAT.log; then
        grep "PRE_C_SIZE" out_JBAT.log | awk '{print $2" "$3}' > jbat_chrom_sizes.txt
    else
        # 备用方案: 使用 scaffold 的 fai
        cut -f1,2 ${REF_FA}.fai > jbat_chrom_sizes.txt
    fi
    
    java -Xmx120G -Xms32G -jar ${JUICER_JAR} pre \
        out_JBAT.txt \
        out_JBAT.hic \
        jbat_chrom_sizes.txt 2>&1 | tee juicer_jbat.log
        
    if [ -f "out_JBAT.hic" ]; then
        echo "✅ Generated JBAT files for manual curation"
    else
        echo "⚠️  Warning: JBAT .hic generation failed"
    fi
else
    echo "⚠️  Warning: Failed to generate JBAT text file."
fi

# =============================================================================
# --- (E) 📈 统计信息 ---
# =============================================================================
echo ""
echo "📊 Assembly Statistics:"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

# 统计 scaffold 数量和长度
if [ -f "${OUT_PREFIX}_scaffolds_final.fa" ]; then
    n_scaffolds=$(grep -c "^>" ${OUT_PREFIX}_scaffolds_final.fa)
    total_length=$(awk '/^>/ {next} {sum += length($0)} END {print sum}' ${OUT_PREFIX}_scaffolds_final.fa)
    
    echo "Scaffold 数量: ${n_scaffolds}"
    echo "总长度: ${total_length} bp"
    
    # 计算 N50
    awk '/^>/ {if (seq) print length(seq); seq=""; next} {seq=seq$0} END {if (seq) print length(seq)}' \
        ${OUT_PREFIX}_scaffolds_final.fa | \
        sort -rn | \
        awk -v total=${total_length} '{
            len[NR]=$1; sum+=$1
        } END {
            for (i=1; i<=NR; i++) {
                cum+=len[i]
                if (cum >= total/2) {
                    print "N50: " len[i] " bp"
                    break
                }
            }
        }'
fi

# =============================================================================
# --- (F) 🏁 结束总结 ---
# =============================================================================
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "🎉 Pipeline Completed Successfully!"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "主要输出文件:"
echo "  📄 Final Fasta: ${OUT_PREFIX}_scaffolds_final.fa"
echo "  📄 AGP File:    ${OUT_PREFIX}_scaffolds_final.agp"
echo "  📄 Hi-C Map:    ${OUT_PREFIX}_final.hic"
echo "  📄 JBAT Files:  out_JBAT.hic & out_JBAT.assembly"
echo ""
echo "日志文件:"
echo "  📋 bwa_mem.log, samtools_*.log, yahs.log"
echo "  📋 alignment_stats.txt"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "ℹ️  INFO: 作业结束于: $(date)"