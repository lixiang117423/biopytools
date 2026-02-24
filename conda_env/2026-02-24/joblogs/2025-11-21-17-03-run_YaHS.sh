#!/bin/bash
# =============================================================================
#  🧬 YaHS 基因组组装全流程 - 终极修正版 (v4.0)
#  包含: Mapping -> Scaffolding -> Visualization (.hic) -> Curation (JBAT)
# =============================================================================

# =============================================================================
# ---  💻 环境与参数配置 ---
# =============================================================================
set -e 
set -o pipefail

echo "ℹ️  INFO: 作业开始于: $(date)"
echo "🖥️  INFO: 运行于计算节点: $(hostname)"

# 1. 📂 路径设置 (已填入您提供的路径)
WORK_DIR="/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/04.YaHS/re_do_yahs_2"
REF_FA="/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/04.YaHS/OV53_1.primary.fa"
R1_FQ="/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/04.YaHS/fastq/OV53_1-hic_R1.fastq.gz"
R2_FQ="/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/04.YaHS/fastq/OV53_1-hic_R2.fastq.gz"

# 2. ⚙️ 软件路径 & 参数
JUICER_JAR="/share/org/YZWL/yzwl_lixg/software/juicer/scripts/common/juicer_tools.jar"
ENZYME_SEQ="GATC" # MboI
THREADS=64

# 环境变量
export PATH="/share/org/YZWL/yzwl_lixg/miniforge3/envs/yahs_v.1.2.2/bin:$PATH"

# =============================================================================
# --- (C) 🚀 流程开始 ---
# =============================================================================
mkdir -p ${WORK_DIR}
cd ${WORK_DIR}

# 0. 检查输入
for file in "${REF_FA}" "${R1_FQ}" "${R2_FQ}"; do
    if [ ! -f "${file}" ]; then echo "❌ 错误: 文件不存在 - ${file}"; exit 1; fi
done

# --- 步骤 1: 建立索引 (Indexing) ---
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "Step 1: Checking/Building Indexes..."
if [ ! -f "${REF_FA}.bwt" ]; then bwa index ${REF_FA}; fi
if [ ! -f "${REF_FA}.fai" ]; then samtools faidx ${REF_FA}; fi

# --- 步骤 2: Hi-C 比对与处理 (Mapping & Processing) ---
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "Step 2: Mapping Hi-C reads..."

FINAL_BAM="aligned_sorted_dedup.bam"

if [ -f "${FINAL_BAM}" ]; then
    echo "✅ BAM 文件已存在，跳过比对。"
else
    mkdir -p tmp_sort
    echo "   > 正在运行 BWA mem, sorting, markdup..."
    bwa mem -5SP -t ${THREADS} ${REF_FA} ${R1_FQ} ${R2_FQ} | \
    samtools sort -n -@ ${THREADS} -m 4G -T tmp_sort/nsort - | \
    samtools fixmate -m -@ ${THREADS} - - | \
    samtools sort -@ ${THREADS} -m 4G -T tmp_sort/csort - | \
    samtools markdup -r -@ ${THREADS} - ${FINAL_BAM}
    
    samtools index -@ ${THREADS} ${FINAL_BAM}
    rm -rf tmp_sort
    echo "✅ 比对完成。"
fi

# --- 步骤 3: 运行 YaHS 组装 (Scaffolding) ---
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "Step 3: Running YaHS..."

OUT_PREFIX="yahs_out"

yahs -e ${ENZYME_SEQ} -q 10 -o ${OUT_PREFIX} ${REF_FA} ${FINAL_BAM} 2>&1 | tee yahs.log

if [ ! -f "${OUT_PREFIX}_scaffolds_final.fa" ]; then echo "❌ YaHS 运行失败"; exit 1; fi

# =============================================================================
# --- (D) 结果处理: 生成可视化文件 (.hic) ---
# =============================================================================

# --- 步骤 4: 生成标准热图 (Standard .hic) ---
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "Step 4: Generating Standard Visualization Map..."

if [ -f "${OUT_PREFIX}.bin" ]; then
    # 4.1 转换 bin 为 text
    juicer pre ${OUT_PREFIX}.bin ${OUT_PREFIX}_scaffolds_final.agp ${REF_FA}.fai 2> juicer_pre_vis.log | \
    sort -k2,2d -k6,6d -T ./ --parallel=${THREADS} -S32G | \
    awk 'NF' > alignments_sorted.txt

    # 4.2 [关键] 基于"最终Scaffold"生成 chrom.sizes
    echo "   > Generating correct chrom.sizes from final scaffolds..."
    samtools faidx ${OUT_PREFIX}_scaffolds_final.fa
    cut -f1,2 ${OUT_PREFIX}_scaffolds_final.fa.fai > chrom.sizes.final

    # 4.3 生成 .hic
    java -Xmx120G -Xms32G -jar ${JUICER_JAR} pre \
        alignments_sorted.txt \
        ${OUT_PREFIX}_final.hic \
        chrom.sizes.final
    
    if [ -s "${OUT_PREFIX}_final.hic" ]; then
        rm alignments_sorted.txt
        echo "✅ 标准热图生成成功: ${OUT_PREFIX}_final.hic"
    fi
else
    echo "❌ Error: BIN file missing."
    exit 1
fi

# --- 步骤 5: 生成 JBAT 纠错文件 (Manual curation) ---
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "Step 5: Generating JBAT Files for Curation..."

# 5.1 运行 juicer pre -a (生成 out_JBAT.txt)
echo "   > Generating JBAT text file (can be huge)..."
juicer pre -a -o out_JBAT \
    ${OUT_PREFIX}.bin \
    ${OUT_PREFIX}_scaffolds_final.agp \
    ${REF_FA}.fai > out_JBAT.log 2>&1

# 5.2 [核心修正] 计算 Assembly 总长度
# 强制生成格式为 "assembly <total_length>" 的文件
echo "   > Calculating total assembly size..."
if [ -f "${REF_FA}.fai" ]; then
    awk '{sum+=$2} END {print "assembly", sum}' ${REF_FA}.fai > jbat_chrom_sizes.txt
else
    echo "❌ 错误: 找不到 REF_FA.fai，无法计算长度。"
    exit 1
fi

# 5.3 生成 JBAT .hic
if [ -s "out_JBAT.txt" ]; then
    echo "   > Running Juicer Tools (JBAT Mode)..."
    java -Xmx120G -Xms32G -jar ${JUICER_JAR} pre \
        out_JBAT.txt \
        out_JBAT.hic \
        jbat_chrom_sizes.txt
        
    if [ -s "out_JBAT.hic" ] && [ $(stat -c%s "out_JBAT.hic") -gt 1000000 ]; then
        echo "✅ JBAT 纠错文件生成成功！"
    else
        echo "❌ JBAT .hic 生成失败，文件过小。"
        exit 1
    fi
fi

# =============================================================================
# --- (E) 自动生成后续处理脚本 ---
# =============================================================================
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "Step 6: Preparing Post-Curation Script..."

cat << EOF > run_post_curation.sh
#!/bin/bash
# -----------------------------------------------------
# YaHS 后续处理脚本
# 用途: 在 Juicebox 完成手动纠错并上传 review.assembly 后运行此脚本
# -----------------------------------------------------
JUICER_CMD="juicer" # 确保这是 YaHS 安装包里的 juicer
LIFTOVER="out_JBAT.liftover.agp"
ORIGINAL_FA="${REF_FA}"
REVIEW_ASM="review.assembly"
OUT_NAME="${OUT_PREFIX}_Curated"

# 确保环境变量包含 yahs
export PATH="${PATH}"

if [ ! -f "\${REVIEW_ASM}" ]; then
    echo "❌ 错误: 请先将 Juicebox 保存的 review.assembly 上传到当前目录！"
    exit 1
fi

echo "🚀 正在根据 review.assembly 重构基因组..."
\${JUICER_CMD} post -o \${OUT_NAME} \${REVIEW_ASM} \${LIFTOVER} \${ORIGINAL_FA}

echo "✅ 完成！最终基因组: \${OUT_NAME}.FINAL.fa"
EOF
chmod +x run_post_curation.sh

# =============================================================================
# --- (F) 结束报告 ---
# =============================================================================
echo ""
echo "🎉 流程完美结束！"
echo "📂 请下载以下文件进行人工纠错:"
echo "  1. out_JBAT.hic"
echo "  2. out_JBAT.assembly"
echo ""
echo "⚠️  Juicebox 缩放提示 (Scale Factor):"
SCALE_FACTOR=$(grep "scale factor" out_JBAT.log | head -n 1)
if [ ! -z "$SCALE_FACTOR" ]; then
    echo "  >>> ${SCALE_FACTOR} <<<"
    echo "  (如果 Juicebox 打开是白的，请在 Assembly > Set Scale 中输入上面的数字)"
else
    echo "  无需设置 Scale。"
fi
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"