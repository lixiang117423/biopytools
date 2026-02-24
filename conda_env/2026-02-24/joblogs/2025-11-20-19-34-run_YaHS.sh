#!/bin/bash
# =============================================================================
#  🧬 YaHS 基因组组装全流程 - 官方文档标准版 (v3.0)
#  功能: BWA比对 -> YaHS挂载 -> 自动修复热图 -> 生成JBAT编辑文件
# =============================================================================

# =============================================================================
# ---  💻 环境与参数配置 ---
# =============================================================================
set -e 
set -o pipefail

echo "ℹ️  INFO: 作业开始于: $(date)"
echo "🖥️  INFO: 运行于计算节点: $(hostname)"

# 1. 📂 路径设置 (请根据实际情况修改)
WORK_DIR="/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/04.YaHS/re_do_yahs"
REF_FA="/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/04.YaHS/OV53_1.primary.fa"
R1_FQ="/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/04.YaHS/fastq/OV53_1-hic_R1.fastq.gz"
R2_FQ="/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/04.YaHS/fastq/OV53_1-hic_R2.fastq.gz"

# 2. ⚙️ 软件路径 & 参数
# Java版 Juicer Tools (用于生成 .hic)
JUICER_JAR="/share/org/YZWL/yzwl_lixg/software/juicer/scripts/common/juicer_tools.jar"
# 限制性内切酶 (MboI/DpnII -> GATC, Arima -> GATC,GANTC)
ENZYME_SEQ="GATC" 
THREADS=64

# 环境变量 (确保 yahs, bwa, samtools 在路径中)
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
    # 使用临时目录防止排序爆内存
    mkdir -p tmp_sort
    # 文档建议: 过滤未比对、补充比对、PCR重复
    bwa mem -5SP -t ${THREADS} ${REF_FA} ${R1_FQ} ${R2_FQ} | \
    samtools sort -n -@ ${THREADS} -m 4G -T tmp_sort/nsort - | \
    samtools fixmate -m -@ ${THREADS} - - | \
    samtools sort -@ ${THREADS} -m 4G -T tmp_sort/csort - | \
    samtools markdup -r -@ ${THREADS} - ${FINAL_BAM}
    
    samtools index -@ ${THREADS} ${FINAL_BAM}
    rm -rf tmp_sort
fi

# --- 步骤 3: 运行 YaHS 组装 (Scaffolding) ---
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "Step 3: Running YaHS..."

OUT_PREFIX="yahs_out"

# -e 酶切位点, -q 质量过滤
yahs -e ${ENZYME_SEQ} -q 10 -o ${OUT_PREFIX} ${REF_FA} ${FINAL_BAM} 2>&1 | tee yahs.log

if [ ! -f "${OUT_PREFIX}_scaffolds_final.fa" ]; then echo "❌ YaHS 运行失败"; exit 1; fi

# =============================================================================
# --- (D) 结果处理 (严格遵循文档 Generate HiC contact maps 部分) ---
# =============================================================================

# --- 步骤 4: 生成可视化热图 (.hic) ---
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "Step 4: Generating Visualization Hi-C Map..."

if [ -f "${OUT_PREFIX}.bin" ]; then
    # 4.1 将 BIN 转为 Text (juicer pre 需要原始 contig 的 fai)
    # 注意：这里的 'juicer' 是 YaHS 自带的 C 程序
    juicer pre ${OUT_PREFIX}.bin ${OUT_PREFIX}_scaffolds_final.agp ${REF_FA}.fai 2> juicer_pre_vis.log | \
    sort -k2,2d -k6,6d -T ./ --parallel=${THREADS} -S32G | \
    awk 'NF' > alignments_sorted.txt

    # 4.2 [关键修正] 生成最终 Scaffold 的 chrom.sizes
    # 文档说: "The file for scaffold sizes... can be taken from the first two columns of the FASTA index file."
    # 指的是 FINAL scaffolds 的 index，不是原始 contig 的。
    samtools faidx ${OUT_PREFIX}_scaffolds_final.fa
    cut -f1,2 ${OUT_PREFIX}_scaffolds_final.fa.fai > chrom.sizes.final

    # 4.3 生成 .hic
    java -Xmx120G -Xms32G -jar ${JUICER_JAR} pre \
        alignments_sorted.txt \
        ${OUT_PREFIX}_final.hic \
        chrom.sizes.final
    
    # 清理大文件
    if [ -s "${OUT_PREFIX}_final.hic" ]; then rm alignments_sorted.txt; fi
else
    echo "❌ Error: BIN file missing."
    exit 1
fi

# --- 步骤 5: 生成 JBAT 纠错文件 (Manual curation with Juicebox) ---
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "Step 5: Generating JBAT Files for Curation..."

# 5.1 运行 juicer pre -a (Assembly Mode)
# 文档: "juicer pre -a -o out_JBAT ..."
juicer pre -a -o out_JBAT \
    ${OUT_PREFIX}.bin \
    ${OUT_PREFIX}_scaffolds_final.agp \
    ${REF_FA}.fai > out_JBAT.log 2>&1

# 5.2 获取 JBAT 特有的 Assembly Size
# 文档: "<(cat out_JBAT.log | grep PRE_C_SIZE | awk '{print $2" "$3}')"
if grep -q "PRE_C_SIZE" out_JBAT.log; then
    grep "PRE_C_SIZE" out_JBAT.log | awk '{print $2" "$3}' > jbat_chrom_sizes.txt
else
    # 如果 grep 失败，使用总长度备用
    echo "assembly $(grep -v '>' ${REF_FA} | wc -c)" > jbat_chrom_sizes.txt
fi

# 5.3 生成 JBAT .hic
if [ -f "out_JBAT.txt" ]; then
    java -Xmx120G -Xms32G -jar ${JUICER_JAR} pre \
        out_JBAT.txt \
        out_JBAT.hic \
        jbat_chrom_sizes.txt
fi

# =============================================================================
# --- (E) 自动生成后续处理脚本 ---
# =============================================================================
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "Step 6: Preparing Post-Curation Script..."

cat << EOF > run_post_curation.sh
#!/bin/bash
# -----------------------------------------------------
# YaHS 后续处理脚本 (自动生成)
# 用途: 手动纠错完成后，运行此脚本生成最终基因组
# -----------------------------------------------------
JUICER_CMD="juicer" # 确保指向 YaHS 的 juicer 工具
LIFTOVER="out_JBAT.liftover.agp"
ORIGINAL_FA="${REF_FA}"
REVIEW_ASM="review.assembly"
OUT_NAME="${OUT_PREFIX}_Curated"

if [ ! -f "\${REVIEW_ASM}" ]; then
    echo "错误: 请先上传 review.assembly 文件！"
    exit 1
fi

echo "正在根据 review.assembly 重构基因组..."
\${JUICER_CMD} post -o \${OUT_NAME} \${REVIEW_ASM} \${LIFTOVER} \${ORIGINAL_FA}

echo "完成！最终文件: \${OUT_NAME}.FINAL.fa"
EOF
chmod +x run_post_curation.sh

# =============================================================================
# --- (F) 结束报告 ---
# =============================================================================
echo ""
echo "🎉 流程结束！"
echo "📊 重要输出文件："
echo "  1. 查看用热图: ${OUT_PREFIX}_final.hic"
echo "  2. 纠错用文件: out_JBAT.hic & out_JBAT.assembly"
echo ""
echo "⚠️  注意 (参照文档 NOTE 3):"
SCALE_FACTOR=$(grep "scale factor" out_JBAT.log | head -n 1)
if [ ! -z "$SCALE_FACTOR" ]; then
    echo "  发现大基因组缩放因子: ${SCALE_FACTOR}"
    echo "  👉 在 Juicebox 中，请点击 'Assembly > Set Scale' 并输入上述数值的数字部分。"
else
    echo "  未检测到特殊的缩放因子，无需设置 Set Scale。"
fi
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"