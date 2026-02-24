#!/bin/bash

# =============================================================================
# 🧬 ALLHiC Results to Juicebox/JBAT Pipeline
# =============================================================================
# 功能：将 ALLHiC 的结果转换为 Juicebox 可视化格式以进行手动调整
# 作者：AI Assistant
# 日期：$(date +%Y-%m-%d)
# =============================================================================

# Stop on error
set -e

# ============================ 1. 参数与路径配置 ============================
# 📂 输入文件路径
AGP_FILE="/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/05.allhic_2/groups.agp"
REF_FASTA="/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/05.allhic_2/groups.asm.fasta"
R1_FASTQ="/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/04.YaHS/fastq/OV53_1-hic_R1.fastq.gz"
R2_FASTQ="/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/04.YaHS/fastq/OV53_1-hic_R2.fastq.gz"

# 🛠️ 工具脚本路径
SCRIPT_VISUALIZER="/share/org/YZWL/yzwl_lixg/software/3d-dna/visualize/run-asm-visualizer.sh"
SCRIPT_AGP2ASM="/share/org/YZWL/yzwl_lixg/software/3d-dna/utils/agp2assembly.py"
JUICER_TOOLS="/share/org/YZWL/yzwl_lixg/software/juicer/CPU/common/juicer_tools.jar"

# ⚙️ 系统资源配置
THREADS=88
MEMORY="800G" #用于Sort的内存
JAVA_MEM="-Xmx800g" #用于Java的堆内存

# 📂 输出目录
OUT_DIR="06.juicebox_curation"
mkdir -p ${OUT_DIR}
cd ${OUT_DIR}

echo "🚀 [Start] Pipeline started at $(date)"
echo "📂 Working Directory: $(pwd)"

# ============================ 2. 建库与比对 (Mapping) ============================
echo "🧬 [Step 1] Indexing Reference FASTA..."
if [ ! -f "${REF_FASTA}.bwt" ]; then
    ln -s ${REF_FASTA} ref.fasta
    bwa index ref.fasta
else
    echo "   Index exists, linking..."
    ln -sf ${REF_FASTA} ref.fasta
    ln -sf ${REF_FASTA}.* .
fi

echo "⚔️ [Step 2] Aligning Hi-C reads to Reference (Scaffolds)..."
# 使用 BWA MEM -SP5 (Hi-C 推荐参数)
# 管道直接处理：BWA -> Samtools View (Filter) -> Sort by Name (for mnd conversion)
# 注意：这里我们生成临时的 name-sorted BAM 用于后续处理
bwa mem -SP5 -t ${THREADS} ref.fasta ${R1_FASTQ} ${R2_FASTQ} | \
    samtools view -@ 10 -Shb -F 2316 - | \
    samtools sort -@ 20 -n -o mapped.namesorted.bam -

echo "✅ Mapping complete."

# ============================ 3. 生成 MND 文件 (Format Conversion) ============================
echo "🔄 [Step 3] Converting BAM to Merged_Nodups.txt format..."

# 3d-dna 需要的格式 (short format): <str1> <pos1> <frag1> <str2> <pos2> <frag2> <mapq1> <mapq2>
# str: 0 for forward, 16 for reverse (converted to 0/1)
# 使用 Samtools view + awk 进行流式转换
# 注意：run-asm-visualizer 需要文件按 chr1, pos1, chr2, pos2 排序

samtools view mapped.namesorted.bam | \
awk '
function get_strand(flag) { return and(flag, 16) ? 1 : 0 }
NR%2==1 { 
    r1_ref=$3; r1_pos=$4; r1_str=get_strand($2); r1_mapq=$5 
} 
NR%2==0 { 
    r2_ref=$3; r2_pos=$4; r2_str=get_strand($2); r2_mapq=$5;
    # 仅输出同一染色体或不同染色体都比对上的 reads (虽然 filter 已经做了一部分)
    if (r1_ref != "*" && r2_ref != "*") {
        # 为了排序规范，通常确保 ref1 <= ref2
        if (r1_ref > r2_ref || (r1_ref == r2_ref && r1_pos > r2_pos)) {
            print r2_str, r2_ref, r2_pos, 0, r1_str, r1_ref, r1_pos, 1, r2_mapq, r1_mapq
        } else {
            print r1_str, r1_ref, r1_pos, 0, r2_str, r2_ref, r2_pos, 1, r1_mapq, r2_mapq
        }
    }
}' | \
sort -k2,2 -k3,3n -k6,6 -k7,7n --parallel=${THREADS} -S ${MEMORY} > mapped.mnd.txt

echo "✅ MND file generated: mapped.mnd.txt"

# ============================ 4. 运行 3d-dna Visualizer ============================
echo "📊 [Step 4] Running 3d-dna Assembly Visualizer..."

# 这里的 ref.fasta 必须对应 mapped.mnd.txt 里的坐标
# run-asm-visualizer.sh [options] <assembly-fasta> <mnd-text-file>
# -p true: print mnd file (not needed here as input is mnd)
# -q 1: mapq threshold (default 0 or 1 is fine for visualizer)

bash ${SCRIPT_VISUALIZER} -q 10 -j ${JUICER_TOOLS} ref.fasta mapped.mnd.txt

echo "✅ Visualization files generated (.hic and .assembly)."

# ============================ 5. 转换 AGP 为 Assembly (可选/补充) ============================
echo "📑 [Step 5] Converting ALLHiC AGP to Juicebox Assembly format..."

# ⚠️ 注意：这个生成的 assembly 文件是基于 ALLHiC 输入的 Contigs 的。
# 而上面的 .hic 文件是基于 ALLHiC 输出的 Scaffolds 的。
# 如果您想在 Juicebox 中拆分 Contig，您通常需要：
# 1. 用 Draft Contigs 生成 .hic
# 2. 载入这个 agp2assembly 生成的 .assembly
# 但鉴于您提供了 asm.fasta，这里生成仅作为备份参考。

if [ -f "${AGP_FILE}" ]; then
    python3 ${SCRIPT_AGP2ASM} ${AGP_FILE} allhic_groups.assembly
    echo "✅ Converted AGP to allhic_groups.assembly"
else
    echo "⚠️ Warning: AGP file not found at ${AGP_FILE}"
fi

# ============================ 6. 结束 ============================
echo "🎉 [Done] Workflow finished!"
echo "👉 结果文件："
echo "   1. ref.hic (可直接在 Juicebox 加载)"
echo "   2. ref.assembly (对应的组装结构文件)"
echo "   3. allhic_groups.assembly (源自 AGP 的结构文件)"

# 清理临时文件
# rm mapped.namesorted.bam
