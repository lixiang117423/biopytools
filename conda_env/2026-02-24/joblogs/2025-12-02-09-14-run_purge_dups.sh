#!/bin/bash

# ================= 配置区域 (请修改这里) =================
source /share/org/YZWL/yzwl_lixg/miniforge3/etc/profile.d/conda.sh
conda activate purge_dups_v.1.2.6

# 1. 输入的组装文件 (你提供的路径)
PRI_ASM="/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/11.去冗余后重新挂载/01.data/OV53_1.primary.fa"

# 2. 输出目录 (你提供的路径)
OUT_DIR="/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/11.去冗余后重新挂载/02.去冗余"

# 3. 【重要】三代测序原始 Reads 路径 (请修改为实际路径)
#    可以是 .fastq.gz, .fasta.gz 或者 .bam
#    如果是 Hifi 数据，建议用 .fastq.gz
READS="/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/01.data/hifi/OV53_1-hifi.fq" 

# 4. 数据类型 (根据你的数据选择)
#    PacBio HiFi 选: map-hifi
#    PacBio CLR  选: map-pb
#    Nanopore    选: map-ont
PRESET="map-hifi"

# 5. 线程数
THREADS=64

# ==========================================================

# 创建输出目录并进入
mkdir -p ${OUT_DIR}
cd ${OUT_DIR}

echo "Start running purge_dups pipeline..."
echo "Input Assembly: ${PRI_ASM}"
echo "Output Directory: ${OUT_DIR}"

# Step 1: 将原始 Reads 比对到组装结果上 (计算覆盖度)
# 如果你已经有了 bam 文件，可以跳过 minimap2 这一步，直接用 bam
echo "[1/5] Mapping reads to assembly for coverage analysis..."
minimap2 -x ${PRESET} -t ${THREADS} ${PRI_ASM} ${READS} | gzip -c - > pb_aln.paf.gz

# Step 2: 统计覆盖度并计算阈值
echo "[2/5] Calculating coverage statistics..."
# pbcstat 会生成 PB.base.cov 和 PB.stat 文件
pbcstat pb_aln.paf.gz 
# calcuts 根据统计结果计算 cutoffs，生成 cutoffs 文件
calcuts PB.stat > cutoffs 2>calcults.log

echo "Coverage cutoffs calculated:"
cat cutoffs

# Step 3: 组装结果自比对 (Self-alignment)
# 这一步是为了发现组装中的冗余片段
echo "[3/5] Running self-alignment..."
# 先切分 fasta (Split assembly by 'N')
split_fa ${PRI_ASM} > split_assembly.fa
# 自比对
minimap2 -x asm5 -DP -t ${THREADS} split_assembly.fa split_assembly.fa | gzip -c - > self_aln.paf.gz

# Step 4: 执行去冗余 (Purge Dups)
echo "[4/5] Identifying duplicates..."
# -2 表示只进行两轮 purge，-T 指定阈值文件，-c 指定覆盖度文件
purge_dups -2 -T cutoffs -c PB.base.cov self_aln.paf.gz > dups.bed

# Step 5: 提取序列 (Get Sequences)
echo "[5/5] Extracting purged sequences..."
# 这会生成 purged.fa (去冗余后的主基因组) 和 hap.fa (剔除下来的单倍型)
get_seqs -e dups.bed ${PRI_ASM} -p OV53_purged

echo "Done!"
echo "------------------------------------------------"
echo "去冗余后的最终文件为: ${OUT_DIR}/OV53_purged.purged.fa"
echo "剔除掉的单倍型文件为: ${OUT_DIR}/OV53_purged.hap.fa"
echo "请使用 OV53_purged.purged.fa 进行后续的 ALLHiC 分析"
