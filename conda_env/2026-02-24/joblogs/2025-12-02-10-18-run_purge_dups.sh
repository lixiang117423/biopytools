#!/bin/bash

# ================= 1. 路径配置 (请核对) =================

# --- 本地软件路径 ---
MINIMAP2_BIN="/share/org/YZWL/yzwl_lixg/.local/bin/minimap2"

# --- Singularity 配置 ---
SINGULARITY_BIN="/share/org/YZWL/yzwl_lixg/miniforge3/envs/singularity_v.3.8.7/bin/singularity"
SIF_IMAGE="/share/org/YZWL/yzwl_lixg/software/singularity/purge_dups_v1.2.6.sif"

# --- 输入/输出文件 ---
PRI_ASM="/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/11.去冗余后重新挂载/01.data/OV53_1.primary.fa"
OUT_DIR="/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/11.去冗余后重新挂载/02.去冗余"

# --- 【必须修改】测序 Reads 数据路径 ---
# 请填入你的 PacBio HiFi 或 CLR 或 Nanopore 的 fastq/fasta 文件路径
# 确保文件也在 /share 目录下，或者修改下面的挂载参数
READS="/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/01.data/hifi/OV53_1-hifi.fq" 

# --- 参数设置 ---
# 数据类型: map-hifi (PacBio HiFi), map-pb (PacBio CLR), map-ont (Nanopore)
PRESET="map-hifi"
THREADS=64

# ================= 2. 环境检查与函数 =================

# 检查 Minimap2
if [ ! -x "$MINIMAP2_BIN" ]; then
    echo "Error: 本地 minimap2 不存在或无执行权限: $MINIMAP2_BIN"
    exit 1
fi

# 检查 Singularity
if [ ! -f "$SINGULARITY_BIN" ] || [ ! -f "$SIF_IMAGE" ]; then
    echo "Error: Singularity 程序或镜像文件不存在。"
    exit 1
fi

# 定义容器运行函数 (挂载 /share 目录)
run_sif() {
    "$SINGULARITY_BIN" exec -B /share "$SIF_IMAGE" "$@"
}

mkdir -p ${OUT_DIR}
cd ${OUT_DIR}

echo "=================================================="
echo "Start Purge_dups Pipeline (Hybrid Mode)"
echo "Assembly: ${PRI_ASM}"
echo "Minimap2: ${MINIMAP2_BIN}"
echo "Image:    ${SIF_IMAGE}"
echo "=================================================="

# ================= 3. 执行流程 =================

# Step 1: Reads 比对 (使用本地 Minimap2)
echo "[1/5] Mapping reads to assembly (Local Minimap2)..."
$MINIMAP2_BIN -x ${PRESET} -t ${THREADS} "${PRI_ASM}" "${READS}" | gzip -c - > pb_aln.paf.gz

# Step 2: 统计覆盖度 (使用容器内的 pbcstat/calcuts)
echo "[2/5] Calculating coverage statistics (Singularity)..."
run_sif pbcstat pb_aln.paf.gz 
# 这会生成 PB.base.cov 和 PB.stat

echo "Calculating cutoffs..."
run_sif calcuts PB.stat > cutoffs 2>calcults.log

echo "--- Cutoffs Content ---"
cat cutoffs
echo "-----------------------"

# Step 3: 组装自比对 (混合模式)
echo "[3/5] Running self-alignment..."

# 3.1 切分 fasta (使用容器内的 split_fa)
echo "Splitting assembly..."
run_sif split_fa "${PRI_ASM}" > split_assembly.fa

# 3.2 自比对 (使用本地 Minimap2)
echo "Self-mapping with Minimap2..."
$MINIMAP2_BIN -x asm5 -DP -t ${THREADS} split_assembly.fa split_assembly.fa | gzip -c - > self_aln.paf.gz

# Step 4: 执行去冗余 (使用容器内的 purge_dups)
echo "[4/5] Identifying duplicates (Singularity)..."
run_sif purge_dups -2 -T cutoffs -c PB.base.cov self_aln.paf.gz > dups.bed

# Step 5: 提取序列 (使用容器内的 get_seqs)
echo "[5/5] Extracting purged sequences (Singularity)..."
run_sif get_seqs -e dups.bed "${PRI_ASM}" -p OV53_purged

echo "=================================================="
echo "Done!"
echo "最终用于 ALLHiC 的去冗余文件: ${OUT_DIR}/OV53_purged.purged.fa"
echo "被剔除的单倍型序列:           ${OUT_DIR}/OV53_purged.hap.fa"
echo "=================================================="