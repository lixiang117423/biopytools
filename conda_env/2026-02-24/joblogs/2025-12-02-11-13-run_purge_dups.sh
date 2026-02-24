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

# --- 测序 Reads 数据路径 ---
READS="/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/01.data/hifi/OV53_1-hifi.fq" 

# --- 参数设置 ---
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
echo "Start Purge_dups Pipeline (Manual Threshold Mode)"
echo "Assembly: ${PRI_ASM}"
echo "Minimap2: ${MINIMAP2_BIN}"
echo "Reads:    ${READS}"
echo "=================================================="

# ================= 3. 执行流程 =================

# Step 1: Reads 比对 (智能跳过)
# 如果 pb_aln.paf.gz 已经存在且大于 10MB，则跳过
if [ -f "pb_aln.paf.gz" ] && [ $(stat -c%s "pb_aln.paf.gz") -gt 10000000 ]; then
    echo "[1/5] Mapping reads found (pb_aln.paf.gz). Skipping Minimap2..."
else
    echo "[1/5] Mapping reads to assembly (Local Minimap2)..."
    $MINIMAP2_BIN -x ${PRESET} -t ${THREADS} "${PRI_ASM}" "${READS}" | gzip -c - > pb_aln.paf.gz
fi

# Step 2: 统计覆盖度
echo "[2/5] Calculating coverage statistics (Singularity)..."
run_sif pbcstat pb_aln.paf.gz 
# 生成 PB.base.cov 和 PB.stat (虽然手动指定参数不需要calcuts的结果，但pbcstat必须跑)

# (跳过 calcuts，因为自动计算不准确，我们将手动指定)
echo "Skipping automatic calcuts..."

# Step 3: 组装自比对 (智能跳过)
if [ -f "self_aln.paf.gz" ] && [ $(stat -c%s "self_aln.paf.gz") -gt 1000000 ]; then
    echo "[3/5] Self-alignment found (self_aln.paf.gz). Skipping self-mapping..."
else
    echo "[3/5] Running self-alignment..."
    echo "Splitting assembly..."
    run_sif split_fa "${PRI_ASM}" > split_assembly.fa
    
    echo "Self-mapping with Minimap2..."
    $MINIMAP2_BIN -x asm5 -DP -t ${THREADS} split_assembly.fa split_assembly.fa | gzip -c - > self_aln.paf.gz
fi

# Step 4: 执行去冗余 (【核心修改：使用手动阈值】)
echo "[4/5] Identifying duplicates with MANUAL thresholds..."
echo "Applying settings: Low=5, Mid=70, High=100"
echo "Note: Mid=70 意味着深度小于70的区域都被视为正常基因组，只有极高深度的才会被切除。"

# 注意：这里不再使用 -T cutoffs，而是直接指定 -l -m -u
run_sif purge_dups -2 -l 5 -m 70 -u 100 -c PB.base.cov self_aln.paf.gz > dups_manual.bed

# Step 5: 提取序列
echo "[5/5] Extracting purged sequences (Singularity)..."
# 输出前缀改为 OV53_manual，以免覆盖之前的文件，方便对比
run_sif get_seqs -e dups_manual.bed "${PRI_ASM}" -p OV53_manual

echo "=================================================="
echo "Done!"
echo "--------------------------------------------------"
echo "请检查新生成的文件大小 (预期 purged.fa 应在 1.4G 左右):"
ls -lh ${OUT_DIR}/OV53_manual.purged.fa
ls -lh ${OUT_DIR}/OV53_manual.hap.fa
echo "--------------------------------------------------"
echo "如果文件大小正常，请使用以下文件进行 ALLHiC:"
echo "${OUT_DIR}/OV53_manual.purged.fa"
echo "=================================================="