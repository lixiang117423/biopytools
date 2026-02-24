#!/bin/bash

# ================= 1. 路径配置 (请修改这里) =================

# --- 本地软件路径 ---
MINIMAP2_BIN="/share/org/YZWL/yzwl_lixg/.local/bin/minimap2"

# --- Singularity 配置 ---
SINGULARITY_BIN="/share/org/YZWL/yzwl_lixg/miniforge3/envs/singularity_v.3.8.7/bin/singularity"
SIF_IMAGE="/share/org/YZWL/yzwl_lixg/software/singularity/purge_dups_v1.2.6.sif"

# --- 输入/输出文件 ---
# 你的原始组装文件
PRI_ASM="/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/11.去冗余后重新挂载/01.data/OV53_1.primary.fa"
# 输出目录
OUT_DIR="/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/11.去冗余后重新挂载/02.去冗余"

# --- 【必须修改】测序 Reads 数据路径 ---
READS="/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/01.data/hifi/OV53_1-hifi.fq" 

# --- 参数设置 ---
PRESET="map-hifi"
THREADS=64

# ================= 2. 环境检查与函数 =================

# 检查本地软件
if [ ! -x "$MINIMAP2_BIN" ]; then
    echo "Error: 本地 minimap2 不存在或无执行权限: $MINIMAP2_BIN"
    exit 1
fi

# 检查 Singularity
if [ ! -f "$SINGULARITY_BIN" ] || [ ! -f "$SIF_IMAGE" ]; then
    echo "Error: Singularity 程序或镜像文件不存在。"
    exit 1
fi

# 定义容器运行函数 (挂载 /share)
run_sif() {
    "$SINGULARITY_BIN" exec -B /share "$SIF_IMAGE" "$@"
}

# 创建并进入输出目录
mkdir -p ${OUT_DIR}
cd ${OUT_DIR}

echo "=================================================="
echo "Start Purge_dups Pipeline (Manual Cutoffs Mode)"
echo "Work Dir: ${OUT_DIR}"
echo "=================================================="

# ================= 3. 执行流程 =================

# Step 1: Reads 比对 (智能跳过：如果文件存在且大于10MB)
if [ -f "pb_aln.paf.gz" ] && [ $(stat -c%s "pb_aln.paf.gz") -gt 10000000 ]; then
    echo "[1/6] Mapping reads found (pb_aln.paf.gz). Skipping Minimap2..."
else
    echo "[1/6] Mapping reads to assembly (Local Minimap2)..."
    $MINIMAP2_BIN -x ${PRESET} -t ${THREADS} "${PRI_ASM}" "${READS}" | gzip -c - > pb_aln.paf.gz
fi

# Step 2: 统计覆盖度
echo "[2/6] Calculating coverage statistics (Singularity)..."
run_sif pbcstat pb_aln.paf.gz 
# 这一步必须跑，因为它生成 PB.base.cov

# Step 3: 【关键】生成手动阈值文件
# 之前的 calcuts 算得太低，导致删多了。我们手动写一个文件。
# 格式: min low mid high max_cov max_cov
# 我们设定: Low=5, Mid=70 (保护单倍体), High=120
echo "[3/6] Creating manual cutoffs file..."
echo "5	5	70	70	120	120" > manual_cutoffs

echo "Manual cutoffs set to:"
cat manual_cutoffs

# Step 4: 组装自比对 (智能跳过)
if [ -f "self_aln.paf.gz" ] && [ $(stat -c%s "self_aln.paf.gz") -gt 1000000 ]; then
    echo "[4/6] Self-alignment found (self_aln.paf.gz). Skipping self-mapping..."
else
    echo "[4/6] Running self-alignment..."
    echo "Splitting assembly..."
    run_sif split_fa "${PRI_ASM}" > split_assembly.fa
    
    echo "Self-mapping with Minimap2..."
    $MINIMAP2_BIN -x asm5 -DP -t ${THREADS} split_assembly.fa split_assembly.fa | gzip -c - > self_aln.paf.gz
fi

# Step 5: 执行去冗余 (引用手动生成的 manual_cutoffs)
echo "[5/6] Identifying duplicates with MANUAL thresholds..."
# 注意：这里使用 -T 指定我们刚才创建的文件
run_sif purge_dups -2 -T manual_cutoffs -c PB.base.cov self_aln.paf.gz > dups_manual.bed

# Step 6: 提取序列
echo "[6/6] Extracting purged sequences (Singularity)..."
run_sif get_seqs -e dups_manual.bed "${PRI_ASM}" -p OV53_manual

echo "=================================================="
echo "Done!"
echo "--------------------------------------------------"
echo "请检查结果文件大小 (预期 OV53_manual.purged.fa 应在 1.4G 左右):"
ls -lh ${OUT_DIR}/OV53_manual.purged.fa
ls -lh ${OUT_DIR}/OV53_manual.hap.fa
echo "--------------------------------------------------"
echo "如果大小正常，请使用此文件进行 ALLHiC:"
echo "${OUT_DIR}/OV53_manual.purged.fa"
echo "=================================================="