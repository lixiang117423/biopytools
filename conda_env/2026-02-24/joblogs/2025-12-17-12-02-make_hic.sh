#!/bin/bash

# ================= 变量定义 =================
# 输出目录
OUT_DIR="/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/21.ALLHiC_hap/test_pipeline/output/13.asmkit"

# 软件路径
ASMKIT="/share/org/YZWL/yzwl_lixg/software/asmkit/asmkit"
VISUALIZER="/share/org/YZWL/yzwl_lixg/software/3d-dna/visualize/run-assembly-visualizer.sh"

# 输入文件路径
INPUT_BAM="/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/21.ALLHiC_hap/test_pipeline/output/01_mapping/sample.clean.bam"
INPUT_AGP="/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/21.ALLHiC_hap/test_pipeline/output/08_build/groups.agp"

# ================= 执行逻辑 =================

# 1. 准备目录
if [ ! -d "$OUT_DIR" ]; then
    echo "创建输出目录: $OUT_DIR"
    mkdir -p "$OUT_DIR"
fi

echo "进入输出目录..."
# cd "$OUT_DIR" || exit 1

# 2. 生成 links 文件 (BAM -> links)
# 这一步可能会比较耗时，取决于 BAM 文件的大小
echo ">>> Step 1: Generating out.links from BAM..."
$ASMKIT bam2links "$INPUT_BAM" out.links

if [ $? -ne 0 ]; then echo "Error: bam2links failed"; exit 1; fi


# 3. 生成 assembly 文件 (AGP -> assembly)
echo ">>> Step 2: Generating groups.assembly from AGP..."
$ASMKIT agp2assembly "$INPUT_AGP" groups.assembly

if [ $? -ne 0 ]; then echo "Error: agp2assembly failed"; exit 1; fi


# 4. 生成 hic 文件 (assembly + links -> hic)
# 这一步需要 Java 环境，生成的文件通常名为 groups.hic
echo ">>> Step 3: Running Visualizer to generate .hic file..."
bash "$VISUALIZER" groups.assembly out.links

if [ $? -ne 0 ]; then echo "Error: 3d-dna visualizer failed"; exit 1; fi

echo "=========================================="
echo "运行完成！文件已生成在: $OUT_DIR"
echo "请下载以下两个文件到本地导入 Juicebox:"
echo "1. $OUT_DIR/groups.assembly"
echo "2. $OUT_DIR/groups.hic"
echo "=========================================="
