#!/bin/bash

# ==============================================================================
# HapHiC 染色体挂载全流程自动化脚本 (2025版)
# ==============================================================================

# 设置出错即停止
set -e

# --- 1. 用户自定义配置区域 ---
FASTA="OV53_1.primary.fa"              # 组装的 fasta 文件
R1="OV53_1-hic_R1.fastq.gz"              # Hi-C Read 1
R2="OV53_1-hic_R2.fastq.gz"              # Hi-C Read 2
NCHRS=12                          # 预期的染色体/单倍型组数 (必填)
RE="GATC"                         # 酶切位点 (MboI/DpnII: GATC; Arima: GATC,GANTC)
THREADS=64                        # 线程/进程数
OUTDIR="haphic_out"               # 输出结果目录
CORRECT_ROUNDS=2                  # 纠错轮次 (若组装质量极高可设为0)

# --- 2. 软件绝对路径 (已根据您的环境固定) ---
HAPHIC="/share/org/YZWL/yzwl_lixg/miniforge3/envs/haphic/bin/haphic"
SAMBLASTER="/share/org/YZWL/yzwl_lixg/miniforge3/envs/haphic/bin/samblaster"
FILTER_BAM="/share/org/YZWL/yzwl_lixg/miniforge3/envs/haphic/bin/filter_bam"

# 检查工具是否存在
for tool in "$HAPHIC" "$SAMBLASTER" "$FILTER_BAM"; do
    if [ ! -x "$tool" ]; then
        echo "错误: 工具 $tool 未找到或无执行权限。"
        exit 1
    fi
done

echo "任务开始时间: $(date)"

# --- 3. 建立索引 ---
echo "[Step 1] Indexing reference genome..."
if [ ! -f "${FASTA}.bwt" ]; then
    bwa index "$FASTA"
fi

# --- 4. 比对、去重并初步过滤 ---
# -5SP: Hi-C 推荐参数
# -F 3340: 过滤 secondary, supplementary, unmapped 和 PCR duplicates
echo "[Step 2] Aligning Hi-C reads with BWA and marking duplicates..."
bwa mem -5SP -t "$THREADS" "$FASTA" "$R1" "$R2" | \
    "$SAMBLASTER" | \
    samtools view - -@ "$THREADS" -S -h -b -F 3340 -o raw_HiC.bam

# --- 5. 质量精滤 (MAPQ & 编辑距离) ---
# 使用 HapHiC 提供的 filter_bam，标准: MAPQ >= 1, NM < 3
echo "[Step 3] Filtering BAM (MAPQ 1, NM 3) using HapHiC utils..."
"$FILTER_BAM" raw_HiC.bam 1 --nm 3 --threads "$THREADS" | \
    samtools view - -b -@ "$THREADS" -o HiC.filtered.bam

# 清理中间大文件
rm raw_HiC.bam

# --- 6. 运行 HapHiC Pipeline ---
echo "[Step 4] Running HapHiC assembly pipeline..."
# 包含：纠错(若开启)、聚类、重分配、排序、构建
"$HAPHIC" pipeline "$FASTA" HiC.filtered.bam "$NCHRS" \
    --RE "$RE" \
    --outdir "$OUTDIR" \
    --threads "$THREADS" \
    --processes "$THREADS" \
    --correct_nrounds "$CORRECT_ROUNDS" \
    --verbose

echo "=============================================================================="
echo "HapHiC 流程圆满完成!"
echo "任务结束时间: $(date)"
echo "------------------------------------------------------------------------------"
echo "最终 Scaffold 序列: $OUTDIR/04.build/scaffolds.fa"
echo "YaHS 风格 AGP (最推荐): $OUTDIR/04.build/scaffolds.raw.agp"
echo "手动调整提示: 您可以直接运行 bash $OUTDIR/04.build/juicebox.sh 生成绘图文件"
echo "=============================================================================="
