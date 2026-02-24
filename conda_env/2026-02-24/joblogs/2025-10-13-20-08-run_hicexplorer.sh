#!/bin/bash
# Hi-C数据质控完整流程 - 最终修正版 v4 (解决samtools header问题)
# 工作目录：/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/62.hic/20251011_test
# 线程数：88

set -e


#==========================================
# 配置参数
#==========================================

WORK_DIR="/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/62.hic/20251011_test"
cd ${WORK_DIR}

SAMPLE="ov53-1-HIC1"
READ1="E250928002_L01_ov53-1-HIC1_1.fq.gz"
READ2="E250928002_L01_ov53-1-HIC1_2.fq.gz"
GENOME_FASTA="OV53_1.primary.fasta"
GENOME_NAME="OV53_1"

RESTRICTION_SITE="GATC"
DANGLING_SEQ="GATC"
THREADS=88
BIN_SIZE=10000

echo "=========================================="
echo "Hi-C数据质控分析 - 最终修正版 v4 (解决samtools header问题)"
echo "=========================================="
echo "样品: ${SAMPLE}"
echo "限制性酶: MboI (${RESTRICTION_SITE})"
echo "线程数: ${THREADS}"
echo "工作目录: ${WORK_DIR}"
echo ""
date
echo ""

#==========================================
# 步骤1：创建目录
#==========================================

echo "步骤1: 创建目录结构..."
mkdir -p mapping qc_results matrix reports
echo "  ✓ 目录创建完成"

#==========================================
# 步骤2：生成BED文件
#==========================================

echo ""
echo "步骤2: 生成限制性酶切位点BED文件..."
BED_FILE="${GENOME_NAME}_MboI_sites.bed"

if [ ! -f "${BED_FILE}" ]; then
    echo "  使用 hicFindRestSite 生成BED文件..."
    hicFindRestSite --fasta ${GENOME_FASTA} --searchPattern ${RESTRICTION_SITE} -o ${BED_FILE}
    echo "  ✓ BED文件已生成: ${BED_FILE}"
else
    echo "  ✓ BED文件已存在"
fi

#==========================================
# 步骤3：构建Bowtie2索引
#==========================================

echo ""
echo "步骤3: 构建Bowtie2索引..."

if [ ! -f "mapping/${GENOME_NAME}.1.bt2" ]; then
    echo "  构建索引中..."
    bowtie2-build --threads ${THREADS} ${GENOME_FASTA} mapping/${GENOME_NAME}
    echo "  ✓ 索引构建完成"
else
    echo "  ✓ 索引已存在"
fi

#==========================================
# 步骤4：比对reads并清理read名
# <--- 最终修正：samtools view 必须加 -h 参数以保留header！---
#==========================================

echo ""
echo "步骤4: 比对Hi-C reads并清理read名..."

# 运行前，必须清理掉旧的、有问题的BAM文件
rm -f mapping/${SAMPLE}_R1.bam mapping/${SAMPLE}_R2.bam

# Read 1
echo "  比对 Read 1..."
bowtie2 -p ${THREADS} \
        --very-sensitive \
        --reorder \
        -x mapping/${GENOME_NAME} \
        -U ${READ1} 2> mapping/${SAMPLE}_R1.log | \
samtools view -hS - | \
sed 's/\/[12]$//' | \
samtools sort -n -@ 32 -o mapping/${SAMPLE}_R1.bam -
echo "  ✓ R1比对完成"

# Read 2
echo "  比对 Read 2..."
bowtie2 -p ${THREADS} \
        --very-sensitive \
        --reorder \
        -x mapping/${GENOME_NAME} \
        -U ${READ2} 2> mapping/${SAMPLE}_R2.log | \
samtools view -hS - | \
sed 's/\/[12]$//' | \
samtools sort -n -@ 32 -o mapping/${SAMPLE}_R2.bam -
echo "  ✓ R2比对完成"

#==========================================
# 步骤5：构建Hi-C矩阵（最关键！）
#==========================================
# (这部分及之后的部分不需要修改)
echo ""
echo "=========================================="
echo "步骤5: 构建Hi-C矩阵（最关键的一步）"
echo "=========================================="
echo "  这将需要 1-3 小时..."
echo ""

# 清理旧文件
rm -f matrix/${SAMPLE}_${BIN_SIZE}.h5 2>/dev/null
rm -rf qc_results/* 2>/dev/null

# 运行
hicBuildMatrix \
    --samFiles mapping/${SAMPLE}_R1.bam mapping/${SAMPLE}_R2.bam \
    --outFileName matrix/${SAMPLE}_${BIN_SIZE}.h5 \
    --QCfolder qc_results/ \
    --restrictionCutFile ${BED_FILE} \
    --restrictionSequence ${RESTRICTION_SITE} \
    --danglingSequence ${DANGLING_SEQ} \
    --binSize ${BIN_SIZE} \
    --threads ${THREADS} \
    --inputBufferSize 400000 2>&1 | tee matrix/build_matrix.log

# (后面的代码都是完整的，直接使用即可)
# ... [脚本的其余部分保持不变] ...
# 检查结果
echo ""
if [ -f "matrix/${SAMPLE}_${BIN_SIZE}.h5" ] && [ -s "matrix/${SAMPLE}_${BIN_SIZE}.h5" ]; then
    H5_SIZE=$(du -h matrix/${SAMPLE}_${BIN_SIZE}.h5 | cut -f1)
    echo "  ✓ Hi-C矩阵构建成功"
    echo "  ✓ 矩阵文件大小: ${H5_SIZE}"
else
    echo "  ✗ 矩阵构建失败"
    echo "  查看日志: matrix/build_matrix.log"
    tail -30 matrix/build_matrix.log
    exit 1
fi

if ls qc_results/*_QC.log 1> /dev/null 2>&1; then
    QC_COUNT=$(ls qc_results/*_QC.log | wc -l)
    echo "  ✓ 生成了 ${QC_COUNT} 个QC日志文件"
else
    echo "  ⚠ 警告：未生成QC日志文件"
fi

#==========================================
# 步骤6：生成QC报告
#==========================================

echo ""
echo "步骤6: 生成HTML质控报告..."

if ls qc_results/*_QC.log 1> /dev/null 2>&1; then
    hicQC \
        --logfiles qc_results/*_QC.log \
        --outputFolder reports/ \
        --labels ${SAMPLE} 2>&1 | tee reports/hicQC.log

    if [ -f "reports/hicQC.html" ]; then
        echo "  ✓ HTML报告: reports/hicQC.html"
    fi
fi

#==========================================
# 步骤7：矩阵诊断
#==========================================

echo ""
echo "步骤7: 矩阵质量诊断..."

if [ -f "matrix/${SAMPLE}_${BIN_SIZE}.h5" ]; then
    hicCorrectMatrix diagnostic_plot \
        -m matrix/${SAMPLE}_${BIN_SIZE}.h5 \
        -o reports/diagnostic_plot.png 2>&1 | tee reports/diagnostic.log || \
        echo "  ⚠ 诊断图生成警告（可能数据量较小）"

    if [ -f "reports/diagnostic_plot.png" ]; then
        echo "  ✓ 诊断图: reports/diagnostic_plot.png"
    fi
fi

#==========================================
# 步骤8：绘制Contact Map
#==========================================

echo ""
echo "步骤8: 绘制Contact Map..."

if [ -f "matrix/${SAMPLE}_${BIN_SIZE}.h5" ]; then
    # 获取第一条染色体名
    FIRST_CHR=$(python3 -c "
import tables
try:
    h5 = tables.open_file('matrix/${SAMPLE}_${BIN_SIZE}.h5', 'r')
    chroms = h5.root.intervals.chrs[:]
    print(chroms[0].decode('utf-8'))
    h5.close()
except Exception:
    print('chr1')
" 2>/dev/null || echo "chr1")

    echo "  绘制 ${FIRST_CHR} 的contact map..."

    hicPlotMatrix \
        --matrix matrix/${SAMPLE}_${BIN_SIZE}.h5 \
        --outFileName reports/contact_map_${FIRST_CHR}.png \
        --log1p \
        --dpi 200 \
        --title "${SAMPLE} - ${FIRST_CHR}" \
        --colorMap RdYlBu_r \
        --chromosomeOrder ${FIRST_CHR} 2>&1 || \
        echo "  ⚠ 绘图警告"

    if [ -f "reports/contact_map_${FIRST_CHR}.png" ]; then
        echo "  ✓ Contact map: reports/contact_map_${FIRST_CHR}.png"
    fi
fi

#==========================================
# 步骤9：提取质控指标
#==========================================

echo ""
echo "=========================================="
echo "步骤9: 质控统计报告"
echo "=========================================="

python3 << 'PYEOF'
import os
import re
import sys

sample = "ov53-1-HIC1"
qc_dir = "qc_results"

print("\n" + "="*70)
print(" "*20 + "Hi-C 质控报告")
print("="*70)
print(f"\n样品: {sample}")

log_files = []
if os.path.exists(qc_dir):
    for f in os.listdir(qc_dir):
        if f.endswith('_QC.log'):
            log_files.append(os.path.join(qc_dir, f))

if not log_files:
    print("\n⚠ 未找到QC日志文件")
    sys.exit(0)

print(f"找到 {len(log_files)} 个QC日志文件\n")

stats = {}
for log_file in sorted(log_files):
    with open(log_file) as f:
        content = f.read()

        patterns = {
            'Total reads': r'Total reads:\s+([\d,]+)',
            'Unmapped reads': r'Unmapped reads:\s+([\d,]+)',
            'Mapped reads': r'Mapped reads:\s+([\d,]+)',
            'Valid pairs': r'Valid pairs:\s+([\d,]+)',
            'Same fragment': r'Same fragment filter:\s+([\d,]+)',
            'Self circles': r'Self circles filter:\s+([\d,]+)',
            'Dangling ends': r'Dangling ends filter:\s+([\d,]+)',
            'Self ligation': r'Self ligation filter:\s+([\d,]+)',
        }

        for key, pattern in patterns.items():
            match = re.search(pattern, content)
            if match:
                value = int(match.group(1).replace(',', ''))
                stats[key] = stats.get(key, 0) + value

if stats:
    print("="*70)
    print("关键质量指标")
    print("="*70)

    total = stats.get('Total reads', 0)
    mapped = stats.get('Mapped reads', 0)
    valid = stats.get('Valid pairs', 0)

    if total > 0:
        map_rate = mapped / total * 100
        valid_rate = valid / total * 100

        print(f"\n总reads数:     {total:>15,}")
        print(f"比对reads:     {mapped:>15,}  ({map_rate:>6.2f}%)")
        print(f"有效配对:      {valid:>15,}  ({valid_rate:>6.2f}%)")

        print("\n过滤统计:")
        for key in ['Same fragment', 'Self circles', 'Dangling ends', 'Self ligation']:
            if key in stats:
                count = stats[key]
                rate = count / total * 100
                print(f"  {key:.<25} {count:>12,}  ({rate:>5.2f}%)")

        print("\n" + "="*70)
        print("质量评估")
        print("="*70)

        print(f"\n比对率: {map_rate:.2f}%")
        print("  " + ("✓ 优秀" if map_rate >= 80 else "✓ 良好" if map_rate >= 70 else "⚠ 一般"))

        print(f"\n有效配对率: {valid_rate:.2f}%")
        if valid_rate >= 60:
            print("  ✓✓✓ 优秀！")
            result = "优秀"
        elif valid_rate >= 40:
            print("  ✓✓ 良好！")
            result = "良好"
        elif valid_rate >= 30:
            print("  ✓ 合格")
            result = "合格"
        else:
            print("  ⚠ 偏低")
            result = "偏低"

        print("\n" + "="*70)
        print("结论")
        print("="*70)

        if valid_rate >= 40:
            print(f"\n✅ 数据质量{result}！")
            if total < 100e6:
                print(f"\n📧 反馈测序公司：")
                print(f"【小测验证OK，有效配对率{valid_rate:.1f}%，请安排大测】")
            else:
                print("\n✅ 数据量充足，可进行Hi-C辅助基因组组装")
        elif valid_rate >= 30:
            print("\n⚠️ 数据质量一般，建议咨询测序公司")
        else:
            print("\n❌ 数据质量较差，建议检查实验流程")

print("\n" + "="*70 + "\n")
PYEOF

#==========================================
# 完成
#==========================================

echo ""
echo "=========================================="
echo "🎉 分析完成！"
echo "=========================================="
echo ""
echo "重要输出文件:"
echo "  📊 HTML报告: reports/hicQC.html"
echo "  🖼️  Contact Map: reports/contact_map_*.png"
echo "  📈 诊断图: reports/diagnostic_plot.png"
echo "  💾 矩阵文件: matrix/${SAMPLE}_${BIN_SIZE}.h5"
echo ""
echo "日志文件:"
echo "  - 构建日志: matrix/build_matrix.log"
echo "  - QC日志: qc_results/*_QC.log"
echo "  - 比对日志: mapping/${SAMPLE}_*.log"
echo ""
date
echo ""
echo "=========================================="