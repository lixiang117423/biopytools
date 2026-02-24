#!/bin/bash
# Hi-C数据质控 - 使用HiCExplorer（修复版）
# 工作目录：/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/62.hic/20251011_test

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

# 限制性内切酶配置
RESTRICTION_ENZYME="MboI"
RESTRICTION_SITE="GATC"  # MboI的识别位点

THREADS=16
BIN_SIZE=10000  # 小测用10kb bin

echo "=========================================="
echo "HiCExplorer Hi-C质控分析（修复版）"
echo "=========================================="
echo "样品: ${SAMPLE}"
echo "限制性酶: ${RESTRICTION_ENZYME} (${RESTRICTION_SITE})"
echo ""

#==========================================
# 步骤1：创建目录
#==========================================

echo "步骤1: 创建工作目录..."
mkdir -p mapping
mkdir -p qc_results
mkdir -p matrix
mkdir -p reports
mkdir -p annotation

#==========================================
# 步骤2：生成限制性酶切位点BED文件
#==========================================

echo ""
echo "步骤2: 生成限制性酶切位点文件..."

RESTRICTION_BED="annotation/${GENOME_NAME}_${RESTRICTION_ENZYME}_sites.bed"

if [ ! -f "${RESTRICTION_BED}" ]; then
    echo "  扫描基因组中的 ${RESTRICTION_SITE} 位点..."
    
    # 使用Python脚本生成BED文件
    python3 << EOF
import re
import sys

genome_file = "${GENOME_FASTA}"
restriction_site = "${RESTRICTION_SITE}"
output_file = "${RESTRICTION_BED}"

print(f"正在扫描 {restriction_site} 位点...")

site_count = 0
with open(output_file, 'w') as out:
    chrom = ""
    sequence = ""
    
    with open(genome_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # 处理上一条染色体
                if chrom and sequence:
                    seq_upper = sequence.upper()
                    
                    # 找到所有识别位点
                    for match in re.finditer(restriction_site, seq_upper):
                        # 切割位点在识别序列之后
                        cut_pos = match.end()
                        # BED格式：chr start end
                        out.write(f"{chrom}\t{cut_pos}\t{cut_pos}\n")
                        site_count += 1
                    
                    print(f"  {chrom}: {site_count} 个位点", file=sys.stderr)
                    site_count = 0
                
                # 开始新染色体
                chrom = line[1:].split()[0]
                sequence = ""
            else:
                sequence += line
        
        # 处理最后一条染色体
        if chrom and sequence:
            seq_upper = sequence.upper()
            for match in re.finditer(restriction_site, seq_upper):
                cut_pos = match.end()
                out.write(f"{chrom}\t{cut_pos}\t{cut_pos}\n")
                site_count += 1
            print(f"  {chrom}: {site_count} 个位点", file=sys.stderr)

print(f"✓ 限制性酶切位点文件已生成: {output_file}", file=sys.stderr)
EOF
    
    echo "  ✓ BED文件生成完成: ${RESTRICTION_BED}"
    
    # 查看前几行验证
    echo "  前5个酶切位点："
    head -5 ${RESTRICTION_BED}
else
    echo "  ✓ BED文件已存在: ${RESTRICTION_BED}"
fi

#==========================================
# 步骤3：构建基因组索引
#==========================================

echo ""
echo "步骤3: 构建Bowtie2索引..."

if [ ! -f "mapping/${GENOME_NAME}.1.bt2" ]; then
    echo "  构建索引中..."
    bowtie2-build ${GENOME_FASTA} mapping/${GENOME_NAME}
    echo "  ✓ 索引构建完成"
else
    echo "  ✓ 索引已存在"
fi

#==========================================
# 步骤4：比对reads
#==========================================

echo ""
echo "步骤4: 比对测序数据..."

# 比对R1
if [ ! -f "mapping/${SAMPLE}_R1.bam" ]; then
    echo "  比对 Read 1..."
    bowtie2 -p ${THREADS} \
            --very-sensitive \
            -x mapping/${GENOME_NAME} \
            -U ${READ1} 2> mapping/${SAMPLE}_R1.log | \
    samtools view -bS - | \
    samtools sort -@ ${THREADS} -o mapping/${SAMPLE}_R1.bam -
    
    samtools index mapping/${SAMPLE}_R1.bam
    echo "  ✓ R1比对完成"
else
    echo "  ✓ R1已比对"
fi

# 比对R2
if [ ! -f "mapping/${SAMPLE}_R2.bam" ]; then
    echo "  比对 Read 2..."
    bowtie2 -p ${THREADS} \
            --very-sensitive \
            -x mapping/${GENOME_NAME} \
            -U ${READ2} 2> mapping/${SAMPLE}_R2.log | \
    samtools view -bS - | \
    samtools sort -@ ${THREADS} -o mapping/${SAMPLE}_R2.bam -
    
    samtools index mapping/${SAMPLE}_R2.bam
    echo "  ✓ R2比对完成"
else
    echo "  ✓ R2已比对"
fi

#==========================================
# 步骤5：构建Hi-C矩阵（修复版）
#==========================================

echo ""
echo "步骤5: 构建Hi-C矩阵和质控..."

# 删除可能损坏的旧文件
if [ -f "matrix/${SAMPLE}_${BIN_SIZE}.h5" ]; then
    echo "  删除旧的矩阵文件..."
    rm -f matrix/${SAMPLE}_${BIN_SIZE}.h5
fi

echo "  构建矩阵中（这需要一些时间）..."

# 使用正确的参数：需要三个限制性酶参数！
hicBuildMatrix \
    --samFiles mapping/${SAMPLE}_R1.bam mapping/${SAMPLE}_R2.bam \
    --binSize ${BIN_SIZE} \
    --restrictionCutFile ${RESTRICTION_BED} \
    --restrictionSequence ${RESTRICTION_SITE} \
    --danglingSequence ${RESTRICTION_SITE} \
    --outFileName matrix/${SAMPLE}_${BIN_SIZE}.h5 \
    --outBam matrix/${SAMPLE}.bam \
    --QCfolder qc_results/ \
    --threads ${THREADS} 2>&1 | tee matrix/build_matrix.log

if [ $? -eq 0 ]; then
    echo "  ✓ 矩阵构建完成"
else
    echo "  ✗ 矩阵构建出错，请查看 matrix/build_matrix.log"
    exit 1
fi

#==========================================
# 步骤6：生成QC报告
#==========================================

echo ""
echo "步骤6: 生成质控报告..."

if [ -f "qc_results/${SAMPLE}_R1.bam_QC.log" ]; then
    echo "  生成HTML质控报告..."
    
    hicQC \
        --logfiles qc_results/*.log \
        --outputFolder reports/ \
        --labels ${SAMPLE} 2>&1 | tee reports/qc.log
    
    echo ""
    echo "  ✓ 质控报告生成完成！"
    echo "  📊 查看报告: reports/hicQC.html"
else
    echo "  ⚠ 未找到QC日志文件"
fi

#==========================================
# 步骤7：矩阵质控和过滤
#==========================================

echo ""
echo "步骤7: 矩阵质控..."

if [ -f "matrix/${SAMPLE}_${BIN_SIZE}.h5" ]; then
    echo "  生成诊断图..."
    
    # 检查h5文件是否有效
    python3 -c "import tables; f = tables.open_file('matrix/${SAMPLE}_${BIN_SIZE}.h5', 'r'); f.close()" 2>/dev/null
    
    if [ $? -eq 0 ]; then
        hicCorrectMatrix diagnostic_plot \
            -m matrix/${SAMPLE}_${BIN_SIZE}.h5 \
            -o reports/diagnostic_plot.png 2>&1 | tee reports/diagnostic.log || echo "  ⚠ 诊断图生成可能有警告"
        
        if [ -f "reports/diagnostic_plot.png" ]; then
            echo "  ✓ 诊断图: reports/diagnostic_plot.png"
        fi
    else
        echo "  ⚠ H5文件可能损坏，跳过诊断图"
    fi
else
    echo "  ⚠ 矩阵文件不存在"
fi

#==========================================
# 步骤8：可视化contact map
#==========================================

echo ""
echo "步骤8: 生成contact map可视化..."

if [ -f "matrix/${SAMPLE}_${BIN_SIZE}.h5" ]; then
    echo "  绘制contact map..."
    
    hicPlotMatrix \
        --matrix matrix/${SAMPLE}_${BIN_SIZE}.h5 \
        --outFileName reports/contact_map.png \
        --log1p \
        --dpi 300 \
        --title "${SAMPLE} Contact Map" \
        --colorMap RdYlBu_r 2>&1 | tee reports/plot.log || echo "  ⚠ 绘图可能有警告"
    
    if [ -f "reports/contact_map.png" ]; then
        echo "  ✓ Contact map: reports/contact_map.png"
    fi
fi

#==========================================
# 步骤9：提取关键质控指标
#==========================================

echo ""
echo "=========================================="
echo "步骤9: 提取关键质控指标"
echo "=========================================="

python3 << 'PYEOF'
import os
import re

sample = "ov53-1-HIC1"
qc_dir = "qc_results"

print("\n" + "="*70)
print(" "*20 + "Hi-C 质控报告")
print("="*70)
print(f"\n样品: {sample}")

# 查找QC日志文件
log_files = []
if os.path.exists(qc_dir):
    for f in os.listdir(qc_dir):
        if f.endswith('_QC.log'):
            log_files.append(os.path.join(qc_dir, f))

if not log_files:
    print("\n⚠ 未找到QC日志文件")
    exit(1)

print(f"\n找到 {len(log_files)} 个QC日志文件")

stats = {}

for log_file in log_files:
    with open(log_file) as f:
        content = f.read()
        
        # 提取所有关键统计
        patterns = {
            'Total reads pairs': r'Total reads pairs:\s+([\d,]+)',
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
                if key in stats:
                    stats[key] += value
                else:
                    stats[key] = value

if stats:
    print("\n" + "-"*70)
    print("关键质量指标")
    print("-"*70)
    
    total_reads = stats.get('Total reads', 0)
    mapped_reads = stats.get('Mapped reads', 0)
    valid_pairs = stats.get('Valid pairs', 0)
    unmapped = stats.get('Unmapped reads', 0)
    
    if total_reads > 0:
        mapping_rate = (mapped_reads / total_reads) * 100
        valid_rate = (valid_pairs / total_reads) * 100
        unmapped_rate = (unmapped / total_reads) * 100
        
        print(f"\n总reads数:         {total_reads:>15,}")
        print(f"未比对reads:       {unmapped:>15,}  ({unmapped_rate:>6.2f}%)")
        print(f"比对reads:         {mapped_reads:>15,}  ({mapping_rate:>6.2f}%)")
        print(f"有效配对:          {valid_pairs:>15,}  ({valid_rate:>6.2f}%)")
        
        # 过滤统计
        print(f"\n过滤统计:")
        for key in ['Same fragment', 'Self circles', 'Dangling ends', 'Self ligation']:
            if key in stats:
                count = stats[key]
                rate = (count / total_reads) * 100
                print(f"  {key:.<30} {count:>12,}  ({rate:>5.2f}%)")
        
        # 质量评估
        print("\n" + "="*70)
        print("质量评估")
        print("="*70)
        
        print(f"\n比对率: {mapping_rate:.2f}%")
        if mapping_rate >= 80:
            print("  ✓ 优秀")
        elif mapping_rate >= 70:
            print("  ✓ 良好")
        else:
            print("  ⚠ 偏低")
        
        print(f"\n有效配对率: {valid_rate:.2f}%")
        if valid_rate >= 60:
            print("  ✓✓✓ 优秀！")
        elif valid_rate >= 40:
            print("  ✓✓ 良好！")
        elif valid_rate >= 30:
            print("  ✓ 合格")
        else:
            print("  ⚠ 偏低")
        
        # 结论
        print("\n" + "="*70)
        print("结论")
        print("="*70)
        
        if valid_rate >= 40:
            print(f"\n✅ 数据质量合格！")
            print(f"\n反馈测序公司：")
            print(f"【小测验证OK，有效配对率{valid_rate:.1f}%，请安排大测】")
        elif valid_rate >= 30:
            print(f"\n⚠️ 数据质量一般，建议咨询测序公司")
        else:
            print(f"\n❌ 数据质量较差，建议检查实验流程")

print("\n" + "="*70 + "\n")
PYEOF

#==========================================
# 完成
#==========================================

echo ""
echo "=========================================="
echo "分析完成！"
echo "=========================================="
echo ""
echo "重要输出文件："
echo "  1. 📊 HTML质控报告: reports/hicQC.html"
echo "  2. 🖼️  Contact Map:   reports/contact_map.png"
echo "  3. 📈 诊断图:        reports/diagnostic_plot.png"
echo "  4. 📋 详细日志:      qc_results/*.log"
echo ""