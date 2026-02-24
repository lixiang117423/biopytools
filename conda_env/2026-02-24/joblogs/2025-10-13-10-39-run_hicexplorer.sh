#!/bin/bash
# Hi-C数据质控 - 使用HiCExplorer
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

# ⚠️ 重要：确认限制性内切酶
RESTRICTION_ENZYME="MboI"
RESTRICTION_SITE="GATC"  # MboI的识别位点

THREADS=16
BIN_SIZE=10000  # 小测用10kb bin

echo "=========================================="
echo "HiCExplorer Hi-C质控分析"
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

#==========================================
# 步骤2：构建基因组索引
#==========================================

echo ""
echo "步骤2: 构建Bowtie2索引..."

if [ ! -f "mapping/${GENOME_NAME}.1.bt2" ]; then
    echo "  构建索引中..."
    bowtie2-build ${GENOME_FASTA} mapping/${GENOME_NAME}
    echo "  ✓ 索引构建完成"
else
    echo "  ✓ 索引已存在"
fi

#==========================================
# 步骤3：比对reads
#==========================================

echo ""
echo "步骤3: 比对测序数据..."

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
# 步骤4：构建Hi-C矩阵（最关键！）
#==========================================

echo ""
echo "步骤4: 构建Hi-C矩阵和质控..."

if [ ! -f "matrix/${SAMPLE}_${BIN_SIZE}.h5" ]; then
    echo "  构建矩阵中（这需要一些时间）..."
    
    hicBuildMatrix \
        --samFiles mapping/${SAMPLE}_R1.bam mapping/${SAMPLE}_R2.bam \
        --binSize ${BIN_SIZE} \
        --restrictionSequence ${RESTRICTION_SITE} \
        --danglingSequence ${RESTRICTION_SITE} \
        --outFileName matrix/${SAMPLE}_${BIN_SIZE}.h5 \
        --outBam matrix/${SAMPLE}.bam \
        --QCfolder qc_results/ \
        --threads ${THREADS} 2>&1 | tee matrix/build_matrix.log
    
    echo "  ✓ 矩阵构建完成"
else
    echo "  ✓ 矩阵已存在"
fi

#==========================================
# 步骤5：生成QC报告（最重要！）
#==========================================

echo ""
echo "步骤5: 生成质控报告..."

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
# 步骤6：矩阵质控和过滤
#==========================================

echo ""
echo "步骤6: 矩阵质控和过滤..."

if [ -f "matrix/${SAMPLE}_${BIN_SIZE}.h5" ]; then
    # 检查矩阵质量
    hicCorrectMatrix diagnostic_plot \
        -m matrix/${SAMPLE}_${BIN_SIZE}.h5 \
        -o reports/diagnostic_plot.png
    
    echo "  ✓ 诊断图生成: reports/diagnostic_plot.png"
fi

#==========================================
# 步骤7：可视化contact map
#==========================================

echo ""
echo "步骤7: 生成contact map可视化..."

if [ -f "matrix/${SAMPLE}_${BIN_SIZE}.h5" ]; then
    echo "  绘制contact map..."
    
    hicPlotMatrix \
        --matrix matrix/${SAMPLE}_${BIN_SIZE}.h5 \
        --outFileName reports/contact_map.png \
        --log1p \
        --dpi 300 \
        --title "${SAMPLE} Contact Map" \
        --colorMap RdYlBu_r
    
    echo "  ✓ Contact map: reports/contact_map.png"
fi

#==========================================
# 步骤8：提取关键质控指标
#==========================================

echo ""
echo "=========================================="
echo "步骤8: 提取关键质控指标"
echo "=========================================="

python3 << 'PYEOF'
import os
import re
import json

sample = "ov53-1-HIC1"
qc_dir = "qc_results"
work_dir = "/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/62.hic/20251011_test"

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
    print("请检查hicBuildMatrix是否成功运行")
else:
    print(f"\n找到 {len(log_files)} 个QC日志文件")
    
    stats = {}
    
    for log_file in log_files:
        print(f"\n处理: {os.path.basename(log_file)}")
        with open(log_file) as f:
            content = f.read()
            
            # 提取关键统计
            patterns = {
                'Total_reads': r'Total reads:\s+([\d,]+)',
                'Mapped_reads': r'Mapped reads:\s+([\d,]+)',
                'Valid_pairs': r'Valid pairs:\s+([\d,]+)',
                'Min_distance_pairs': r'Min distance filter:\s+([\d,]+)',
                'Self_circles': r'Self circles filter:\s+([\d,]+)',
                'Dangling_ends': r'Dangling ends filter:\s+([\d,]+)',
                'Same_fragment': r'Same fragment filter:\s+([\d,]+)',
                'Self_ligation': r'Self ligation filter:\s+([\d,]+)',
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
        
        total_reads = stats.get('Total_reads', 0)
        mapped_reads = stats.get('Mapped_reads', 0)
        valid_pairs = stats.get('Valid_pairs', 0)
        
        if total_reads > 0:
            mapping_rate = (mapped_reads / total_reads) * 100
            valid_rate = (valid_pairs / total_reads) * 100
            
            print(f"\n总reads数:         {total_reads:>15,}")
            print(f"比对reads数:       {mapped_reads:>15,}  ({mapping_rate:>6.2f}%)")
            print(f"有效配对数:        {valid_pairs:>15,}  ({valid_rate:>6.2f}%)")
            
            # 过滤统计
            print("\n过滤统计:")
            for key in ['Min_distance_pairs', 'Same_fragment', 'Self_circles', 
                       'Dangling_ends', 'Self_ligation']:
                if key in stats:
                    print(f"  {key:.<30} {stats[key]:>15,}")
            
            # 质量评估
            print("\n" + "-"*70)
            print("质量评估")
            print("-"*70)
            
            print(f"\n📊 有效配对率: {valid_rate:.2f}%")
            if valid_rate >= 60:
                print("   ✓✓✓ 优秀！数据质量非常好")
                result = "优秀"
            elif valid_rate >= 40:
                print("   ✓✓ 良好！达到质控标准")
                result = "良好"
            elif valid_rate >= 30:
                print("   ✓ 合格，可以使用")
                result = "合格"
            else:
                print("   ✗ 质量偏低，需要检查")
                result = "偏低"
            
            print(f"\n📊 比对率: {mapping_rate:.2f}%")
            if mapping_rate >= 80:
                print("   ✓ 优秀")
            elif mapping_rate >= 70:
                print("   ✓ 良好")
            else:
                print("   ⚠ 偏低，检查参考基因组")
            
            # 建议
            print("\n" + "="*70)
            print("结论和建议")
            print("="*70)
            
            if valid_rate >= 40:
                print("\n✓ 数据质量合格！")
                print("\n📧 反馈测序公司：")
                print("   【小测验证OK，有效配对率 {:.1f}%，请安排大测】".format(valid_rate))
            elif valid_rate >= 30:
                print("\n⚠ 数据质量一般")
                print("\n建议：")
                print("   - 咨询测序公司是否可接受")
                print("   - 检查限制性酶设置是否正确")
            else:
                print("\n✗ 数据质量较差")
                print("\n建议：")
                print("   - 联系测序公司讨论")
                print("   - 确认实验protocol")
                print("   - 可能需要重新建库")
    else:
        print("\n⚠ 无法提取统计数据")

print("\n" + "="*70)
print("详细结果文件")
print("="*70)
print(f"\n📁 QC报告:     {work_dir}/reports/hicQC.html")
print(f"📊 Contact Map: {work_dir}/reports/contact_map.png")
print(f"📈 诊断图:      {work_dir}/reports/diagnostic_plot.png")
print(f"💾 矩阵文件:    {work_dir}/matrix/{sample}_10000.h5")
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
echo "  1. 📊 HTML质控报告: reports/hicQC.html (用浏览器打开)"
echo "  2. 🖼️  Contact Map:   reports/contact_map.png"
echo "  3. 📈 诊断图:        reports/diagnostic_plot.png"
echo "  4. 📋 详细日志:      qc_results/*.log"
echo ""
echo "下一步："
echo "  1. 打开 reports/hicQC.html 查看详细质控报告"
echo "  2. 根据有效配对率决定是否安排大测"
echo "  3. 如有问题，查看日志文件排查"
echo ""
