#!/bin/bash
# Hi-C数据质控完整流程 - 使用HiCExplorer
# 方法2：不使用BED文件，自动识别酶切位点
# 工作目录：/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/62.hic/20251011_test

set -e

#==========================================
# 配置参数
#==========================================

WORK_DIR="/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/62.hic/20251011_test"
cd ${WORK_DIR}

# 输入文件
SAMPLE="ov53-1-HIC1"
READ1="E250928002_L01_ov53-1-HIC1_1.fq.gz"
READ2="E250928002_L01_ov53-1-HIC1_2.fq.gz"
GENOME_FASTA="OV53_1.primary.fasta"
GENOME_NAME="OV53_1"

# 限制性内切酶配置（MboI）
RESTRICTION_SITE="GATC"
DANGLING_SEQ="GATC"

# 计算参数
THREADS=88
BIN_SIZE=10000

echo "=========================================="
echo "Hi-C数据质控分析（方法2：无需BED文件）"
echo "=========================================="
echo "样品: ${SAMPLE}"
echo "限制性酶: MboI (${RESTRICTION_SITE})"
echo "线程数: ${THREADS}"
echo "Bin size: ${BIN_SIZE}"
echo ""

#==========================================
# 步骤1：创建目录结构
#==========================================

echo "步骤1: 创建工作目录..."
mkdir -p mapping
mkdir -p qc_results
mkdir -p matrix
mkdir -p reports

echo "  ✓ 目录创建完成"

#==========================================
# 步骤2：构建Bowtie2索引
#==========================================

echo ""
echo "步骤2: 构建Bowtie2索引..."

if [ ! -f "mapping/${GENOME_NAME}.1.bt2" ]; then
    echo "  构建索引中（这可能需要一些时间）..."
    bowtie2-build --threads ${THREADS} ${GENOME_FASTA} mapping/${GENOME_NAME}
    echo "  ✓ Bowtie2索引构建完成"
else
    echo "  ✓ Bowtie2索引已存在，跳过构建"
fi

#==========================================
# 步骤3：比对Hi-C reads
#==========================================

echo ""
echo "步骤3: 比对Hi-C reads..."

# 比对Read 1
if [ ! -f "mapping/${SAMPLE}_R1.bam" ]; then
    echo "  比对 Read 1..."
    bowtie2 -p ${THREADS} \
            --very-sensitive \
            --reorder \
            -x mapping/${GENOME_NAME} \
            -U ${READ1} 2> mapping/${SAMPLE}_R1.bowtie2.log | \
    samtools view -@ ${THREADS} -bS - | \
    samtools sort -@ ${THREADS} -o mapping/${SAMPLE}_R1.bam -
    
    samtools index -@ ${THREADS} mapping/${SAMPLE}_R1.bam
    echo "  ✓ Read 1 比对完成"
    
    # 显示比对率
    echo "  Read 1 比对统计:"
    grep "overall alignment rate" mapping/${SAMPLE}_R1.bowtie2.log || true
else
    echo "  ✓ Read 1 已比对，跳过"
fi

# 比对Read 2
if [ ! -f "mapping/${SAMPLE}_R2.bam" ]; then
    echo "  比对 Read 2..."
    bowtie2 -p ${THREADS} \
            --very-sensitive \
            --reorder \
            -x mapping/${GENOME_NAME} \
            -U ${READ2} 2> mapping/${SAMPLE}_R2.bowtie2.log | \
    samtools view -@ ${THREADS} -bS - | \
    samtools sort -@ ${THREADS} -o mapping/${SAMPLE}_R2.bam -
    
    samtools index -@ ${THREADS} mapping/${SAMPLE}_R2.bam
    echo "  ✓ Read 2 比对完成"
    
    # 显示比对率
    echo "  Read 2 比对统计:"
    grep "overall alignment rate" mapping/${SAMPLE}_R2.bowtie2.log || true
else
    echo "  ✓ Read 2 已比对，跳过"
fi

#==========================================
# 步骤4：构建Hi-C矩阵（关键步骤！）
#==========================================

echo ""
echo "步骤4: 构建Hi-C矩阵和质控..."
echo "  这是最关键的一步，可能需要1-3小时..."
echo ""

# 清理旧文件
rm -f matrix/${SAMPLE}_${BIN_SIZE}.h5 2>/dev/null
rm -rf qc_results/* 2>/dev/null

# 运行hicBuildMatrix（方法2：不用BED文件）
hicBuildMatrix \
    --samFiles mapping/${SAMPLE}_R1.bam mapping/${SAMPLE}_R2.bam \
    --outFileName matrix/${SAMPLE}_${BIN_SIZE}.h5 \
    --QCfolder qc_results/ \
    --restrictionSequence ${RESTRICTION_SITE} \
    --danglingSequence ${DANGLING_SEQ} \
    --binSize ${BIN_SIZE} \
    --threads ${THREADS} \
    --inputBufferSize 400000 2>&1 | tee matrix/build_matrix.log

# 检查是否成功
if [ -f "matrix/${SAMPLE}_${BIN_SIZE}.h5" ] && [ -s "matrix/${SAMPLE}_${BIN_SIZE}.h5" ]; then
    echo ""
    echo "  ✓ Hi-C矩阵构建成功"
    ls -lh matrix/${SAMPLE}_${BIN_SIZE}.h5
else
    echo ""
    echo "  ✗ 矩阵构建失败，请查看日志: matrix/build_matrix.log"
    exit 1
fi

# 检查QC文件
if ls qc_results/*_QC.log 1> /dev/null 2>&1; then
    echo "  ✓ QC日志文件已生成"
    ls -lh qc_results/
else
    echo "  ⚠ 警告：QC日志文件未生成"
fi

#==========================================
# 步骤5：生成HTML质控报告
#==========================================

echo ""
echo "步骤5: 生成HTML质控报告..."

if ls qc_results/*_QC.log 1> /dev/null 2>&1; then
    hicQC \
        --logfiles qc_results/*_QC.log \
        --outputFolder reports/ \
        --labels ${SAMPLE} 2>&1 | tee reports/hicQC.log
    
    if [ -f "reports/hicQC.html" ]; then
        echo "  ✓ HTML质控报告已生成"
        echo "  📊 报告位置: reports/hicQC.html"
    else
        echo "  ⚠ HTML报告生成可能有问题"
    fi
else
    echo "  ⚠ 跳过：没有QC日志文件"
fi

#==========================================
# 步骤6：矩阵质控诊断
#==========================================

echo ""
echo "步骤6: 矩阵质量诊断..."

if [ -f "matrix/${SAMPLE}_${BIN_SIZE}.h5" ]; then
    echo "  生成诊断图..."
    
    hicCorrectMatrix diagnostic_plot \
        -m matrix/${SAMPLE}_${BIN_SIZE}.h5 \
        -o reports/diagnostic_plot.png 2>&1 | tee reports/diagnostic.log || \
        echo "  ⚠ 诊断图生成有警告（可能是小测数据量较小导致）"
    
    if [ -f "reports/diagnostic_plot.png" ]; then
        echo "  ✓ 诊断图: reports/diagnostic_plot.png"
    fi
fi

#==========================================
# 步骤7：绘制Contact Map
#==========================================

echo ""
echo "步骤7: 绘制Contact Map..."

if [ -f "matrix/${SAMPLE}_${BIN_SIZE}.h5" ]; then
    echo "  绘制contact map..."
    
    hicPlotMatrix \
        --matrix matrix/${SAMPLE}_${BIN_SIZE}.h5 \
        --outFileName reports/contact_map.png \
        --log1p \
        --dpi 300 \
        --title "${SAMPLE} Contact Map (${BIN_SIZE}bp bins)" \
        --colorMap RdYlBu_r \
        --region chr1:1-50000000 2>&1 | tee reports/plot.log || \
        echo "  ⚠ 绘图可能有警告（尝试调整region参数）"
    
    # 如果指定区域失败，尝试全基因组
    if [ ! -f "reports/contact_map.png" ]; then
        echo "  尝试绘制全基因组contact map..."
        hicPlotMatrix \
            --matrix matrix/${SAMPLE}_${BIN_SIZE}.h5 \
            --outFileName reports/contact_map.png \
            --log1p \
            --dpi 200 \
            --title "${SAMPLE} Contact Map" \
            --colorMap RdYlBu_r 2>&1 || echo "  ⚠ 可视化可能需要更多数据"
    fi
    
    if [ -f "reports/contact_map.png" ]; then
        echo "  ✓ Contact map: reports/contact_map.png"
    fi
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
import sys

sample = "ov53-1-HIC1"
qc_dir = "qc_results"

print("\n" + "="*70)
print(" "*20 + "Hi-C 质控报告")
print("="*70)
print(f"\n样品: {sample}")
print(f"分析时间: {os.popen('date').read().strip()}")

# 查找QC日志文件
log_files = []
if os.path.exists(qc_dir):
    for f in os.listdir(qc_dir):
        if f.endswith('_QC.log'):
            log_files.append(os.path.join(qc_dir, f))

if not log_files:
    print("\n⚠ 未找到QC日志文件")
    print("  检查路径: qc_results/")
    sys.exit(1)

print(f"\n找到 {len(log_files)} 个QC日志文件")

# 汇总统计
stats = {}

for log_file in sorted(log_files):
    print(f"\n处理: {os.path.basename(log_file)}")
    
    with open(log_file) as f:
        content = f.read()
        
        # 提取所有关键统计
        patterns = {
            'Total reads pairs': r'Total reads pairs:\s+([\d,]+)',
            'Total reads': r'Total reads:\s+([\d,]+)',
            'Unmapped reads': r'Unmapped reads:\s+([\d,]+)',
            'Low quality': r'Low quality:\s+([\d,]+)',
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
    print("\n" + "="*70)
    print("汇总统计")
    print("="*70)
    
    total_reads = stats.get('Total reads', 0)
    mapped_reads = stats.get('Mapped reads', 0)
    valid_pairs = stats.get('Valid pairs', 0)
    unmapped = stats.get('Unmapped reads', 0)
    
    if total_reads > 0:
        mapping_rate = (mapped_reads / total_reads) * 100
        valid_rate = (valid_pairs / total_reads) * 100
        unmapped_rate = (unmapped / total_reads) * 100
        
        print(f"\n📊 基础统计:")
        print(f"  总reads数:         {total_reads:>15,}")
        print(f"  未比对reads:       {unmapped:>15,}  ({unmapped_rate:>6.2f}%)")
        print(f"  比对reads:         {mapped_reads:>15,}  ({mapping_rate:>6.2f}%)")
        print(f"  有效配对:          {valid_pairs:>15,}  ({valid_rate:>6.2f}%)")
        
        # 过滤统计
        print(f"\n📋 过滤统计:")
        filter_keys = ['Same fragment', 'Self circles', 'Dangling ends', 'Self ligation']
        for key in filter_keys:
            if key in stats:
                count = stats[key]
                rate = (count / total_reads) * 100 if total_reads > 0 else 0
                print(f"  {key:.<30} {count:>12,}  ({rate:>5.2f}%)")
        
        # 质量评估
        print("\n" + "="*70)
        print("质量评估")
        print("="*70)
        
        print(f"\n1️⃣ 比对率: {mapping_rate:.2f}%")
        if mapping_rate >= 80:
            print("   ✓ 优秀")
        elif mapping_rate >= 70:
            print("   ✓ 良好")
        elif mapping_rate >= 60:
            print("   ⚠ 一般，检查参考基因组")
        else:
            print("   ✗ 偏低，参考基因组可能不匹配")
        
        print(f"\n2️⃣ 有效配对率: {valid_rate:.2f}%")
        if valid_rate >= 60:
            print("   ✓✓✓ 优秀！数据质量非常好")
            result = "优秀"
        elif valid_rate >= 40:
            print("   ✓✓ 良好！达到质控标准")
            result = "良好"
        elif valid_rate >= 30:
            print("   ✓ 合格，可以使用")
            result = "合格"
        elif valid_rate >= 20:
            print("   ⚠ 偏低，建议检查")
            result = "偏低"
        else:
            print("   ✗ 很低，可能有问题")
            result = "异常"
        
        # 数据量评估
        print(f"\n3️⃣ 数据量评估:")
        if total_reads < 50e6:
            print(f"   这是小测数据 ({total_reads/1e6:.1f}M reads)")
            print("   用途: 质控验证")
        elif total_reads < 200e6:
            print(f"   这是中等测序量 ({total_reads/1e6:.1f}M reads)")
            print("   用途: 初步分析")
        else:
            print(f"   这是大测数据 ({total_reads/1e6:.1f}M reads)")
            print("   用途: 完整分析")
        
        # 结论和建议
        print("\n" + "="*70)
        print("结论和建议")
        print("="*70)
        
        if valid_rate >= 40 and mapping_rate >= 70:
            print("\n✅ 数据质量合格！")
            if total_reads < 100e6:
                print("\n📧 反馈测序公司：")
                print(f"   【小测验证OK，有效配对率{valid_rate:.1f}%，比对率{mapping_rate:.1f}%，请安排大测】")
                print("\n💡 大测建议:")
                print("   - 建议测序深度: ≥500M reads（用于基因组组装scaffold）")
                print("   - 建议测序深度: ≥1B reads（用于精细TAD分析）")
            else:
                print("\n✅ 数据量充足，可以进行:")
                print("   - Hi-C辅助基因组scaffold（使用Juicer + 3D-DNA）")
                print("   - TAD分析")
                print("   - A/B compartment分析")
        elif valid_rate >= 30:
            print("\n⚠️ 数据质量一般")
            print("\n建议:")
            print("   - 咨询测序公司是否可接受")
            print("   - 确认限制性酶设置正确（当前：MboI/GATC）")
            print("   - 检查参考基因组是否匹配")
        else:
            print("\n❌ 数据质量较差")
            print("\n建议:")
            print("   - 联系测序公司讨论问题")
            print("   - 可能的问题:")
            print("     • 限制性酶识别错误？")
            print("     • 参考基因组不匹配？")
            print("     • 建库质量问题？")
    else:
        print("\n⚠ 无法计算比率")
        print("原始统计：")
        for key, value in stats.items():
            print(f"  {key}: {value:,}")
else:
    print("\n⚠ 未能提取统计数据")

print("\n" + "="*70 + "\n")
PYEOF

#==========================================
# 步骤9：生成比对统计摘要
#==========================================

echo ""
echo "=========================================="
echo "步骤9: 比对统计摘要"
echo "=========================================="

if [ -f "mapping/${SAMPLE}_R1.bowtie2.log" ]; then
    echo ""
    echo "Read 1 比对统计:"
    grep -E "reads|aligned" mapping/${SAMPLE}_R1.bowtie2.log | head -5
fi

if [ -f "mapping/${SAMPLE}_R2.bowtie2.log" ]; then
    echo ""
    echo "Read 2 比对统计:"
    grep -E "reads|aligned" mapping/${SAMPLE}_R2.bowtie2.log | head -5
fi

#==========================================
# 完成
#==========================================

echo ""
echo "=========================================="
echo "🎉 分析全部完成！"
echo "=========================================="
echo ""
echo "📁 重要输出文件:"
echo ""
echo "1. 质控报告:"
echo "   - HTML报告: reports/hicQC.html ⭐"
echo "   - 详细日志: qc_results/*_QC.log"
echo ""
echo "2. 可视化:"
echo "   - Contact Map: reports/contact_map.png"
echo "   - 诊断图: reports/diagnostic_plot.png"
echo ""
echo "3. 数据文件:"
echo "   - Hi-C矩阵: matrix/${SAMPLE}_${BIN_SIZE}.h5"
echo "   - BAM文件: mapping/${SAMPLE}_R1.bam, mapping/${SAMPLE}_R2.bam"
echo ""
echo "4. 日志:"
echo "   - 构建日志: matrix/build_matrix.log"
echo "   - 比对日志: mapping/${SAMPLE}_R*.bowtie2.log"
echo ""
echo "=========================================="
echo "下一步操作:"
echo "=========================================="
echo ""
echo "📊 查看质控报告:"
echo "   在浏览器中打开: reports/hicQC.html"
echo ""
echo "📧 如果是小测数据且质量OK:"
echo "   反馈测序公司安排大测"
echo ""
echo "🧬 如果是大测数据且质量OK:"
echo "   可以进行Hi-C辅助基因组组装（Juicer + 3D-DNA）"
echo ""
echo "=========================================="
echo ""