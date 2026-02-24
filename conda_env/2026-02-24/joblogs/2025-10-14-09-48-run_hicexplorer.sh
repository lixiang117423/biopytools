#!/bin/bash
#=========================================
# Juicer Hi-C质控完整流程
# 适合基因组组装辅助的Hi-C数据质控
#=========================================

set -e

#=========================================
# 配置参数
#=========================================

WORK_DIR="/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/62.hic/20251011_test"
cd ${WORK_DIR}

SAMPLE="ov53-1-HIC1"
READ1="E250928002_L01_ov53-1-HIC1_1.fq.gz"
READ2="E250928002_L01_ov53-1-HIC1_2.fq.gz"
GENOME_FASTA="OV53_1.primary.fasta"
GENOME_NAME="OV53_1"

# 限制性酶识别序列（MboI = GATC）
RESTRICTION_SITE="GATC"
THREADS=88

# Juicer路径（需要根据实际安装位置修改）
JUICER_DIR="/path/to/juicer"  # 修改为你的Juicer安装路径
# 如果conda安装，通常不需要设置JUICER_DIR

echo "=========================================="
echo "Juicer Hi-C 质控分析"
echo "=========================================="
echo "样品: ${SAMPLE}"
echo "限制性酶: MboI (${RESTRICTION_SITE})"
echo "线程数: ${THREADS}"
echo "工作目录: ${WORK_DIR}"
echo ""
date
echo ""

#=========================================
# 步骤1：创建Juicer标准目录结构
#=========================================

echo "步骤1: 创建Juicer目录结构..."

mkdir -p juicer_work
cd juicer_work

mkdir -p fastq references restriction_sites splits aligned HIC_tmp

echo "  ✓ 目录创建完成"

#=========================================
# 步骤2：准备参考基因组
#=========================================

echo ""
echo "步骤2: 准备参考基因组..."

# 复制基因组到references目录
cp ../${GENOME_FASTA} references/${GENOME_NAME}.fasta

# 生成染色体大小文件
if [ ! -f "references/${GENOME_NAME}.chrom.sizes" ]; then
    echo "  生成染色体大小文件..."
    samtools faidx references/${GENOME_NAME}.fasta
    cut -f1,2 references/${GENOME_NAME}.fasta.fai > references/${GENOME_NAME}.chrom.sizes
    echo "  ✓ 染色体大小文件已生成"
else
    echo "  ✓ 染色体大小文件已存在"
fi

# 显示染色体信息
echo ""
echo "  染色体信息:"
head -5 references/${GENOME_NAME}.chrom.sizes
TOTAL_SIZE=$(awk '{sum+=$2} END {printf "%.2f Mb", sum/1000000}' references/${GENOME_NAME}.chrom.sizes)
echo "  总长度: ${TOTAL_SIZE}"

#=========================================
# 步骤3：构建BWA索引
#=========================================

echo ""
echo "步骤3: 构建BWA索引..."

if [ ! -f "references/${GENOME_NAME}.fasta.bwt" ]; then
    echo "  构建索引中（这可能需要10-30分钟）..."
    bwa index references/${GENOME_NAME}.fasta
    echo "  ✓ BWA索引构建完成"
else
    echo "  ✓ BWA索引已存在"
fi

#=========================================
# 步骤4：生成限制性酶切位点文件
#=========================================

echo ""
echo "步骤4: 生成限制性酶切位点文件..."

SITE_FILE="restriction_sites/${GENOME_NAME}_${RESTRICTION_SITE}.txt"

if [ ! -f "${SITE_FILE}" ]; then
    echo "  查找 ${RESTRICTION_SITE} 位点..."
    
    # 使用python脚本查找酶切位点
    python3 << 'PYEOF'
import sys
from Bio import SeqIO
from Bio.Seq import Seq

genome_file = "references/OV53_1.fasta"
site = "GATC"
output_file = "restriction_sites/OV53_1_GATC.txt"

print(f"  扫描基因组查找 {site} 位点...")

with open(output_file, 'w') as out:
    for record in SeqIO.parse(genome_file, "fasta"):
        chrom = record.id
        seq = str(record.seq).upper()
        
        # 查找正向位点
        pos = 0
        count = 0
        while True:
            pos = seq.find(site, pos)
            if pos == -1:
                break
            out.write(f"{chrom} {pos}\n")
            count += 1
            pos += 1
        
        print(f"    {chrom}: 找到 {count} 个位点")

print(f"  ✓ 位点文件已保存: {output_file}")
PYEOF
    
    echo "  ✓ 限制性酶切位点文件已生成"
else
    echo "  ✓ 限制性酶切位点文件已存在"
fi

# 统计位点数量
SITE_COUNT=$(wc -l < ${SITE_FILE})
echo "  总位点数: ${SITE_COUNT}"

#=========================================
# 步骤5：准备fastq文件
#=========================================

echo ""
echo "步骤5: 准备fastq文件..."

# Juicer要求特定的文件命名格式
ln -sf ../../${READ1} fastq/${SAMPLE}_R1.fastq.gz
ln -sf ../../${READ2} fastq/${SAMPLE}_R2.fastq.gz

echo "  ✓ fastq文件链接完成"

#=========================================
# 步骤6：运行Juicer Pipeline
#=========================================

echo ""
echo "=========================================="
echo "步骤6: 运行Juicer Pipeline"
echo "=========================================="
echo "  这将需要 1-4 小时，请耐心等待..."
echo ""

# Juicer运行命令
# 注意：根据你的Juicer安装方式，命令可能略有不同

if command -v juicer.sh &> /dev/null; then
    # 如果juicer.sh在PATH中
    juicer.sh \
        -g ${GENOME_NAME} \
        -s ${RESTRICTION_SITE} \
        -z references/${GENOME_NAME}.fasta \
        -p references/${GENOME_NAME}.chrom.sizes \
        -y ${SITE_FILE} \
        -d $(pwd) \
        -D $(pwd) \
        -t ${THREADS} \
        2>&1 | tee juicer_run.log
else
    # 如果需要指定Juicer路径
    echo "  请确保已安装Juicer，或修改脚本中的JUICER_DIR路径"
    echo "  下载Juicer: git clone https://github.com/aidenlab/juicer.git"
    exit 1
fi

#=========================================
# 步骤7：提取质控统计
#=========================================

echo ""
echo "=========================================="
echo "步骤7: Hi-C质控报告"
echo "=========================================="

# Juicer会生成aligned/inter.txt等统计文件
if [ -f "aligned/inter.txt" ]; then
    python3 << 'PYEOF'
import os
import re

print("\n" + "="*70)
print(" "*20 + "Juicer 质控报告")
print("="*70)

# 读取统计文件
stats_file = "aligned/inter.txt"
if os.path.exists(stats_file):
    with open(stats_file) as f:
        lines = f.readlines()
    
    # 解析统计信息
    for line in lines:
        line = line.strip()
        if line.startswith("Sequenced Read Pairs:"):
            total_pairs = int(line.split(":")[1].strip())
            print(f"\n总read pairs: {total_pairs:,}")
        elif line.startswith("Normal Paired:"):
            normal = int(line.split(":")[1].strip().split()[0])
            pct = float(line.split("(")[1].split("%")[0])
            print(f"正常配对:     {normal:,}  ({pct:.2f}%)")
        elif line.startswith("Chimeric Paired:"):
            chimeric = int(line.split(":")[1].strip().split()[0])
            pct = float(line.split("(")[1].split("%")[0])
            print(f"嵌合配对:     {chimeric:,}  ({pct:.2f}%)")
        elif line.startswith("Unmapped:"):
            unmapped = int(line.split(":")[1].strip().split()[0])
            pct = float(line.split("(")[1].split("%")[0])
            print(f"未比对:       {unmapped:,}  ({pct:.2f}%)")

# 读取inter_30.txt（更详细的统计）
inter30_file = "aligned/inter_30.txt"
if os.path.exists(inter30_file):
    with open(inter30_file) as f:
        content = f.read()
    
    print("\n" + "="*70)
    print("详细质量指标")
    print("="*70)
    
    # 提取关键指标
    patterns = {
        'Unique Valid Pairs': r'Unique.*?Valid.*?:\s*([\d,]+)',
        'PCR Duplicates': r'PCR.*?[Dd]uplicate.*?:\s*([\d,]+)',
        'Intra-chromosomal': r'Intra.*?chromosomal.*?:\s*([\d,]+)',
        'Inter-chromosomal': r'Inter.*?chromosomal.*?:\s*([\d,]+)',
        'Short Range (<20Kb)': r'[Ss]hort.*?[Rr]ange.*?:\s*([\d,]+)',
        'Long Range (>20Kb)': r'[Ll]ong.*?[Rr]ange.*?:\s*([\d,]+)',
    }
    
    for key, pattern in patterns.items():
        match = re.search(pattern, content, re.IGNORECASE)
        if match:
            value = int(match.group(1).replace(',', ''))
            print(f"\n{key:.<30} {value:>15,}")

print("\n" + "="*70)
print("质量评估")
print("="*70)

# 读取基本统计
if os.path.exists(stats_file):
    with open(stats_file) as f:
        content = f.read()
    
    # 计算质量指标
    total_match = re.search(r'Sequenced Read Pairs:\s*([\d,]+)', content)
    normal_match = re.search(r'Normal Paired:\s*([\d,]+)', content)
    
    if total_match and normal_match:
        total = int(total_match.group(1).replace(',', ''))
        normal = int(normal_match.group(1).replace(',', ''))
        
        if total > 0:
            normal_rate = normal / total * 100
            print(f"\n正常配对率: {normal_rate:.2f}%")
            
            if normal_rate >= 80:
                print("  ✓✓✓ 优秀！")
                result = "优秀"
            elif normal_rate >= 70:
                print("  ✓✓ 良好！")
                result = "良好"
            elif normal_rate >= 60:
                print("  ✓ 合格")
                result = "合格"
            else:
                print("  ⚠ 偏低")
                result = "偏低"

print("\n" + "="*70)
print("结论")
print("="*70)

if os.path.exists(inter30_file):
    with open(inter30_file) as f:
        content = f.read()
    
    # 查找有效配对数
    valid_match = re.search(r'Unique.*?Valid.*?:\s*([\d,]+)', content, re.IGNORECASE)
    if valid_match:
        valid_pairs = int(valid_match.group(1).replace(',', ''))
        
        if valid_pairs >= 50_000_000:
            print("\n✅ 数据量充足，质量优秀！")
            print("✅ 适合用于Hi-C辅助基因组组装")
        elif valid_pairs >= 10_000_000:
            print("\n✅ 数据质量良好，可用于组装辅助")
        elif valid_pairs >= 1_000_000:
            print("\n⚠️ 数据量较少，建议增加测序深度")
        else:
            print("\n❌ 数据量不足，建议重新测序")

print("\n" + "="*70 + "\n")
PYEOF

else
    echo "  ⚠ 未找到Juicer统计文件"
    echo "  检查 juicer_run.log 查看详细信息"
fi

#=========================================
# 步骤8：生成.hic文件（用于可视化）
#=========================================

echo ""
echo "步骤8: 生成.hic文件（可选）..."

if [ -f "aligned/merged_nodups.txt" ]; then
    echo "  生成.hic文件用于Juicebox可视化..."
    
    # 需要juicer_tools
    if command -v juicer_tools &> /dev/null; then
        juicer_tools pre \
            aligned/merged_nodups.txt \
            ${SAMPLE}.hic \
            references/${GENOME_NAME}.chrom.sizes \
            2>&1 | tee hic_generation.log
        
        if [ -f "${SAMPLE}.hic" ]; then
            HIC_SIZE=$(du -h ${SAMPLE}.hic | cut -f1)
            echo "  ✓ .hic文件已生成: ${SAMPLE}.hic (${HIC_SIZE})"
            echo "  可使用Juicebox打开查看: https://aidenlab.org/juicebox/"
        fi
    else
        echo "  ⚠ 需要安装juicer_tools来生成.hic文件"
        echo "  下载: https://github.com/aidenlab/juicer/wiki/Download"
    fi
else
    echo "  ⚠ 未找到merged_nodups.txt文件"
fi

#=========================================
# 完成
#=========================================

echo ""
echo "=========================================="
echo "🎉 Juicer分析完成！"
echo "=========================================="
echo ""
echo "重要输出文件:"
echo "  📊 主统计文件: aligned/inter.txt"
echo "  📈 详细统计: aligned/inter_30.txt"
echo "  💾 去重配对: aligned/merged_nodups.txt"
if [ -f "${SAMPLE}.hic" ]; then
    echo "  🖼️  可视化文件: ${SAMPLE}.hic"
fi
echo ""
echo "日志文件:"
echo "  - Juicer运行日志: juicer_run.log"
echo "  - 比对日志: aligned/*.out"
echo ""
echo "下一步（基因组组装）:"
echo "  1. 使用3D-DNA进行scaffolding"
echo "  2. 使用Juicebox进行手动校正"
echo ""
date
echo ""
echo "=========================================="