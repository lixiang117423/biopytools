#!/bin/bash

# ============================================
# 📊 BAM 文件覆盖度和比对率统计脚本
# ============================================

# --- 📁 配置路径 ---
# PanDepth 软件的完整路径
PANDEPTH_PATH="/share/org/YZWL/yzwl_lixg/.local/bin/pandepth"
# 存放 BAM 文件的文件夹
BAM_DIR="/share/org/YZWL/yzwl_lixg/project/19.大豆疫霉菌/02.mapping/bam"
# 输出结果的文件夹
OUTPUT_DIR="/share/org/YZWL/yzwl_lixg/project/19.大豆疫霉菌/05.比对结果统计"

# --- ✅ 检查工具和路径 ---
echo "🔍 正在检查运行环境..."

# 检查 PanDepth 是否可执行
if [ ! -x "$PANDEPTH_PATH" ]; then
    echo "❌ 错误：PanDepth 不存在或没有执行权限: $PANDEPTH_PATH"
    exit 1
fi
echo "✅ PanDepth 工具检查通过"

# 检查 samtools 是否在环境中
if ! command -v samtools &> /dev/null; then
    echo "❌ 错误：samtools 未安装或未在环境路径中。"
    echo "💡 请先安装 samtools (例如: conda install samtools)。"
    exit 1
fi
echo "✅ samtools 工具检查通过"

# 检查 BAM 文件夹是否存在
if [ ! -d "$BAM_DIR" ]; then
    echo "❌ 错误：BAM 文件夹不存在: $BAM_DIR"
    exit 1
fi
echo "✅ BAM 文件夹存在: $BAM_DIR"

# 统计 BAM 文件数量
bam_count=$(find "$BAM_DIR" -maxdepth 1 -name "*.bam" -type f | wc -l)
if [ "$bam_count" -eq 0 ]; then
    echo "⚠️  警告：在 $BAM_DIR 中未找到任何 BAM 文件"
    exit 1
fi
echo "📦 找到 $bam_count 个 BAM 文件待处理"

# 如果输出文件夹不存在，则创建它
mkdir -p "$OUTPUT_DIR"
echo "📂 输出目录已准备: $OUTPUT_DIR"

# --- 🚀 主分析流程 ---
echo ""
echo "============================================"
echo "🚀 开始处理 BAM 文件..."
echo "============================================"
echo ""

# 创建一个总的统计结果文件
SUMMARY_FILE="$OUTPUT_DIR/summary_statistics.txt"
echo -e "Sample\tAverage_Coverage\tAlignment_Rate(%)" > "$SUMMARY_FILE"

# 计数器
processed=0
failed=0

# 遍历所有 bam 文件
for bam_file in "$BAM_DIR"/*.bam; do
    if [ -f "$bam_file" ]; then
        # 从文件路径中提取样本名 (例如: test.bam -> test)
        sample_name=$(basename "$bam_file" .bam)
        
        echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
        echo "🧬 正在处理样本: $sample_name"
        echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
        
        # --- 1. 使用 PanDepth 统计覆盖度 ---
        echo "  📈 [1/2] 正在使用 PanDepth 计算覆盖度..."
        output_prefix="$OUTPUT_DIR/$sample_name"
        
        avg_coverage="N/A"
        if $PANDEPTH_PATH -i "$bam_file" -o "$output_prefix" 2>/dev/null; then
            echo "      ✅ PanDepth 覆盖度统计完成"
            echo "      💾 结果保存: ${output_prefix}.chr.stat.gz"
            
            # 提取平均覆盖度（从 .chr.stat.gz 文件的第三列计算均值）
            if [ -f "${output_prefix}.chr.stat.gz" ]; then
                avg_coverage=$(zcat "${output_prefix}.chr.stat.gz" | awk 'NR>1 {sum+=$3; count++} END {if(count>0) printf "%.2f", sum/count; else print "N/A"}')
                echo "      📊 平均覆盖度: ${avg_coverage}X"
            fi
        else
            echo "      ⚠️  PanDepth 运行出现问题"
            ((failed++))
        fi
        
        # --- 2. 使用 samtools flagstat 统计比对率 ---
        echo "  📊 [2/2] 正在使用 samtools flagstat 统计比对率..."
        flagstat_output_file="$OUTPUT_DIR/${sample_name}.flagstat.txt"
        
        mapped_rate="N/A"
        if samtools flagstat "$bam_file" > "$flagstat_output_file" 2>/dev/null; then
            # 从 flagstat 结果中提取整体比对率
            # 'mapped' 行的格式通常是: "X + Y mapped (Z% : N/A)"
            mapped_rate=$(grep "mapped (" "$flagstat_output_file" | head -n 1 | awk -F '[()%]' '{print $2}')
            
            if [ -n "$mapped_rate" ]; then
                echo "      ✅ 比对率统计完成: ${mapped_rate}%"
            else
                echo "      ⚠️  无法从 flagstat 输出中提取比对率"
                mapped_rate="N/A"
            fi
        else
            echo "      ❌ samtools flagstat 运行失败"
            ((failed++))
        fi
        
        # 将样本名、平均覆盖度和比对率追加到汇总文件中
        echo -e "${sample_name}\t${avg_coverage}\t${mapped_rate}" >> "$SUMMARY_FILE"
        
        ((processed++))
        echo ""
    fi
done

# --- 📋 总结报告 ---
echo "============================================"
echo "✨ 所有样本处理完毕！"
echo "============================================"
echo ""
echo "📊 处理统计："
echo "   ✅ 成功处理: $processed 个样本"
if [ "$failed" -gt 0 ]; then
    echo "   ⚠️  部分失败: $failed 个样本存在问题"
fi
echo ""
echo "📁 输出文件："
echo "   📈 覆盖度详细结果: $OUTPUT_DIR/*.chr.stat.gz"
echo "   📊 比对率统计结果: $OUTPUT_DIR/*.flagstat.txt"
echo "   📋 汇总统计文件: $SUMMARY_FILE"
echo ""
echo "🎉 分析完成！"