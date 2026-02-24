#!/bin/bash
# =====================================================
# 🛠️ JBAT .hic 修复工具
# 功能: 利用现有的 out_JBAT.txt 重新生成可用的 .hic 文件
# =====================================================

# --- 1. 设置环境 (请根据你的实际情况修改路径) ---
JUICER_JAR="/share/org/YZWL/yzwl_lixg/software/juicer/scripts/common/juicer_tools.jar"
# 原始的 Contig 索引文件 (用于万一 log 里找不到大小时手动计算)
REF_FAI="/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/04.YaHS/OV53_1.primary.fa.fai"

# --- 2. 检查输入文件 ---
if [ ! -s "out_JBAT.txt" ]; then
    echo "❌ 错误: out_JBAT.txt 不存在或为空！无法修复。"
    exit 1
fi

echo "✅ 发现输入文件: out_JBAT.txt ($(du -h out_JBAT.txt | cut -f1))"

# --- 3. 生成正确的 chrom.sizes 文件 ---
echo "🔍 正在获取基因组总长度..."

# 方法 A: 尝试从 YaHS 日志中获取 (最准确)
if [ -f "out_JBAT.log" ] && grep -q "PRE_C_SIZE" out_JBAT.log; then
    grep "PRE_C_SIZE" out_JBAT.log | awk '{print $2" "$3}' > jbat.sizes.fix
    echo "   > 从日志中提取成功: $(cat jbat.sizes.fix)"
else
    # 方法 B: 如果日志没提取到，手动计算 (求和所有 contig 长度)
    echo "   ⚠️ 日志中未找到长度，正在根据 FAI 索引计算..."
    if [ -f "${REF_FAI}" ]; then
        awk '{sum+=$2} END {print "assembly", sum}' ${REF_FAI} > jbat.sizes.fix
        echo "   > 手动计算成功: $(cat jbat.sizes.fix)"
    else
        echo "❌ 错误: 找不到 FAI 索引文件，无法计算长度。"
        exit 1
    fi
fi

# --- 4. 检查 out_JBAT.txt 的格式 ---
# 确保第一行包含 'assembly' (防止意外)
# 读取前几行检查
echo "🔍 检查文本文件格式..."
HEAD_CHECK=$(head -n 3 out_JBAT.txt | grep "assembly")
if [ -z "$HEAD_CHECK" ]; then
    echo "⚠️  警告: out_JBAT.txt 的前几行似乎不包含 'assembly' 关键字。"
    echo "这可能导致生成为空。程序将继续尝试，但请留意结果。"
fi

# --- 5. 运行 Juicer Tools ---
echo "🚀 正在重新生成 .hic 文件 (这可能需要几分钟到几十分钟)..."

# 这里的关键是:
# 1. 内存给足 (-Xmx)
# 2. 确保 size 文件只有一行: assembly <长度>
# 3. 输入文件 42G 很大，不要用排序管道，直接喂给它

java -Xmx128G -jar ${JUICER_JAR} pre \
    out_JBAT.txt \
    out_JBAT_fixed.hic \
    jbat.sizes.fix

# --- 6. 结果检查 ---
if [ -s "out_JBAT_fixed.hic" ]; then
    FILE_SIZE=$(stat -c%s "out_JBAT_fixed.hic")
    if [ $FILE_SIZE -gt 1000000 ]; then
        echo "🎉 修复成功！"
        echo "生成的有效文件: out_JBAT_fixed.hic"
        echo "大小: $(du -h out_JBAT_fixed.hic | cut -f1)"
        echo "👉 请下载这个文件和 out_JBAT.assembly 进行手动纠错。"
    else
        echo "❌ 修复失败: 文件生成了，但依然很小 ($FILE_SIZE bytes)。"
        echo "可能是 Java 内存溢出或输入文件格式被破坏。"
    fi
else
    echo "❌ 修复失败: 未生成文件。"
fi
