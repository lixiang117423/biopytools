#!/bin/bash
# =====================================================
# 🚑 JBAT 越界强制修复工具
# 原理: 扫描文本文件找到最大坐标，确保 chrom.sizes 绝对够大
# =====================================================

JUICER_JAR="/share/org/YZWL/yzwl_lixg/software/juicer/scripts/common/juicer_tools.jar"
INPUT_TXT="out_JBAT.txt"
OUTPUT_HIC="out_JBAT_final.hic"

# --- 1. 检查输入 ---
if [ ! -f "$INPUT_TXT" ]; then
    echo "❌ 找不到 $INPUT_TXT"; exit 1;
fi

echo "🕵️‍♂️ 正在扫描 42G 文件以寻找最大坐标 (这需要几分钟)..."

# --- 2. 找出文件中的最大坐标 ---
# 这是一个极其暴力的做法：读取第3列和第7列，找出最大值
# 42G文件可能需要跑 5-10 分钟，请耐心等待
MAX_POS=$(awk '
    BEGIN { max = 0 }
    {
        if ($3 > max) max = $3
        if ($7 > max) max = $7
    }
    END { print max }
' $INPUT_TXT)

if [ -z "$MAX_POS" ] || [ "$MAX_POS" -eq 0 ]; then
    echo "❌ 错误: 无法从文件中读取有效坐标。"
    exit 1
fi

echo "   > 监测到的最大坐标是: ${MAX_POS}"

# --- 3. 增加安全余量 (Buffer) ---
# 给它加 1000 bp，确保绝对不越界
SAFE_SIZE=$(($MAX_POS + 1000))
echo "   > 设置的安全边界是: ${SAFE_SIZE}"

# --- 4. 生成新的 sizes 文件 ---
echo "assembly ${SAFE_SIZE}" > jbat_safe.sizes
echo "📝 生成了尺寸文件: jbat_safe.sizes"

# --- 5. 运行 Juicer ---
echo "🚀 正在重新生成 .hic..."

java -Xmx120G -Xms32G -jar ${JUICER_JAR} pre \
    $INPUT_TXT \
    $OUTPUT_HIC \
    jbat_safe.sizes

# --- 6. 验证 ---
if [ -s "$OUTPUT_HIC" ]; then
    SIZE=$(stat -c%s "$OUTPUT_HIC")
    if [ $SIZE -gt 1000000 ]; then
        echo "🎉🎉🎉 终于成功了！"
        echo "文件: $OUTPUT_HIC"
        echo "大小: $(du -h $OUTPUT_HIC | cut -f1)"
    else
        echo "💀 依然很小... 请检查 Java 版本或尝试更换 Juicer Tools jar 包版本。"
    fi
else
    echo "❌ 生成失败。"
fi
