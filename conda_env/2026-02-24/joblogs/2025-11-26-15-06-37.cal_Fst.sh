#!/bin/bash

# ==============================================================================
# 🧬 自动化 Fst 计算脚本 (Python 强力清洗版)
# ==============================================================================

# 1️⃣ 定义路径 (请确认路径无误)
INPUT_VCF="/share/org/YZWL/yzwl_lixg/project/19.大豆疫霉菌/13.Fst/variation.filtered.snp.vcf.gz"
INPUT_SETS="/share/org/YZWL/yzwl_lixg/project/19.大豆疫霉菌/13.Fst/55.D检验的分组信息.txt"
OUT_DIR="/share/org/YZWL/yzwl_lixg/project/19.大豆疫霉菌/13.Fst/Fst_Output"

# 创建目录
mkdir -p "$OUT_DIR"
mkdir -p "${OUT_DIR}/pop_lists"
mkdir -p "${OUT_DIR}/logs"

# 定义清洗后的文件路径
CLEAN_SETS="${OUT_DIR}/clean_sets_final.txt"

# ==============================================================================
# 2️⃣ 使用 Python 进行智能格式清洗 (解决编码和分隔符问题)
# ------------------------------------------------------------------------------
echo "🚀 开始处理..."
echo "🧹 正在调用 Python 进行文件清洗和转码..."

python3 -c "
import sys

input_file = '$INPUT_SETS'
output_file = '$CLEAN_SETS'

# 尝试读取文件，处理编码问题
lines = []
try:
    # 先尝试 UTF-8
    with open(input_file, 'r', encoding='utf-8') as f:
        lines = f.readlines()
except UnicodeDecodeError:
    print('⚠️  检测到非 UTF-8 编码，尝试使用 GB18030 (Windows中文) 读取...')
    try:
        # 如果失败，尝试 GB18030 (兼容 GBK)
        with open(input_file, 'r', encoding='gb18030') as f:
            lines = f.readlines()
    except Exception as e:
        print(f'❌ 读取文件失败: {e}')
        sys.exit(1)

valid_count = 0
groups = set()

with open(output_file, 'w', encoding='utf-8') as f_out:
    for line in lines:
        # 去除首尾空白
        line = line.strip()
        if not line: continue
        
        # 按空白符(空格或Tab)分割
        parts = line.split()
        
        # 必须至少有两列 (SampleID, GroupID)
        if len(parts) >= 2:
            sample = parts[0]
            group = parts[1]
            
            # 过滤掉可能的脏数据 (例如 Group 名不能纯粹是数字，除非你确认它是)
            # 这里我们只做标准化输出: Sample [TAB] Group
            f_out.write(f'{sample}\t{group}\n')
            groups.add(group)
            valid_count += 1
        else:
            print(f'⚠️  跳过无效行 (列数不足): {line}')

print(f'✅ 清洗完成！保留了 {valid_count} 个样本。')
print(f'📋 识别到的分组 ({len(groups)} 个): {sorted(list(groups))}')
"

# 检查 Python 是否运行成功
if [ $? -ne 0 ]; then
    echo "❌ Python 清洗脚本运行失败，请检查报错。"
    exit 1
fi

echo "📄 清洗后文件预览 ($CLEAN_SETS):"
head -n 5 "$CLEAN_SETS"
echo "------------------------------------------------------------------"

# ==============================================================================
# 3️⃣ 拆分群体并过滤样本数少于2的群体 (改用 Python 处理)
# ------------------------------------------------------------------------------
echo "🧩 正在拆分群体列表并过滤..."

# 定义群体列表文件
GROUPS_FILE="${OUT_DIR}/valid_groups.txt"

python3 -c "
import os

clean_file = '$CLEAN_SETS'
pop_list_dir = '${OUT_DIR}/pop_lists'
groups_file = '$GROUPS_FILE'
min_samples = 2  # 最小样本数阈值

# 读取数据
groups = {}
with open(clean_file, 'r', encoding='utf-8') as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 2:
            sample, group = parts[0], parts[1]
            if group not in groups:
                groups[group] = []
            groups[group].append(sample)

# 过滤并写入各群体的样本列表
valid_groups = []
excluded_groups = []

for group, samples in sorted(groups.items()):
    sample_count = len(samples)
    
    if sample_count >= min_samples:
        # 样本数足够，写入文件
        target_file = os.path.join(pop_list_dir, f'{group}.txt')
        with open(target_file, 'w', encoding='utf-8') as f:
            for sample in samples:
                f.write(f'{sample}\n')
        print(f'   ✓ 群体: [{group}] (样本数: {sample_count})')
        valid_groups.append(group)
    else:
        # 样本数不足，剔除
        print(f'   ✗ 群体: [{group}] (样本数: {sample_count}) - 剔除 (少于 {min_samples} 个样本)')
        excluded_groups.append(group)

# 将有效群体名写入文件
with open(groups_file, 'w', encoding='utf-8') as f:
    for group in valid_groups:
        f.write(f'{group}\n')

print(f'')
print(f'✅ 有效群体数量: {len(valid_groups)}')
if excluded_groups:
    print(f'⚠️  已剔除 {len(excluded_groups)} 个样本数不足的群体: {excluded_groups}')
"

# 检查 Python 是否运行成功
if [ $? -ne 0 ]; then
    echo "❌ Python 拆分脚本运行失败，请检查报错。"
    exit 1
fi

# 从文件读取群体列表到数组
count=0
GROUP_ARRAY=()
while IFS= read -r group; do
    if [ -n "$group" ]; then
        GROUP_ARRAY[$count]="$group"
        let count++
    fi
done < "$GROUPS_FILE"

echo ""
echo "📝 已加载 $count 个有效群体到数组"

if [ "$count" -lt 2 ]; then
    echo "❌ 错误: 有效群体少于 2 个，无法进行 Fst 计算。"
    echo "   需要至少 2 个群体，且每个群体至少 2 个样本。"
    exit 1
fi

# 调试：显示数组内容
echo "🔍 群体数组内容："
for (( i=0; i<count; i++ )); do
    echo "   [$i] ${GROUP_ARRAY[$i]}"
done

# ==============================================================================
# 4️⃣ 计算 Fst
# ------------------------------------------------------------------------------
echo "------------------------------------------------------------------"
echo "⚔️  开始计算 Fst (共 $((count*(count-1)/2)) 个配对)..."
SUMMARY_FILE="${OUT_DIR}/Final_Fst_Summary.txt"
echo -e "Group1\tGroup2\tMean_Fst\tWeighted_Fst" > "$SUMMARY_FILE"

pair_count=0
for (( i=0; i<count; i++ )); do
    for (( j=i+1; j<count; j++ )); do
        POP1="${GROUP_ARRAY[$i]}"
        POP2="${GROUP_ARRAY[$j]}"
        OUT_PREFIX="${OUT_DIR}/logs/Fst_${i}_vs_${j}"
        
        let pair_count++
        echo "   👉 [$pair_count] ${POP1} vs ${POP2}"
        
        vcftools --gzvcf "$INPUT_VCF" \
            --weir-fst-pop "${OUT_DIR}/pop_lists/${POP1}.txt" \
            --weir-fst-pop "${OUT_DIR}/pop_lists/${POP2}.txt" \
            --out "$OUT_PREFIX" \
            --remove-indels --min-alleles 2 --max-alleles 2 \
            > "${OUT_PREFIX}.log" 2>&1
        
        MEAN=$(grep "Weir and Cockerham mean Fst estimate" "${OUT_PREFIX}.log" | awk '{print $7}')
        WEIGHTED=$(grep "Weir and Cockerham weighted Fst estimate" "${OUT_PREFIX}.log" | awk '{print $7}')
        
        # 如果为空给个默认值
        : ${MEAN:="NaN"}
        : ${WEIGHTED:="NaN"}
        
        echo -e "${POP1}\t${POP2}\t${MEAN}\t${WEIGHTED}" >> "$SUMMARY_FILE"
    done
done

echo "------------------------------------------------------------------"
echo "🎉 计算完成！"
echo "📊 结果文件: $SUMMARY_FILE"
echo "📁 详细日志: ${OUT_DIR}/logs/"