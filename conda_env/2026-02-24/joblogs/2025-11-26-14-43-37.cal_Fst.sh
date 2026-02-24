#!/bin/bash

# ==============================================================================
# 🧬 自动化 Fst 计算脚本 (基于 VCFtools)
# ==============================================================================

# 1️⃣ 定义输入输出路径
# ------------------------------------------------------------------------------
# 输入文件
INPUT_VCF="/share/org/YZWL/yzwl_lixg/project/19.大豆疫霉菌/13.Fst/variation.filtered.snp.vcf.gz"
INPUT_SETS="/share/org/YZWL/yzwl_lixg/project/19.大豆疫霉菌/13.Fst/55.D检验的分组信息.txt"

# 输出主目录
OUT_DIR="/share/org/YZWL/yzwl_lixg/project/19.大豆疫霉菌/13.Fst/Fst_Output"

# 临时存放群体样本列表的文件夹
POP_LIST_DIR="${OUT_DIR}/pop_lists"
# 存放 VCFtools 原始日志的文件夹
LOG_DIR="${OUT_DIR}/logs"

# ==============================================================================
# 2️⃣ 环境检查
# ------------------------------------------------------------------------------
echo "🚀 开始 Fst 分析流程..."

# 检查 vcftools
if ! command -v vcftools &> /dev/null; then
    echo "❌ 错误: 未找到 vcftools 命令。请先安装或加载它 (module load vcftools)！"
    exit 1
fi

# 创建目录
mkdir -p "$OUT_DIR"
mkdir -p "$POP_LIST_DIR"
mkdir -p "$LOG_DIR"

echo "📄 VCF 文件: $INPUT_VCF"
echo "📂 输出目录: $OUT_DIR"

# ==============================================================================
# 3️⃣ 自动化拆分群体 (生成样本列表)
# ------------------------------------------------------------------------------
echo "------------------------------------------------------------------"
echo "🧩 第一步: 正在从分组文件中拆分群体列表..."

# 1. 提取所有不重复的群体名称 (排除 xxx 标记的个体)
# 假设第二列是群体名
GROUPS=$(awk '$2 != "xxx" {print $2}' "$INPUT_SETS" | sort | uniq)

# 2. 为每个群体生成单独的 .txt 文件
count=0
for group in $GROUPS; do
    # 提取属于该 group 的第一列(SampleID)，写入文件
    awk -v g="$group" '$2 == g {print $1}' "$INPUT_SETS" > "${POP_LIST_DIR}/${group}.txt"
    echo "   - 已生成群体列表: ${group} (包含 $(wc -l < "${POP_LIST_DIR}/${group}.txt") 个样本)"
    
    # 将群体名存入数组，方便后面做循环
    GROUP_ARRAY[$count]=$group
    let count++
done

echo "✅ 群体列表拆分完成！共发现 $count 个有效群体。"

# ==============================================================================
# 4️⃣ 循环计算 Pairwise Fst
# ------------------------------------------------------------------------------
echo "------------------------------------------------------------------"
echo "⚔️  第二步: 开始两两计算 Fst (基于 Weir and Cockerham 1984)..."

# 初始化汇总结果文件
SUMMARY_FILE="${OUT_DIR}/Final_Fst_Summary.txt"
echo -e "Group1\tGroup2\tMean_Fst\tWeighted_Fst" > "$SUMMARY_FILE"

# 双重循环遍历所有组合
for (( i=0; i<count; i++ )); do
    for (( j=i+1; j<count; j++ )); do
        
        POP1=${GROUP_ARRAY[$i]}
        POP2=${GROUP_ARRAY[$j]}
        
        # 定义输出前缀
        OUT_PREFIX="${LOG_DIR}/${POP1}_vs_${POP2}"
        
        echo "   👉 正在计算: $POP1 vs $POP2 ..."
        
        # 运行 VCFtools
        # --weir-fst-pop: 指定群体文件
        # --fst-window-size: 也可以指定窗口，这里我们先看不加窗口的位点Fst，并取均值
        # 2>&1 | tee : 将屏幕输出同时保存到 log 文件
        
        vcftools --gzvcf "$INPUT_VCF" \
            --weir-fst-pop "${POP_LIST_DIR}/${POP1}.txt" \
            --weir-fst-pop "${POP_LIST_DIR}/${POP2}.txt" \
            --out "$OUT_PREFIX" \
            --remove-indels --min-alleles 2 --max-alleles 2 \
            > "${OUT_PREFIX}.log" 2>&1
        
        # 提取 Mean Fst 值 (从 Log 文件中 grep 出来)
        # VCFtools 的日志里有一行写着: "Weir and Cockerham mean Fst estimate: 0.12345"
        MEAN_FST=$(grep "Weir and Cockerham mean Fst estimate" "${OUT_PREFIX}.log" | awk '{print $7}')
        
        # 也是从 log 中提取加权 Fst (Weighted Fst)，有时这个更准
        WEIGHTED_FST=$(grep "Weir and Cockerham weighted Fst estimate" "${OUT_PREFIX}.log" | awk '{print $7}')
        
        # 如果计算失败（比如没有共有的SNP），Fst可能是 NaN
        if [ -z "$MEAN_FST" ]; then MEAN_FST="NaN"; fi
        if [ -z "$WEIGHTED_FST" ]; then WEIGHTED_FST="NaN"; fi
        
        # 写入汇总表
        echo -e "${POP1}\t${POP2}\t${MEAN_FST}\t${WEIGHTED_FST}" >> "$SUMMARY_FILE"
        
    done
done

# ==============================================================================
# 5️⃣ 结束
# ------------------------------------------------------------------------------
echo "------------------------------------------------------------------"
echo "🎉 Fst 计算全部完成！"
echo "📊 最终汇总结果请查看: $SUMMARY_FILE"
echo "📂 详细 Log 和 单位点 Fst 数据在: $LOG_DIR"
echo "------------------------------------------------------------------"