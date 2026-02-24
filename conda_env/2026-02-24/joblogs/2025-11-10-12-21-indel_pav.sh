#!/bin/bash

# --- 设置变量 ---

# 输入的 Indel VCF 文件
VCF_FILE="variation.filtered.indel.vcf.gz"

# 输出的 0/1 矩阵文件名
OUTPUT_MATRIX="indel_pa_matrix_min50perc.01.tsv"

# 定义群体前缀
POP1_PREFIX="BS"
POP2_PREFIX="GC"

# --- 脚本开始 ---
echo "--- 生成带群体频率过滤的 0/1 矩阵 ---"
echo "输入 VCF 文件: ${VCF_FILE}"
echo "筛选标准: Indel 必须在 ${POP1_PREFIX} 或 ${POP2_PREFIX} 组中至少 50% 的样本里存在。"

# --- 步骤 1: 准备样本信息和计算阈值 ---
echo "步骤 1: 分析群体构成并计算阈值..."

# 创建临时文件来存储样本列表
bcftools query -l ${VCF_FILE} > samples.all.tmp
grep "^${POP1_PREFIX}" samples.all.tmp > samples.pop1.tmp
grep "^${POP2_PREFIX}" samples.all.tmp > samples.pop2.tmp

# 计算每个群体的样本数
POP1_COUNT=$(wc -l < samples.pop1.tmp)
POP2_COUNT=$(wc -l < samples.pop2.tmp)

# 检查样本数是否为0
if [ "$POP1_COUNT" -eq 0 ] || [ "$POP2_COUNT" -eq 0 ]; then
    echo "错误：一个或两个群体没有找到任何样本。请检查您的 VCF 文件和群体前缀。"
    rm -f *.tmp
    exit 1
fi

# 计算阈值（至少50%）。使用 awk 来确保正确处理奇偶数。
# 例如，10个样本的50%是5；9个样本的50%向上取整也是5。
POP1_THRESHOLD=$(awk -v n=${POP1_COUNT} 'BEGIN{print int(n/2 + 0.99)}')
POP2_THRESHOLD=$(awk -v n=${POP2_COUNT} 'BEGIN{print int(n/2 + 0.99)}')

echo "群体 ${POP1_PREFIX}: ${POP1_COUNT} 个样本, 保留阈值 = ${POP1_THRESHOLD}"
echo "群体 ${POP2_PREFIX}: ${POP2_COUNT} 个样本, 保留阈值 = ${POP2_THRESHOLD}"


# --- 步骤 2: 创建矩阵表头 ---
echo "步骤 2: 创建输出文件表头..."
# 将表头写入新文件（覆盖旧文件）
(echo -n -e "CHROM\tPOS\tREF\tALT\t" && \
 cat samples.all.tmp | tr '\n' '\t' | sed 's/\t$//g'
) > ${OUTPUT_MATRIX}


# --- 步骤 3: 处理 VCF 并生成过滤后的矩阵 ---
echo "步骤 3: 正在处理 VCF 文件，请稍候..."
# 将所有样本名作为一个空格分隔的字符串传递给 awk
ALL_SAMPLES_STR=$(cat samples.all.tmp | tr '\n' ' ')

bcftools query -f '\n%CHROM\t%POS\t%REF\t%ALT[\t%GT]' ${VCF_FILE} | \
awk -v p1_prefix="${POP1_PREFIX}" \
    -v p2_prefix="${POP2_PREFIX}" \
    -v p1_thresh="${POP1_THRESHOLD}" \
    -v p2_thresh="${POP2_THRESHOLD}" \
    -v samples_str="${ALL_SAMPLES_STR}" \
'
BEGIN {
    # 设置输出分隔符
    OFS="\t";
    # 将传入的样本名字符串分割成一个数组，用于后续匹配
    split(samples_str, sample_names, " ");
}
{
    # 为每一行（每个 Indel）重置计数器
    p1_count = 0;
    p2_count = 0;
    
    # 创建一个数组来暂存这一行转换后的 0/1 值
    delete row_values;

    # 从第5列（第一个样本）开始循环到最后一列
    for (i=5; i<=NF; i++) {
        # 获取当前列对应的样本名 (数组索引是 i-4, 因为awk数组从1开始)
        current_sample = sample_names[i-4];
        
        # 判断基因型是否存在变异 (0/1 或 1/1)
        is_present = 0;
        if ($i == "0/1" || $i == "1/1" || $i == "0|1" || $i == "1|0" || $i == "1|1") {
            is_present = 1;
        }
        
        # 存储转换后的 0 或 1
        row_values[i-4] = is_present;

        # 如果存在变异，则根据样本名增加对应群体的计数
        if (is_present == 1) {
            if (index(current_sample, p1_prefix) == 1) {
                p1_count++;
            } else if (index(current_sample, p2_prefix) == 1) {
                p2_count++;
            }
        }
    }

    # 检查是否满足任一群体的阈值条件
    if (p1_count >= p1_thresh || p2_count >= p2_thresh) {
        # 如果满足，则打印这一行的所有信息
        printf "%s\t%s\t%s\t%s", $1, $2, $3, $4;
        # 打印所有样本的 0/1 值
        for (j=1; j<=length(row_values); j++) {
            printf "\t%d", row_values[j];
        }
        printf "\n";
    }
}
' >> ${OUTPUT_MATRIX} # 将满足条件的行追加到输出文件中

# --- 步骤 4: 清理 ---
echo "步骤 4: 清理临时文件..."
rm -f samples.all.tmp samples.pop1.tmp samples.pop2.tmp

echo "处理完成！"
echo "经过频率过滤的 0/1 矩阵已保存到: ${OUTPUT_MATRIX}"