#!/bin/bash

# --- 使用说明 ---
# 将此脚本与你的两个基因组文件放在同一个目录下
# 然后在终端运行: bash find_large_indels.sh
# ----------------

# --- 用户需配置的变量 ---
REF_GENOME="120.sorted.genome.fa"  # 参考基因组 (Reference)
QUERY_GENOME="119.sorted.genome.fa" # 查询基因组 (Query)
PREFIX="120_vs_119"               # 输出文件的前缀
MIN_INDEL_LENGTH=15               # 你想筛选的INDEL最小长度

# --- 脚本主体 ---

echo "步骤 1: 使用 nucmer 进行全基因组比对..."
# --maxmatch: 寻找最长的唯一匹配序列
# -c 100: 设定一个合理的最小匹配簇长度
# -p ${PREFIX}: 指定输出文件的前缀
nucmer --maxmatch -t 80 -c 100 -p ${PREFIX} ${REF_GENOME} ${QUERY_GENOME}

# 检查 nucmer 是否成功运行
if [ ! -f "${PREFIX}.delta" ]; then
    echo "错误: Nucmer 比对失败，未生成 .delta 文件。"
    exit 1
fi

echo "步骤 2: 使用 show-snps 提取 INDELs..."
# -C: 只报告INDEL，不报告INDEL内部的SNP
# -T: 输出为制表符分隔的表格，方便后续处理
# -l: 在输出中包含长度信息
show-snps -C -T -l ${PREFIX}.delta > ${PREFIX}.snps

# 检查 show-snps 是否成功运行
if [ ! -s "${PREFIX}.snps" ]; then
    echo "错误: show-snps 未能提取差异位点。"
    exit 1
fi

echo "步骤 3: 使用 awk 筛选并格式化大于 ${MIN_INDEL_LENGTH}bp 的 INDELs..."
# awk 脚本解释:
# BEGIN{...}: 在处理文件前，先打印表头。OFS="\t"表示输出字段以tab分隔。
# NR > 4: show-snps输出的前4行是头信息，跳过它们。
# $2 == ".": 表示在参考基因组的这个位置是个空位，对应查询基因组的插入(Insertion)。
# $3 == ".": 表示在查询基因组的这个位置是个空位，对应查询基因组的缺失(Deletion)。
# len_ref=$7; len_qry=$8: 获取参考和查询序列中INDEL的长度。
# if (len_qry > MIN_INDEL_LENGTH): 如果是插入且长度大于阈值，则输出。
# if (len_ref > MIN_INDEL_LENGTH): 如果是缺失且长度大于阈值，则输出。
# 输出格式: 参考染色体, 起始位置, 终止位置, 类型, 长度, 查询染色体
awk -v min_len="${MIN_INDEL_LENGTH}" '
BEGIN {
    OFS="\t";
    print "Ref_Chr\tStart\tEnd\tType\tLength\tQuery_Chr";
}
NR > 4 {
    ref_chr = $10;
    query_chr = $11;
    ref_pos = $1;
    len_ref = $7;
    len_qry = $8;

    # 类型为插入 (相对于参考基因组120)
    if ($2 == "." && len_qry > min_len) {
        # 插入事件在参考基因组上只是一个点，我们用 pos 和 pos+1 表示这个间隙
        print ref_chr, ref_pos, ref_pos + 1, "INS", len_qry, query_chr;
    }
    # 类型为缺失 (相对于参考基因组120)
    else if ($3 == "." && len_ref > min_len) {
        # 缺失事件在参考基因组上是一个区间
        print ref_chr, ref_pos, ref_pos + len_ref, "DEL", len_ref, query_chr;
    }
}' ${PREFIX}.snps > large_indels_gt${MIN_INDEL_LENGTH}bp.tsv

echo "完成！"
echo "结果已保存在文件: large_indels_gt${MIN_INDEL_LENGTH}bp.tsv"
echo "中间文件包括: ${PREFIX}.delta, ${PREFIX}.snps"
