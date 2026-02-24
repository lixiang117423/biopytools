#!/bin/bash

# ================= 配置区域 =================
# 输入数据目录
INPUT_DIR="/share/org/YZWL/yzwl_lixg/project/19.大豆疫霉菌/01.data/clean"
# 输出主目录
OUTPUT_DIR="/share/org/YZWL/yzwl_lixg/project/19.大豆疫霉菌/17.基因组大小估计"
# 软件路径
JELLYFISH="/share/org/YZWL/yzwl_lixg/.local/bin/jellyfish"
GENOMESCOPE_R="/share/org/YZWL/yzwl_lixg/software/scripts/genomescope.R"

# 参数设置
THREADS=88
KMERS=(21 25 31)
# Jellyfish hash大小预估，大豆疫霉菌基因组约95M，设1G足够
HASH_SIZE="1G" 

# 最终汇总表文件
SUMMARY_FILE="${OUTPUT_DIR}/All_Samples_GenomeScope_Summary.csv"

# ===========================================

# 创建输出目录
mkdir -p "$OUTPUT_DIR"

# 初始化汇总表头
echo "Sample,Kmer,Genome_Size_Min,Genome_Size_Max,Heterozygosity_Min,Heterozygosity_Max" > "$SUMMARY_FILE"

echo "开始分析..."
echo "输入目录: $INPUT_DIR"
echo "输出目录: $OUTPUT_DIR"
echo "汇总文件: $SUMMARY_FILE"

# 遍历样品 (寻找 _1.clean.fq.gz 结尾的文件)
for r1_file in ${INPUT_DIR}/*_1.clean.fq.gz; do
    # 检查文件是否存在
    [ -e "$r1_file" ] || continue

    # 获取文件名和样品名
    filename=$(basename "$r1_file")
    # 去除后缀得到样品名 (例如 T30-4)
    sample_name=${filename/_1.clean.fq.gz/}
    
    # 定义 R2 文件路径
    r2_file="${INPUT_DIR}/${sample_name}_2.clean.fq.gz"

    if [ ! -f "$r2_file" ]; then
        echo "警告: 未找到样品 $sample_name 的 R2 文件 ($r2_file)，跳过。"
        continue
    fi

    echo "=========================================="
    echo "正在处理样品: $sample_name"
    echo "R1: $r1_file"
    echo "R2: $r2_file"

    # 对每个 Kmer 进行循环
    for k in "${KMERS[@]}"; do
        echo "  正在分析 K-mer = $k ..."
        
        # 定义该样品该Kmer的输出子目录
        current_out_dir="${OUTPUT_DIR}/${sample_name}_k${k}"
        mkdir -p "$current_out_dir"
        
        jf_file="${current_out_dir}/${sample_name}_k${k}.jf"
        histo_file="${current_out_dir}/${sample_name}_k${k}.histo"

        # 1. 运行 Jellyfish Count
        # 使用 <(zcat ...) 避免解压文件，节省空间
        echo "    1. 运行 Jellyfish Count..."
        $JELLYFISH count -C -m $k -s $HASH_SIZE -t $THREADS -o "$jf_file" \
            <(zcat "$r1_file") <(zcat "$r2_file")

        # 2. 运行 Jellyfish Histo
        echo "    2. 运行 Jellyfish Histo..."
        $JELLYFISH histo -t $THREADS "$jf_file" > "$histo_file"

        # 3. 删除巨大的 .jf 文件以节省空间 (建议保留)
        rm "$jf_file"

        # 4. 运行 GenomeScope2
        echo "    3. 运行 GenomeScope R 脚本..."
        # 注意：这里假设是 GenomeScope 2.0 的参数格式 (-i -o -k)
        # 如果报错，请尝试旧版格式：Rscript $GENOMESCOPE_R "$histo_file" $k 150 "$current_out_dir"
        Rscript "$GENOMESCOPE_R" -i "$histo_file" -o "$current_out_dir" -k $k -p 2 > /dev/null 2>&1

        # 5. 提取结果并汇总
        summary_txt="${current_out_dir}/summary.txt"
        
        if [ -f "$summary_txt" ]; then
            # 提取基因组大小 (Genome Haploid Length) 的范围
            # 格式通常为: Genome Haploid Length   min_bp   max_bp
            gsize_min=$(grep "Genome Haploid Length" "$summary_txt" | awk '{print $4}' | sed 's/,//g')
            gsize_max=$(grep "Genome Haploid Length" "$summary_txt" | awk '{print $5}' | sed 's/,//g')
            
            # 提取杂合度 (Heterozygosity)
            # 格式通常为: Heterozygosity   min%   max%
            het_min=$(grep "Heterozygosity" "$summary_txt" | awk '{print $3}')
            het_max=$(grep "Heterozygosity" "$summary_txt" | awk '{print $4}')

            echo "${sample_name},${k},${gsize_min},${gsize_max},${het_min},${het_max}" >> "$SUMMARY_FILE"
            echo "    -> 结果已写入汇总表。"
        else
            echo "    错误: 未生成 summary.txt，GenomeScope 可能运行失败。"
            echo "${sample_name},${k},NA,NA,NA,NA" >> "$SUMMARY_FILE"
        fi

    done
done

echo "=========================================="
echo "所有任务完成！"
echo "汇总结果已保存至: $SUMMARY_FILE"
