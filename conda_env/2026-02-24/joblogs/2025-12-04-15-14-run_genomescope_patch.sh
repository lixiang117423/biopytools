#!/bin/bash

# ================= 配置区域 =================
# 输入数据目录
INPUT_DIR="/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/01.data/clean"
# 输出主目录
OUTPUT_DIR="/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/16.genomescope评估重测序样品基因组大小"
# 软件路径
JELLYFISH="/share/org/YZWL/yzwl_lixg/.local/bin/jellyfish"
GENOMESCOPE_R="/share/org/YZWL/yzwl_lixg/software/scripts/genomescope.R"

# 参数设置
THREADS=88
# 测序读长 (根据您之前的脚本默认为 150bp)
READ_LENGTH=150
# Kmer 梯度
KMERS=(21 25 31)
# Jellyfish hash大小 (大豆疫霉菌基因组较小，10G足够)
HASH_SIZE="10G" 
# GenomeScope 最大覆盖度限制 (对应旧版脚本的参数)
MAX_KMER_COV=1000

# 最终汇总表文件
SUMMARY_FILE="${OUTPUT_DIR}/All_Samples_GenomeScope_Summary.csv"

# ===========================================

# 1. 准备环境
mkdir -p "$OUTPUT_DIR"
# 初始化汇总表头 (CSV格式)
echo "Sample,Kmer,Genome_Size_Min,Genome_Size_Max,Heterozygosity_Min,Heterozygosity_Max,Model_Fit" > "$SUMMARY_FILE"

echo "=========================================="
echo "开始批量分析..."
echo "输入目录: $INPUT_DIR"
echo "输出目录: $OUTPUT_DIR"
echo "汇总文件: $SUMMARY_FILE"
echo "=========================================="

# 2. 遍历样品 (寻找 _1.clean.fq.gz 结尾的文件)
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
        echo "警告: 未找到样品 $sample_name 的 R2 文件，跳过。"
        continue
    fi

    echo ">>> 正在处理样品: $sample_name"

    # 3. 对每个 Kmer 进行循环
    for k in "${KMERS[@]}"; do
        
        # 定义该样品该Kmer的输出子目录
        current_out_dir="${OUTPUT_DIR}/${sample_name}_k${k}"
        mkdir -p "$current_out_dir"
        
        jf_file="${current_out_dir}/${sample_name}_k${k}.jf"
        histo_file="${current_out_dir}/${sample_name}_k${k}.histo"
        summary_txt="${current_out_dir}/summary.txt"

        # --- 步骤 1: Jellyfish Count ---
        # 如果 histo 文件不存在，才运行 jellyfish (实现断点续传功能)
        if [ ! -f "$histo_file" ]; then
            echo "    [K=$k] 1. 运行 Jellyfish Count & Histo..."
            
            # 使用流式传输 (zcat) 避免解压
            $JELLYFISH count -C -m $k -s $HASH_SIZE -t $THREADS -o "$jf_file" \
                <(zcat "$r1_file") <(zcat "$r2_file")
            
            # 导出 Histo
            $JELLYFISH histo -t $THREADS "$jf_file" > "$histo_file"
            
            # 删除巨大的 .jf 文件以节省空间
            rm "$jf_file"
        else
            echo "    [K=$k] 1. Histo 文件已存在，跳过 Jellyfish 步骤。"
        fi

        # --- 步骤 2: GenomeScope (旧版 V1.0 语法) ---
        echo "    [K=$k] 2. 运行 GenomeScope R..."
        
        # 旧版语法: Rscript 脚本 histo文件 kmer长度 读长 输出目录 最大覆盖度
        Rscript "$GENOMESCOPE_R" "$histo_file" $k $READ_LENGTH "$current_out_dir" $MAX_KMER_COV > /dev/null 2>&1

        # --- 步骤 3: 提取结果并写入汇总表 ---
        if [ -f "$summary_txt" ]; then
            # 提取数据 (自动处理逗号和百分号)
            # Genome Haploid Length (Min/Max) - 通常在第4和6列
            gsize_min=$(grep "Genome Haploid Length" "$summary_txt" | awk '{print $4}' | sed 's/,//g')
            gsize_max=$(grep "Genome Haploid Length" "$summary_txt" | awk '{print $6}' | sed 's/,//g')
            
            # Heterozygosity (Min/Max) - 通常在第2和3列 (例如 1.05% 1.08%)
            het_min=$(grep "Heterozygosity" "$summary_txt" | awk '{print $2}' | sed 's/%//g')
            het_max=$(grep "Heterozygosity" "$summary_txt" | awk '{print $3}' | sed 's/%//g')

            # Model Fit
            model_fit=$(grep "Model Fit" "$summary_txt" | awk '{print $3}' | sed 's/%//g')

            # 写入 CSV
            echo "${sample_name},${k},${gsize_min},${gsize_max},${het_min},${het_max},${model_fit}" >> "$SUMMARY_FILE"
        else
            echo "    [K=$k] 错误: GenomeScope 运行失败 (未生成 summary.txt)。"
            echo "${sample_name},${k},NA,NA,NA,NA,NA" >> "$SUMMARY_FILE"
        fi

    done
    echo "    --------------------------------"
done

echo "=========================================="
echo "所有任务完成！"
echo "汇总结果已保存至: $SUMMARY_FILE"