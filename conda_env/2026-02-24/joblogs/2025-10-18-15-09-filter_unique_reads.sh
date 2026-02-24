#!/bin/bash

# --- 用户配置区 ---

# 包含原始BAM文件的输入文件夹
BAM_DIR="/share/org/YZWL/yzwl_lixg/project/16.荠菜/02.each/bam"

# 筛选后BAM文件的输出文件夹
OUTPUT_DIR="/share/org/YZWL/yzwl_lixg/project/16.荠菜/19.筛选唯一比对上的reads"

# 使用的线程数，可以根据服务器性能调整
THREADS=8

# --- 脚本主体 ---

echo "--- 开始筛选唯一比对的 reads ---"
echo "输入文件夹: $BAM_DIR"
echo "输出文件夹: $OUTPUT_DIR"
echo "-------------------------------------"

# 检查 samtools 是否已安装
if ! command -v samtools &> /dev/null
then
    echo "错误: samtools 未安装或不在您的 PATH 环境变量中。"
    echo "请先安装 samtools (推荐使用 conda: conda install -c bioconda samtools)"
    exit 1
fi

# 创建输出文件夹，如果它不存在的话
# -p 选项可以确保在父目录不存在时也一并创建，且目录已存在时不会报错
mkdir -p "$OUTPUT_DIR"

# 遍历输入文件夹中所有以 .bam 结尾的文件
for input_bam in "$BAM_DIR"/*.bam
do
    # 检查是否有匹配的文件
    if [ -f "$input_bam" ]; then
        # 从完整路径中提取文件名 (例如: 150.sorted.bam)
        base_name=$(basename "$input_bam")
        
        # 构建输出文件的完整路径
        output_bam="$OUTPUT_DIR/$base_name"
        
        echo "正在处理: $base_name"
        
        # 使用 samtools 进行筛选
        # samtools view: 查看和转换BAM/SAM文件
        # -h: 输出结果中包含BAM文件的头信息 (header)，这非常重要！
        # -@ $THREADS: 使用指定数量的线程
        # -d 'NH:1': 这是核心筛选条件，只选择 NH 标签值为 1 的 reads
        # $input_bam: 输入文件
        # -o $output_bam: 指定输出文件名
        # -b: 指定输出格式为 BAM
        samtools view -h -@ "$THREADS" -d 'NH:1' -o "$output_bam" -b "$input_bam"
        
        # 为新生成的BAM文件创建索引，方便后续工具（如IGV）使用
        echo "正在为 $output_bam 创建索引..."
        samtools index -@ "$THREADS" "$output_bam"
        
        echo "完成处理: $base_name -> $output_bam"
        echo "-------------------------------------"
    fi
done

echo "--- 所有文件处理完毕！ ---"
