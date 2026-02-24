#!/bin-bash

# --- 用户配置区 ---

# 包含原始BAM文件的输入文件夹
BAM_DIR="/share/org/YZWL/yzwl_lixg/project/16.荠菜/02.each/bam"

# 筛选后BAM文件的输出文件夹
OUTPUT_DIR="/share/org/YZWL/yzwl_lixg/project/16.荠菜/19.筛选唯一比对上的reads"

# MAPQ 质量值阈值。只有 MAPQ >= 这个值的 reads 会被保留。
# 20 是一个常用的、比较严格的阈值，可以有效过滤掉多重比对。
MAPQ_THRESHOLD=20

# 使用的线程数，可以根据服务器性能调整
THREADS=8

# --- 脚本主体 ---

echo "--- 开始筛选唯一比对的 reads (使用 MAPQ >= $MAPQ_THRESHOLD) ---"
echo "输入文件夹: $BAM_DIR"
echo "输出文件夹: $OUTPUT_DIR"
echo "---------------------------------------------------------"

# 检查 samtools 是否已安装
if ! command -v samtools &> /dev/null
then
    echo "错误: samtools 未安装或不在您的 PATH 环境变量中。"
    echo "请先安装 samtools (推荐使用 conda: conda install -c bioconda samtools)"
    exit 1
fi

# 创建输出文件夹，如果它不存在的话
mkdir -p "$OUTPUT_DIR"

# 遍历输入文件夹中所有以 .bam 结尾的文件
for input_bam in "$BAM_DIR"/*.bam
do
    # 检查是否有匹配的文件
    if [ -f "$input_bam" ]; then
        # 从完整路径中提取文件名
        base_name=$(basename "$input_bam")
        
        # 构建输出文件的完整路径
        output_bam="$OUTPUT_DIR/$base_name"
        
        echo "正在处理: $base_name"
        
        # 使用 samtools view 进行筛选
        # -h: 输出结果中包含BAM文件的头信息 (header)，非常重要！
        # -@ $THREADS: 使用指定数量的线程
        # -q $MAPQ_THRESHOLD: 这是核心筛选条件，只保留 MAPQ >= 阈值的 reads
        # -o $output_bam: 指定输出文件名
        # -b: 指定输出格式为 BAM
        samtools view -h -@ "$THREADS" -q "$MAPQ_THRESHOLD" -o "$output_bam" -b "$input_bam"
        
        # 为新生成的BAM文件创建索引
        echo "正在为 $output_bam 创建索引..."
        samtools index -@ "$THREADS" "$output_bam"
        
        echo "完成处理: $base_name -> $output_bam"
        echo "---------------------------------------------------------"
    fi
done

echo "--- 所有文件处理完毕！ ---"
