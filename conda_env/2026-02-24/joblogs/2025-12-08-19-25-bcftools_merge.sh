#!/bin/bash

# ================= 配置区域 =================
# 输入文件夹 (g.vcf.gz 所在位置)
IN_DIR="/share/org/YZWL/yzwl_lixg/project/21.野生大豆群体/03.glnexus_joint/vcf"

# 参考基因组
REF="/share/org/YZWL/yzwl_lixg/project/21.野生大豆群体/01.data/genome/genome.fa"

# 输出文件夹
OUT_DIR="/share/org/YZWL/yzwl_lixg/project/21.野生大豆群体/03.第一批720个样品/03.bcftools_joint"

# 输出文件名
OUT_FILE="merged_720samples.vcf.gz"

# 使用的线程数 (根据服务器配置调整)
THREADS=10
# ===========================================

# 1. 创建输出目录
if [ ! -d "$OUT_DIR" ]; then
    mkdir -p "$OUT_DIR"
    echo "Created output directory: $OUT_DIR"
fi

# 进入输出目录，方便管理临时文件
cd "$OUT_DIR"

# 2. 生成文件列表
echo "Searching for g.vcf.gz files..."
# 使用 find 获取绝对路径，避免路径错误
find "$IN_DIR" -name "*.g.vcf.gz" > input_file_list.txt

FILE_COUNT=$(wc -l < input_file_list.txt)
echo "Found $FILE_COUNT files to merge."

if [ "$FILE_COUNT" -eq 0 ]; then
    echo "Error: No .g.vcf.gz files found in $IN_DIR"
    exit 1
fi

# 3. 检查并建立索引 (bcftools merge 必须要求输入文件有索引)
echo "Checking indexes for input files..."
# 遍历列表中的文件，检查是否存在 .tbi 或 .csi，如果不存在则建立索引
# 注意：如果文件很多，这一步可能耗时，建议确认上游步骤已经做好了索引
while read file; do
    if [ ! -f "$file.tbi" ] && [ ! -f "$file.csi" ]; then
        echo "Indexing $file ..."
        bcftools index -t "$file" --threads 4
    fi
done < input_file_list.txt

# 4. 调整 ulimit (防止 'Too many open files' 错误)
# 720个文件通常在默认限制内，但为了保险起见尝试提高限制
ulimit -n 4096 2>/dev/null || echo "Warning: Could not increase ulimit, proceeding anyway."

# 5. 执行合并
echo "Starting bcftools merge..."

# --file-list: 读取文件列表
# --output-type z: 输出为压缩的 vcf.gz
# --threads: 多线程压缩
# --merge all: 即使缺失某些信息也强制合并（针对 gVCF 常用）
# -0: (可选) 假设未见位点为 Ref，gVCF合并时视情况开启，这里默认使用标准合并

bcftools merge \
    --file-list input_file_list.txt \
    --output-type z \
    --output "$OUT_DIR/$OUT_FILE" \
    --threads "$THREADS"

# 6. 对合并后的结果建立索引
if [ -f "$OUT_DIR/$OUT_FILE" ]; then
    echo "Merge finished. Indexing output file..."
    bcftools index -t "$OUT_DIR/$OUT_FILE" --threads "$THREADS"
    echo "All done! Output file is at:"
    echo "$OUT_DIR/$OUT_FILE"
else
    echo "Error: Output file was not generated."
    exit 1
fi
