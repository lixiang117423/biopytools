#!/bin/bash

# 设置工作目录
WORK_DIR="/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/45.候选区间的共线性/01.genome"
GENOME_LIST="$WORK_DIR/genome.id.txt"

cd "$WORK_DIR" || exit 1

echo "========================================"
echo "开始处理基因组文件"
echo "工作目录: $WORK_DIR"
echo "========================================"

# 1. 拷贝基因组文件
echo ""
echo "步骤1: 拷贝基因组文件..."
count=0
while read genome; do
  # 跳过空行和注释行
  [[ -z "$genome" ]] && continue
  [[ "$genome" == \#* ]] && continue
  
  if [ -f ~/database/brassicaceae/fasta/${genome}.fa ]; then
    cp ~/database/brassicaceae/fasta/${genome}.fa ./
    echo "  拷贝: ${genome}.fa"
    ((count++))
  else
    echo "  警告: 未找到 ${genome}.fa"
  fi
done < "$GENOME_LIST"

echo "  共拷贝 $count 个基因组文件"

# 2. 重命名序列ID
echo ""
echo "步骤2: 添加样品名前缀到序列ID..."
for file in *.fa; do
  sample=$(basename "$file" .fa)
  
  # 检查第一条序列ID
  first_id=$(grep "^>" "$file" | head -1 | sed 's/^>//' | sed 's/ .*//')
  
  if [[ "$first_id" == "${sample}_"* ]]; then
    echo "  跳过 $file (已有前缀)"
  else
    sed -i "s/^>/>${sample}_/g" "$file"
    echo "  处理: $file"
  fi
done

echo ""
echo "========================================"
echo "处理完成！"
echo "最终文件数: $(ls -1 *.fa 2>/dev/null | wc -l)"
echo "========================================"

cat 01.genome/*fa > all_genome.fa
minimap2 -x asm5 region.fa all_genome.fa > 02.minimap2/all_genome.paf
