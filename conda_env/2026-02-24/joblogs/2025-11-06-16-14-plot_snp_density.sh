#!/bin/bash
set -e # 如果任何命令失败，则立即退出脚本

# --- 1. 用户配置区 ---
# 请根据您的实际文件名修改下面的变量
VCF_FILE="variation.filtered.snp.vcf.gz"
GENOME_FASTA="genome.fasta" # 您的参考基因组文件
WINDOW_SIZE=1000000 # 计算SNP密度的窗口大小，单位: bp (例如 1000000 = 1Mb)
OUTPUT_PREFIX="snp_density" # 输出图片的前缀

# --- 脚本区，通常无需修改 ---

echo "=== 步骤 1: 创建染色体名称映射文件 ==="
cat << EOF > chr.map
NC_081805.1  Chr1
NC_081806.1  Chr2
NC_081807.1  Chr3
NC_081808.1  Chr4
NC_081809.1  Chr5
NC_081810.1  Chr6
NC_081811.1  Chr7
NC_081812.1  Chr8
NC_081813.1  Chr9
NC_081814.1 Chr10
NC_081815.1 Chr11
NC_081816.1 Chr12
NC_081817.1 Chr13
NC_081818.1 Chr14
NC_081819.1 Chr15
NC_081820.1 Chr16
NC_081821.1 Chr17
NC_081822.1 Chr18
NC_081823.1 Chr19
EOF
echo "文件 'chr.map' 创建成功."

echo -e "\n=== 步骤 2: 生成染色体长度文件 ==="
# 检查fasta索引是否存在，如果不存在则创建
if [ ! -f "${GENOME_FASTA}.fai" ]; then
    echo "正在为 ${GENOME_FASTA} 创建索引..."
    samtools faidx ${GENOME_FASTA}
fi
# 使用awk根据映射文件转换染色体名称并提取长度
awk 'FNR==NR{map[$1]=$2; next} {if ($1 in map) print map[$1] "\t" $2}' chr.map ${GENOME_FASTA}.fai > genome.len
echo "文件 'genome.len' 创建成功."

echo -e "\n=== 步骤 3: 从 VCF 文件中提取 SNP 位置并转换为 BED 格式 ==="
# 解压VCF -> 忽略头信息 -> 提取染色体和位置 -> 转换为BED格式 -> 替换染色体名
gunzip -c ${VCF_FILE} | grep -v '^#' | awk -v OFS='\t' '{print $1, $2-1, $2}' | \
awk 'FNR==NR{map[$1]=$2; next} {if ($1 in map) {$1=map[$1]; print}}' chr.map - > ${OUTPUT_PREFIX}.positions.bed
echo "文件 '${OUTPUT_PREFIX}.positions.bed' 创建成功."

echo -e "\n=== 步骤 4: 计算指定窗口大小的 SNP 密度 ==="
# 使用 JCVI 的工具按窗口统计 SNP 数量
python -m jcvi.formats.bed --window ${OUTPUT_PREFIX}.positions.bed ${WINDOW_SIZE} > ${OUTPUT_PREFIX}.density.bed
echo "文件 '${OUTPUT_PREFIX}.density.bed' 创建成功."

echo -e "\n=== 步骤 5: 创建 JCVI 绘图布局文件 ==="
cat << EOF > ${OUTPUT_PREFIX}.layout
# y, x, rotation, height, width
.canvas
8, 5, 0, 1200, 600

# seqid, length
.seqids
$(cat genome.len)

# trackname, color, min, max, height, file
.tracks
heatmap, viridis, 0, db, 80, ${OUTPUT_PREFIX}.density.bed

# options
.config
proportional=no
spacing=20
EOF
echo "文件 '${OUTPUT_PREFIX}.layout' 创建成功."

echo -e "\n=== 步骤 6: 使用 JCVI 绘制染色体密度图 ==="
python -m jcvi.graphics.karyotype ${OUTPUT_PREFIX}.layout
echo -e "\n🎉 绘图完成! 输出文件为 'karyotype.pdf' 和 'karyotype.png'."
