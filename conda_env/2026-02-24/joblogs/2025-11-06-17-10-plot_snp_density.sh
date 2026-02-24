# 旧的、错误的命令
python -m jcvi.formats.bed bins snp_density.positions.bed 1000000 > ...
```根据这个报错，`jcvi` 的 `bins` 动作把第二个参数 (`1000000`) **误认为是一个文件名**，而不是我们想传递的“窗口大小”这个数值。

正确的做法是，我们需要明确地告诉 `bins` 动作，“1000000”是一个窗口大小的值。这通常需要使用一个**参数标志**，比如 `--window`。

### 正确的解决方案

我们必须把步骤 4 的命令格式从：
`python -m jcvi.formats.bed [动作] [输入文件] [窗口大小]`

修改为正确的格式：
`python -m jcvi.formats.bed [动作] [输入文件] --window [窗口大小]`

这个小小的改动，就能让 JCVI 正确理解我们的意图。

---

### 最终的、已修正的完整脚本

之前的反复报错都是因为旧的脚本内容没有被彻底更新。为了确保万无一失，请您**再一次**，也是最后一次，用下面的**完整内容**替换掉您的 `plot_snp_density.sh` 文件。

我已经将步骤 4 修正为正确的命令格式。

```bash
#!/bin/bash
set -e # 如果任何命令失败，则立即退出脚本

# --- 1. 用户配置区 ---
VCF_FILE="variation.filtered.snp.vcf.gz"
WINDOW_SIZE=1000000 
OUTPUT_PREFIX="snp_density"

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

echo -e "\n=== 步骤 2: 生成染色体长度文件 (来自VCF头文件) ==="
gunzip -c ${VCF_FILE} | grep '^##contig' | sed -e 's/.*<ID=//' -e 's/,length=/\t/' -e 's/>.*//' | \
awk 'FNR==NR{map[$1]=$2; next} {if ($1 in map) print map[$1] "\t" $2}' chr.map - > genome.len
echo "文件 'genome.len' 创建成功."

echo -e "\n=== 步骤 3: 从 VCF 文件中提取 SNP 位置并转换为 BED 格式 ==="
gunzip -c ${VCF_FILE} | grep -v '^#' | awk -v OFS='\t' '{print $1, $2-1, $2}' | \
awk 'FNR==NR{map[$1]=$2; next} {if ($1 in map) {$1=map[$1]; print}}' chr.map - > ${OUTPUT_PREFIX}.positions.bed
echo "文件 '${OUTPUT_PREFIX}.positions.bed' 创建成功."

echo -e "\n=== 步骤 4: 计算指定窗口大小的 SNP 密度 (最终修正版) ==="
# 正确的命令格式是使用 --window 标志来指定窗口大小
python -m jcvi.formats.bed bins ${OUTPUT_PREFIX}.positions.bed --window ${WINDOW_SIZE} > ${OUTPUT_PREFIX}.density.bed
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