#!/bin/bash

# ==============================================================================
# 配置区域 (已根据您的 .fai 文件更新)
# ==============================================================================

# 1. 输入与输出
HIC_FILE="/share/org/YZWL/yzwl_lixg/project/22.Est1_hic/03.hic/fithic2/inter_30.hic"
OUT_DIR="/share/org/YZWL/yzwl_lixg/project/22.Est1_hic/03.hic/fithic2/result_5kb"

# 2. 软件与环境路径
FITHIC_BIN="/share/org/YZWL/yzwl_lixg/miniforge3/envs/fithic-v.2.0.8/bin/fithic"
UTILS_DIR="/share/org/YZWL/yzwl_lixg/software/fithic2/utils"
# Juicer Tools Jar
JUICER_JAR=~/software/juicer/scripts/juicer_tools.3.0.0.jar
# Python 解释器
PYTHON_EXE="/share/org/YZWL/yzwl_lixg/miniforge3/envs/fithic-v.2.0.8/bin/python"

# 3. 运行参数
RESOLUTION=5000       # 5kb
LOWER_BOUND=2000      # 2kb
UPPER_BOUND=3000000   # 3Mb
QVALUE_THRESH=0.01

# 4. 染色体列表 (仅分析核基因组 Chr1-Chr5)
CHROMS=("Chr1" "Chr2" "Chr3" "Chr4" "Chr5")

# ==============================================================================
# 脚本逻辑执行区
# ==============================================================================

mkdir -p "$OUT_DIR"
cd "$OUT_DIR" || exit

echo ">>> [Start] FitHiC Pipeline (Custom Genome Sizes)"
echo ">>> Resolution: ${RESOLUTION} bp"

# --------------------------------------------------------
# Step 1: 生成 Chrom Sizes 文件 (严格匹配您的 genome.fa.fai)
# --------------------------------------------------------
echo ">>> Step 1: Creating chrom.sizes..."
CHROM_SIZES="chrom.sizes"

# 注意：这里使用的是您提供的精确长度
# 格式: ChrName <TAB> Length
cat > $CHROM_SIZES <<EOF
Chr1	29412698
Chr2	19420461
Chr3	23119851
Chr4	18582608
Chr5	26592585
EOF
# CP, MT, Chr0 已排除，因为它们不适合做标准的 cis-loop calling

# --------------------------------------------------------
# Step 2: 生成 Fragments 文件
# --------------------------------------------------------
echo ">>> Step 2: Creating Fragments file..."
FRAGS_FILE="fragments_${RESOLUTION}.txt"

if [ ! -f "${FRAGS_FILE}.gz" ]; then
    "$PYTHON_EXE" "$UTILS_DIR/createFitHiCFragments-fixedsize.py" \
        --chrLens "$CHROM_SIZES" \
        --resolution "$RESOLUTION" \
        --outFile "$FRAGS_FILE"
    
    gzip -f "$FRAGS_FILE"
    echo "    -> Generated: ${FRAGS_FILE}.gz"
else
    echo "    -> Fragments file already exists, skipping."
fi

# --------------------------------------------------------
# Step 3: 按染色体循环处理
# --------------------------------------------------------
HALF_RES=$(($RESOLUTION / 2))

for CHR in "${CHROMS[@]}"; do
    echo "--------------------------------------------------"
    echo ">>> Processing Chromosome: ${CHR}"
    
    CHR_DIR="${OUT_DIR}/${CHR}"
    mkdir -p "$CHR_DIR"
    
    # [3.1] 提取互作矩阵
    echo "    [3.1] Dumping Observed Counts from .hic..."
    INTER_RAW="${CHR_DIR}/raw_inter.txt"
    INTER_GZ="${CHR_DIR}/interactions.txt.gz"
    
    # 使用 -Xmx32g 保证内存充足
    java -Xmx32g -jar "$JUICER_JAR" dump observed NONE "$HIC_FILE" "$CHR" "$CHR" BP "$RESOLUTION" > "$INTER_RAW"
    
    # 检查 dump 是否成功
    if [ ! -s "$INTER_RAW" ]; then
        echo "    ERROR: Dump failed for ${CHR}. File is empty. Check chromosome name in .hic file!"
        rm "$INTER_RAW"
        continue
    fi
    
    # 转换坐标 (Start -> Midpoint) 并压缩
    awk -v offset="$HALF_RES" '{printf "%s\t%d\t%s\t%d\t%d\n", $1, $2+offset, $3, $4+offset, $5}' "$INTER_RAW" | gzip > "$INTER_GZ"
    rm "$INTER_RAW"
    
    # [3.2] 提取 KR Bias 向量
    echo "    [3.2] Dumping KR Norm Vector..."
    BIAS_RAW="${CHR_DIR}/raw_kr.txt"
    BIAS_GZ="${CHR_DIR}/bias.txt.gz"
    
    java -Xmx32g -jar "$JUICER_JAR" dump norm KR "$HIC_FILE" "$CHR" BP "$RESOLUTION" > "$BIAS_RAW"
    
    # 提取 Fragment 用于合并
    zcat "${FRAGS_FILE}.gz" | awk -v c="$CHR" '$1==c {print $0}' > "${CHR_DIR}/tmp_frags.txt"
    
    # 合并 Fragment 和 Bias (去除 NaN)
    # 增加额外的列数判断，防止空文件报错
    paste "${CHR_DIR}/tmp_frags.txt" "$BIAS_RAW" | \
    awk '{if(NF>=3 && $NF!="NaN" && $NF!="0" && $NF!="") print $1"\t"$3"\t"$NF}' | gzip > "$BIAS_GZ"
    
    rm "$BIAS_RAW" "${CHR_DIR}/tmp_frags.txt"
    
    # [3.3] 运行 FitHiC
    echo "    [3.3] Running FitHiC..."
    PREFIX="At_${CHR}"
    
    "$FITHIC_BIN" \
        -f "${FRAGS_FILE}.gz" \
        -i "$INTER_GZ" \
        -t "$BIAS_GZ" \
        -o "$CHR_DIR" \
        -r "$RESOLUTION" \
        -L "$LOWER_BOUND" \
        -U "$UPPER_BOUND" \
        -p 2 \
        -l "$PREFIX" \
        -v
        
    # [3.4] 合并与过滤结果
    echo "    [3.4] Merging and Filtering (Q < $QVALUE_THRESH)..."
    SIG_FILE="${CHR_DIR}/${PREFIX}.spline_pass2.significances.txt.gz"
    
    if [ -f "$SIG_FILE" ]; then
        # 显式传递 utils 路径
        bash "$UTILS_DIR/merge-filter.sh" \
            "$SIG_FILE" \
            "$RESOLUTION" \
            "${CHR_DIR}/merged_results" \
            "$QVALUE_THRESH" \
            "$UTILS_DIR"
            
        echo "    -> Done! Check results in: ${CHR_DIR}/merged_results"
    else
        echo "    !!! ERROR: FitHiC output not found: $SIG_FILE"
    fi

done

echo "=================================================="
echo ">>> All Finished."
echo "=================================================="
