#!/bin/bash

# ==============================================================================
# 配置区域
# ==============================================================================

# 1. 输入与输出
HIC_FILE="/share/org/YZWL/yzwl_lixg/project/22.Est1_hic/03.hic/fithic2/inter_30.hic"
OUT_DIR="/share/org/YZWL/yzwl_lixg/project/22.Est1_hic/03.hic/fithic2/result_5kb"

# 2. 软件与环境路径
FITHIC_BIN="/share/org/YZWL/yzwl_lixg/miniforge3/envs/fithic-v.2.0.8/bin/fithic"
UTILS_DIR="/share/org/YZWL/yzwl_lixg/software/fithic2/utils"

# [已更新] 使用您验证通过的 1.9.9 版本 JAR 包
JUICER_JAR=~/software/juicer/scripts/juicer_tools.1.9.9_jcuda.0.8.jar

# Python 解释器
PYTHON_EXE="/share/org/YZWL/yzwl_lixg/miniforge3/envs/fithic-v.2.0.8/bin/python"

# 3. 运行参数 (针对拟南芥优化)
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

echo ">>> [Start] FitHiC Pipeline v6 (Jar 1.9.9 + Plain Text)"
echo ">>> Resolution: ${RESOLUTION} bp"
echo ">>> Using JAR: $JUICER_JAR"

# --------------------------------------------------------
# Step 1: 生成 Chrom Sizes 文件
# --------------------------------------------------------
echo ">>> Step 1: Creating chrom.sizes..."
CHROM_SIZES="chrom.sizes"

if [ ! -f "$CHROM_SIZES" ]; then
cat > $CHROM_SIZES <<EOF
Chr1	29412698
Chr2	19420461
Chr3	23119851
Chr4	18582608
Chr5	26592585
EOF
fi

# --------------------------------------------------------
# Step 2: 生成 Fragments 文件 (纯文本，不压缩)
# --------------------------------------------------------
echo ">>> Step 2: Creating Fragments file..."
FRAGS_FILE="fragments_${RESOLUTION}.txt"

# 关键修改：去掉 .gz 后缀，不压缩，防止 Python 报错
if [ ! -f "${FRAGS_FILE}" ]; then
    "$PYTHON_EXE" "$UTILS_DIR/createFitHiCFragments-fixedsize.py" \
        --chrLens "$CHROM_SIZES" \
        --resolution "$RESOLUTION" \
        --outFile "$FRAGS_FILE"
    
    echo "    -> Generated: ${FRAGS_FILE}"
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
    echo "    [3.1] Dumping Observed Counts..."
    INTER_RAW="${CHR_DIR}/raw_inter.txt"
    INTER_TXT="${CHR_DIR}/interactions.txt" # 使用 .txt
    
    # 尝试运行 dump 命令
    java -Xmx32g -jar "$JUICER_JAR" dump observed NONE "$HIC_FILE" "$CHR" "$CHR" BP "$RESOLUTION" > "$INTER_RAW"
    
    # 检查文件是否生成成功
    if [ ! -s "$INTER_RAW" ]; then
        echo "    !!! CRITICAL ERROR: Dump failed for ${CHR}. File is empty."
        rm "$INTER_RAW"
        continue
    fi
    
    # 转换坐标 (Start -> Midpoint)，输出为纯文本
    awk -v offset="$HALF_RES" '{printf "%s\t%d\t%s\t%d\t%d\n", $1, $2+offset, $3, $4+offset, $5}' "$INTER_RAW" > "$INTER_TXT"
    rm "$INTER_RAW"
    
    # [3.2] 提取 KR Bias 向量
    echo "    [3.2] Dumping KR Norm Vector..."
    BIAS_RAW="${CHR_DIR}/raw_kr.txt"
    BIAS_TXT="${CHR_DIR}/bias.txt" # 使用 .txt
    
    java -Xmx32g -jar "$JUICER_JAR" dump norm KR "$HIC_FILE" "$CHR" BP "$RESOLUTION" > "$BIAS_RAW"
    
    # 提取 Fragment 用于合并 (直接 cat，不 zcat)
    awk -v c="$CHR" '$1==c {print $0}' "$FRAGS_FILE" > "${CHR_DIR}/tmp_frags.txt"
    
    # 合并 Fragment 和 Bias (纯文本)
    paste "${CHR_DIR}/tmp_frags.txt" "$BIAS_RAW" | \
    awk '{if(NF>=3 && $NF!="NaN" && $NF!="0" && $NF!="") print $1"\t"$3"\t"$NF}' > "$BIAS_TXT"
    
    rm "$BIAS_RAW" "${CHR_DIR}/tmp_frags.txt"
    
    # [3.3] 运行 FitHiC
    echo "    [3.3] Running FitHiC..."
    PREFIX="At_${CHR}"
    
    # 传入纯文本文件
    "$FITHIC_BIN" \
        -f "${FRAGS_FILE}" \
        -i "$INTER_TXT" \
        -t "$BIAS_TXT" \
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