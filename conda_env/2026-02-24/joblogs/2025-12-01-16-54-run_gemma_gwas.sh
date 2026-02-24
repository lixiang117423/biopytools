#!/bin/bash

# 🧬 GEMMA GWAS批量分析脚本
# 作者: 自动生成
# 日期: $(date +%Y-%m-%d)

set -e

# 📁 定义文件路径
VCF="/share/org/YZWL/yzwl_lixg/project/94.rice_gas/06.gemma_gwas/variation.filtered.snp.maf5.vcf.gz"
PHENO="/share/org/YZWL/yzwl_lixg/project/94.rice_gas/06.gemma_gwas/16.微生物和甲烷菌用于GWAS的数据.txt"
GEMMA="/share/org/YZWL/yzwl_lixg/.local/bin/gemma"

# 📊 定义参数
N_PCA=10  # PCA主成分数量
OUTDIR="gemma_results"

# ✅ 检查文件是否存在
echo "🔍 检查输入文件..."
if [ ! -f "$VCF" ]; then
    echo "❌ 错误: VCF文件不存在: $VCF"
    exit 1
fi

if [ ! -f "$PHENO" ]; then
    echo "❌ 错误: 表型文件不存在: $PHENO"
    exit 1
fi

if [ ! -f "$GEMMA" ]; then
    echo "❌ 错误: GEMMA程序不存在: $GEMMA"
    exit 1
fi

# 📂 创建输出目录
mkdir -p "$OUTDIR"
cd "$OUTDIR"

echo "📁 工作目录: $(pwd)"

# 🔄 步骤1: 将VCF转换为PLINK格式
echo ""
echo "🔄 步骤1: 转换VCF为PLINK格式..."
plink --vcf "$VCF" \
      --make-bed \
      --out genotype \
      --allow-extra-chr \
      --double-id

# 🔧 修复FAM文件的表型列（替换-9为第一个表型）
echo "   🔧 修复FAM文件的表型列..."
# 提取第一个表型列（第2列，因为第1列是样本ID）
tail -n +2 "$PHENO" | awk '{print $2}' > first_pheno.tmp

# 替换FAM文件的第6列
awk 'NR==FNR{pheno[NR]=$1; next} {$6=pheno[FNR]; print}' first_pheno.tmp genotype.fam > genotype.fam.tmp
mv genotype.fam.tmp genotype.fam
rm first_pheno.tmp
echo "   ✅ FAM文件已更新"

# 📝 步骤2: 准备表型文件(GEMMA格式)
echo ""
echo "📝 步骤2: 准备表型文件..."

# 获取表型列数和列名
N_COLS=$(head -1 "$PHENO" | awk '{print NF}')
echo "   表型文件共有 $N_COLS 列（包括样本ID列）"
echo "   表型数量: $((N_COLS - 1))"

# 获取表型列名
PHENO_NAMES=($(head -1 "$PHENO" | tr '\t' ' '))
echo "   列名: ${PHENO_NAMES[@]}"

# 检查样本数
echo "   表型文件中的样本数: $(($(wc -l < "$PHENO") - 1))"
echo "   基因型文件中的样本数: $(wc -l < genotype.fam)"

# 创建不含表头的表型文件供GEMMA使用
tail -n +2 "$PHENO" > pheno_no_header.txt
echo "   ✅ 已创建无表头表型文件"

# 🧮 步骤3: 计算PCA作为协变量
echo ""
echo "🧮 步骤3: 计算PCA (前${N_PCA}个主成分)..."
plink --bfile genotype \
      --pca ${N_PCA} \
      --out pca \
      --allow-extra-chr

echo "   📊 检查PCA结果的样本数..."
pca_samples=$(wc -l < pca.eigenvec)
fam_samples=$(wc -l < genotype.fam)
echo "   PCA结果样本数: $pca_samples"
echo "   基因型样本数: $fam_samples"

# 如果样本数不一致，需要更新基因型文件
if [ "$pca_samples" -ne "$fam_samples" ]; then
    echo "   ⚠️  样本数不一致，提取PCA中存在的样本..."
    awk '{print $1, $2}' pca.eigenvec > keep_samples.txt
    
    # 重新生成只包含PCA样本的基因型文件
    plink --bfile genotype \
          --keep keep_samples.txt \
          --make-bed \
          --out genotype_matched \
          --allow-extra-chr
    
    # 更新基因型文件路径
    mv genotype_matched.bed genotype.bed
    mv genotype_matched.bim genotype.bim
    mv genotype_matched.fam genotype.fam
    rm genotype_matched.log
    
    echo "   ✅ 基因型文件已更新为匹配的样本"
    
    # 同时更新表型文件，只保留这些样本
    awk 'NR==FNR{keep[$1]; next} ($1 in keep)' keep_samples.txt pheno_no_header.txt > pheno_matched.txt
    mv pheno_matched.txt pheno_no_header.txt
    echo "   ✅ 表型文件已更新为匹配的样本"
fi

# 准备协变量文件（GEMMA格式: 不需要FID/IID列，只需要协变量值）
tail -n +1 pca.eigenvec | awk '{for(i=3; i<=NF; i++) printf "%s%s", $i, (i==NF?"\n":"\t")}' > covariate.txt
echo "   协变量文件已准备"
echo "   协变量文件样本数: $(wc -l < covariate.txt)"

# 🔢 步骤4: 计算亲缘关系矩阵 (只需计算一次)
echo ""
echo "🔢 步骤4: 计算亲缘关系矩阵..."

# 最终验证所有文件的样本数一致性
echo "   🔍 最终样本数检查:"
geno_n=$(wc -l < genotype.fam)
pheno_n=$(wc -l < pheno_no_header.txt)
covar_n=$(wc -l < covariate.txt)
echo "   基因型: $geno_n 个样本"
echo "   表型: $pheno_n 个样本"
echo "   协变量: $covar_n 个样本"

if [ "$geno_n" -ne "$pheno_n" ] || [ "$geno_n" -ne "$covar_n" ]; then
    echo ""
    echo "   ❌ 错误: 样本数不一致！"
    echo "   请检查以下文件的样本ID是否匹配:"
    echo "   - genotype.fam"
    echo "   - pheno_no_header.txt"
    echo "   - covariate.txt"
    exit 1
fi

echo "   ✅ 所有文件样本数一致"

$GEMMA -bfile genotype \
        -gk 1 \
        -o kinship

if [ ! -f "output/kinship.cXX.txt" ]; then
    echo "   ❌ 亲缘关系矩阵计算失败"
    exit 1
fi

echo "   ✅ 亲缘关系矩阵计算完成"

# 🚀 步骤5: 遍历每个表型进行GWAS分析
echo ""
echo "🚀 步骤5: 开始GWAS分析..."
echo "==========================================="

for ((i=2; i<=$N_COLS; i++)); do
    PHENO_NAME="${PHENO_NAMES[$((i-1))]}"
    PHENO_COL=$((i-1))  # 列索引（去掉表头后的列号）
    
    echo ""
    echo "📊 分析表型 $((i-1))/$((N_COLS-1)): ${PHENO_NAME}"
    
    # 运行GEMMA LMM分析
    echo "   🧬 运行线性混合模型 (LMM)..."
    $GEMMA -bfile genotype \
            -k output/kinship.cXX.txt \
            -lmm 4 \
            -p pheno_no_header.txt \
            -n ${PHENO_COL} \
            -c covariate.txt \
            -o ${PHENO_NAME}_lmm \
            -notsnp
    
    echo "   ✅ ${PHENO_NAME} 分析完成"
    echo "   📄 结果文件: output/${PHENO_NAME}_lmm.assoc.txt"
done

echo ""
echo "==========================================="
echo "🎉 所有GWAS分析完成!"
echo "📂 结果目录: $(pwd)/output"
echo ""
echo "📋 结果文件列表:"
ls -lh output/*.assoc.txt 2>/dev/null || echo "   无结果文件"

# 📊 生成汇总统计
echo ""
echo "📊 生成分析汇总..."
echo "表型名称,SNP数量,显著SNP数(P<1e-5),显著SNP数(P<1e-6)" > gwas_summary.txt
for ((i=2; i<=$N_COLS; i++)); do
    PHENO_NAME="${PHENO_NAMES[$((i-1))]}"
    if [ -f "output/${PHENO_NAME}_lmm.assoc.txt" ]; then
        TOTAL_SNPS=$(tail -n +2 "output/${PHENO_NAME}_lmm.assoc.txt" | wc -l)
        SIG_SNPS_5=$(tail -n +2 "output/${PHENO_NAME}_lmm.assoc.txt" | awk '$13 < 1e-5' | wc -l)
        SIG_SNPS_6=$(tail -n +2 "output/${PHENO_NAME}_lmm.assoc.txt" | awk '$13 < 1e-6' | wc -l)
        echo "${PHENO_NAME},${TOTAL_SNPS},${SIG_SNPS_5},${SIG_SNPS_6}" >> gwas_summary.txt
    fi
done

echo ""
echo "📄 分析汇总已保存到: gwas_summary.txt"
cat gwas_summary.txt

echo ""
echo "✨ 全部完成! ✨"