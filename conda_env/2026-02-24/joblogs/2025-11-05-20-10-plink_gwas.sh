#!/bin/bash

# PLINK GWAS分析脚本 - 支持非标准染色体命名
VCF="variation.filtered.snp.vcf.gz"
PHENO="phe.txt"
OUT="DI_gwas_dominant"

echo "==============================================="
echo "GWAS分析 - 显性模型（支持非标准染色体）"
echo "==============================================="

echo "Step 1: 转换VCF到PLINK格式..."
plink --vcf ${VCF} \
      --double-id \
      --allow-extra-chr \
      --set-missing-var-ids @:#:\$1:\$2 \
      --make-bed \
      --out ${OUT}

# 检查转换是否成功
if [ ! -f "${OUT}.bed" ]; then
    echo "错误：VCF转换失败！"
    exit 1
fi

echo "Step 2: 检查数据..."
echo "SNP数量: $(wc -l < ${OUT}.bim)"
echo "样本数量: $(wc -l < ${OUT}.fam)"
echo ""
echo "前5个样本ID（VCF）："
awk '{print $1,$2}' ${OUT}.fam | head -5
echo ""
echo "前5个样本ID（表型文件）："
head -6 ${PHENO}

echo "Step 3: 准备表型文件..."
# PLINK表型格式：FID IID phenotype
# 二分类：1=control, 2=case
awk 'NR==1{next} {
    if($2==0) pheno=1; 
    else if($2==1) pheno=2; 
    else pheno=-9;
    print $1,$1,pheno
}' ${PHENO} > ${OUT}.pheno

echo "表型文件前5行："
head -5 ${OUT}.pheno

echo "Step 4: 运行逻辑回归GWAS（显性模型）..."
plink --bfile ${OUT} \
      --pheno ${OUT}.pheno \
      --logistic dominant hide-covar \
      --ci 0.95 \
      --allow-extra-chr \
      --allow-no-sex \
      --out ${OUT}_logistic

# 检查是否生成了结果
if [ -f "${OUT}_logistic.assoc.logistic" ]; then
    echo "✓ GWAS分析成功！"
    
    echo "Step 5: 基础关联分析..."
    plink --bfile ${OUT} \
          --pheno ${OUT}.pheno \
          --assoc \
          --allow-extra-chr \
          --allow-no-sex \
          --out ${OUT}_assoc
    
    echo "Step 6: Fisher精确检验（显性模型）..."
    plink --bfile ${OUT} \
          --pheno ${OUT}.pheno \
          --model dom \
          --allow-extra-chr \
          --allow-no-sex \
          --out ${OUT}_model
    
    echo "Step 7: 多重检验校正..."
    if [ -f "${OUT}_assoc.assoc" ]; then
        plink --bfile ${OUT} \
              --pheno ${OUT}.pheno \
              --assoc \
              --adjust \
              --allow-extra-chr \
              --allow-no-sex \
              --out ${OUT}_adjust
    fi
    
    echo "Step 8: 提取显性模型结果..."
    grep "DOM" ${OUT}_logistic.assoc.logistic > ${OUT}_dominant_results.txt
    awk '{print $1,$3,$12}' ${OUT}_dominant_results.txt > ${OUT}_plot.txt
    
    # 统计结果
    echo ""
    echo "==============================================="
    echo "分析结果摘要："
    echo "==============================================="
    TOTAL=$(wc -l < ${OUT}_dominant_results.txt)
    SIG5=$(awk '$12 < 1e-5' ${OUT}_dominant_results.txt | wc -l)
    # SIG8=$(awk '$12 < 5e-8' ${OUT}_dominant_results.txt | wc -l)
    BONFERRONI_THRESHOLD=$(echo "0.05 / ${TOTAL}" | bc -l)
    echo "📈 计算得到的BONFERRONI校正阈值为: ${BONFERRONI_THRESHOLD}"
    #    使用 awk 将阈值作为变量传入，并进行筛选和统计
    SIG_BONFERRONI=$(awk -v threshold="${BONFERRONI_THRESHOLD}" '$12 < threshold' ${OUT}_dominant_results.txt | wc -l)
    
    echo "总SNP数: ${TOTAL}"
    echo "Suggestive (P<1e-5): ${SIG5}"
    echo "BONFERRONI significant (0.05/总SNP数): ${SIG_BONFERRONI}"
    echo ""
    echo "Top 10 关联位点："
    echo "CHR         SNP                    POS            P-value"
    echo "---------------------------------------------------------"
    sort -k9 -g ${OUT}_dominant_results.txt | head -10 | \
    awk '{printf "%-10s %-20s %-12s %8.3f\n", $1, $2, $3, $12}'
    
else
    echo "✗ GWAS分析失败！"
    echo ""
    echo "可能的原因和解决方案："
    echo "1. 样本ID不匹配"
    echo "   检查: comm -3 <(awk '{print $1}' ${OUT}.fam | sort) <(tail -n +2 ${PHENO} | cut -f1 | sort)"
    echo ""
    echo "2. 所有样本表型缺失"
    echo "   检查表型分布: awk '{print \$6}' ${OUT}.fam | sort | uniq -c"
    echo ""
    echo "3. 使用--1参数如果表型编码有问题"
fi

# 生成简单绘图脚本
cat > plot_results.R << 'EOF'
# 安装包（如果需要）
if(!require(ggplot2)) install.packages("ggplot2")

# 读取数据
data <- read.table("DI_gwas_dominant_plot.txt", header=F)
colnames(data) <- c("CHR", "BP", "P")

# -log10转换
data$logP <- -log10(data$P)

# 简单曼哈顿图
pdf("manhattan.pdf", width=12, height=6)
plot(1:nrow(data), data$logP, 
     pch=20, col=as.numeric(factor(data$CHR))%%8+1,
     xlab="Genomic Position", ylab="-log10(P)",
     main="GWAS Results - Dominant Model")
abline(h=-log10(5e-8), col="red", lty=2)
abline(h=-log10(1e-5), col="blue", lty=2)
dev.off()

# QQ图
pdf("qqplot.pdf", width=6, height=6)
observed <- sort(-log10(data$P))
expected <- -log10(ppoints(length(observed)))
plot(expected, observed, 
     xlab="Expected -log10(P)", ylab="Observed -log10(P)",
     main="QQ Plot", pch=20)
abline(0, 1, col="red")
dev.off()

cat("图形已生成：manhattan.pdf 和 qqplot.pdf\n")
EOF

echo ""
echo "==============================================="
echo "完成！运行 'Rscript plot_results.R' 生成图形"
echo "==============================================="
