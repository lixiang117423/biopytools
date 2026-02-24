# 先计算 HWE
  bcftools +hwe \
    -O z \
    -o tmp_hwe.vcf.gz \
    variation.filtered.snp.vcf.gz

  # 然后应用所有过滤条件
  bcftools view \
    -i 'TYPE="snp" && N_ALT=1 && INFO/DP>5 && INFO/DP<50 && MAF>0.05 && F_MISSING<=0.05 && HWE>1e-6' \
    -O z \
    -o variation.filtered.snp.strict.vcf.gz \
    tmp_hwe.vcf.gz

  # 索引
  bcftools index variation.filtered.snp.strict.vcf.gz
