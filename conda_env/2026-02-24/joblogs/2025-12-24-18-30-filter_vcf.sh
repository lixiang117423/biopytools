bcftools view \
    -i 'TYPE="snp" && N_ALT=1 && INFO/DP>5 && INFO/DP<50 && MAF>0.05 && F_MISSING<=0.05' \
    -O z \
    -o variation.filtered.snp.strict.vcf.gz \
    variation.filtered.snp.vcf.gz

bcftools index variation.filtered.snp.strict.vcf.gz
