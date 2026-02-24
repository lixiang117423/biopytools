biopytools plink-gwas \
    -i chr.snp.182.samples.maf005.vcf.gz \
    -p 114.第一批369个有表型的样品去冗余做GWAS-样品病害登记.txt \
    -o 01.plink \
    -T quantitative \
    -m all \
    --no-strat-corr \
    -t 64
