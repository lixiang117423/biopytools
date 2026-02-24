biopytools plink-gwas \
    -i chr.snp.182.samples.maf005.vcf.gz \
    -p 111.第一批369个有表型的样品去冗余做GWAS-样品病斑平均长度.txt \
    -o 01.plink \
    -T quantitative \
    -m all \
    --no-strat-corr \
    -t 64
