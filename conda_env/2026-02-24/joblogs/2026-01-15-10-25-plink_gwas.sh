biopytools plink-gwas \
    -i input.vcf.gz \
    -p 122.第一批369个样品的病斑分级表型.txt \
    -o 01.plink \
    -T quantitative \
    -m all \
    --no-strat-corr \
    -t 64
