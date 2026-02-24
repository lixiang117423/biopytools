biopytools vcf-sampler -i ../variation.filtered.snp.vcf.gz -o lab.snp.0.25.vcf.gz
biopytools admixture -i lab.snp.0.25.vcf.gz -o ./ -K 20 -s