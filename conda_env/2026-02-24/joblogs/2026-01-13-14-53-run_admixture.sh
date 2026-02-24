biopytools vcf-sampler -i ../variation.filtered.snp.biallelic.vcf.gz -o lab.snp.025.vcf.gz

tabix -@ 64 -p vcf ./lab.snp.025.vcf.gz
biopytools admixture -i lab.snp.025.vcf.gz -o ./ -K 20 -s