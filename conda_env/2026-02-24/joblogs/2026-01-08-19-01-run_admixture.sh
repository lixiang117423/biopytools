tabix -@ 64 -p vcf ./variation.filtered.snp.vcf.gz
biopytools vcf-sampler -i ./variation.filtered.snp.vcf.gz -o lab.snp.0.1.vcf.gz

tabix -@ 64 -p vcf ./lab.snp.0.1.vcf.gz
biopytools admixture -i lab.snp.0.1.vcf.gz -o ./ -K 20 -s