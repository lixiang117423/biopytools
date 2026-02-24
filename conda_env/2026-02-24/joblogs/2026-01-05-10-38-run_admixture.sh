# biopytools vcf-sampler -i chr.snp.vcf.gz -o chr.snp.10.vcf.gz -r 0.1
biopytools vcf-sampler -i chr.snp.10.vcf.gz -o chr.snp.10.2.vcf.gz -r 0.1
biopytools admixture -i chr.snp.10.2.vcf.gz -o ./ -K 30 -s