biopytools vcf-sampler -i ../chr.snp.lab.vcf.gz -o lab.snp.0.05.vcf.gz -r 0.05
biopytools admixture -i lab.snp.0.05.vcf.gz -o ./ -K 20 -s