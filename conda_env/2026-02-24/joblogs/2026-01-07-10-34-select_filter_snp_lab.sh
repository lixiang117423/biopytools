bcftools view -S 70.实验室测序的样品信息.txt chr.snp.vcf.gz -Oz -o chr.snp.lab.vcf.gz
biopytools filter-snp-indel -i chr.snp.lab.vcf.gz -o ./
