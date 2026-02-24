bcftools view -S 90.第一次根据枝长筛选的样品.txt ../chr.snp.lab.vcf.gz -Oz -o chr.snp.lab.1st.vcf.gz --threads 64
biopytools filter-snp-indel -i chr.snp.lab.1st.vcf.gz -o ./