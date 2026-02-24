bcftools view -@ 64 -S 88.第一批369个样品.txt variation.filtered.snp.vcf.gz -Oz -o filtered_369samples.vcf.gz

tabix -@ 64 -p vcf ./filtered_369samples.vcf.gz

biopytools filter-snp-indel -i ./filtered_369samples.vcf.gz -o ./