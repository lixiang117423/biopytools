bcftools view -S 88.第一批369个样品.txt /share/org/YZWL/yzwl_lixg/project/21.野生大豆群体/05.所有样品/14.实验室的样品/chr.snp.lab.vcf.gz -Oz -o filtered_369samples.vcf.gz --threads 64 --force-samples

tabix -@ 64 -p vcf ./filtered_369samples.vcf.gz

biopytools filter-snp-indel -i ./filtered_369samples.vcf.gz -o ./

biopytools plink-gwas -i ./variation.filtered.snp.biallelic.vcf.gz -p ./88.第一批369个样品表型-病斑长度均值.txt  -T quantitative -m all -o ./01.plink --no-strat-corr -t 64

biopytools tassel-gwas -i ./variation.filtered.snp.biallelic.vcf.gz -p ./88.第一批369个样品表型-病斑长度均值.txt -o ./02.tassel -m BOTH -t 64

biopytools gemma-gwas -i ./variation.filtered.snp.biallelic.vcf.gz -p ./88.第一批369个样品表型-病斑长度均值.txt -o ./03.gemma --no-qc