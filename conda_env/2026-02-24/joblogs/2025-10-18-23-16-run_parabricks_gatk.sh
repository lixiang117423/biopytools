# biopytools fastp \
#     -i /share/org/YZWL/yzwl_lixg/database/yimeijun/resequence/raw \
#     -o /share/org/YZWL/yzwl_lixg/database/yimeijun/resequence/clean \
#     --read1-suffix "_1.fastq.gz" \
#     --read2-suffix "_2.fastq.gz"

biopytools parabricks \
  -i /share/org/YZWL/yzwl_lixg/project/16.荠菜/01.data/clean \
  -o /share/org/YZWL/yzwl_lixg/project/16.荠菜/21.四倍体参数重新比对和检测变异_parabricks \
  -r /share/org/YZWL/yzwl_lixg/project/16.荠菜/01.data/genome/genome.fa \
  --ploidy 4 -t 64 

biopytools gatk-joint \
    -i /share/org/YZWL/yzwl_lixg/project/16.荠菜/21.四倍体参数重新比对和检测变异_parabricks/vcf \
    -o /share/org/YZWL/yzwl_lixg/project/16.荠菜/21.四倍体参数重新比对和检测变异_parabricks/gatk_filtered_snp_indel \
    -r /share/org/YZWL/yzwl_lixg/project/16.荠菜/01.data/genome/genome.fa

# biopytools filter-snp-indel \
#     -i /share/org/YZWL/yzwl_lixg/project/16.荠菜/21.四倍体参数重新比对和检测变异_parabricks/vcf/combined.g.vcf \
#     -o 04.filtered_snp_indel \
#     -t 88