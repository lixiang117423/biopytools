biopytools gatk-joint \
    -i /share/org/YZWL/yzwl_lixg/project/94.rice_gas/02.mapping/vcf \
    -o /share/org/YZWL/yzwl_lixg/project/94.rice_gas/03.gatk_joint \
    -r /share/org/YZWL/yzwl_lixg/project/94.rice_gas/01.data/genome/MSU.fa

biopytools filter-snp-indel \
    -i /share/org/YZWL/yzwl_lixg/project/94.rice_gas/03.gatk_joint/joint_genotyping_merged_filtered.vcf.gz \
    -o /share/org/YZWL/yzwl_lixg/project/94.rice_gas/04.filtered_snp_indel \
    -t 64 \
    --snp-dp 5 \
    --indel-dp 5