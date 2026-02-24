biopytools parabricks \
  -i /share/org/YZWL/yzwl_lixg/project/16.荠菜/15.each_8_samples \
  -o /share/org/YZWL/yzwl_lixg/project/16.荠菜/15.each_8_samples/mapping \
  -r /share/org/YZWL/yzwl_lixg/project/16.荠菜/15.each_8_samples/genome.fa \
  -t 64 \
  --read1-pattern "*_1.fq" \
  --read2-pattern "*_2.fq"

biopytools gatk-joint \
    -i /share/org/YZWL/yzwl_lixg/project/16.荠菜/15.each_8_samples/mapping/vcf \
    -o /share/org/YZWL/yzwl_lixg/project/16.荠菜/15.each_8_samples/04.filtered_snp_indel \
    -r /share/org/YZWL/yzwl_lixg/project/16.荠菜/15.each_8_samples/genome.fa