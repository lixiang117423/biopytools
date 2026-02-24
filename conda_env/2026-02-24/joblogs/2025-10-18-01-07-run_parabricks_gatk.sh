biopytools fastp \
    -i /share/org/YZWL/yzwl_lixg/database/yimeijun/resequence/raw \
    -o /share/org/YZWL/yzwl_lixg/database/yimeijun/resequence/clean \
    --read1-suffix "_1.fastq.gz" \
    --read2-suffix "_2.fastq.gz"

biopytools parabricks \
  -i /share/org/YZWL/yzwl_lixg/database/yimeijun/resequence/clean \
  -o /share/org/YZWL/yzwl_lixg/project/17.大豆疫霉基因变异_李磊/03.mapping \
  -r /share/org/YZWL/yzwl_lixg/project/17.大豆疫霉基因变异_李磊/01.data/GCF_000149755.1_P.sojae_V3.0_genomic_modified.fna \
  -t 64

biopytools gatk-joint \
    -i /share/org/YZWL/yzwl_lixg/project/17.大豆疫霉基因变异_李磊/03.mapping/vcf \
    -o /share/org/YZWL/yzwl_lixg/project/17.大豆疫霉基因变异_李磊/04.filtered_snp_indel \
    -r /share/org/YZWL/yzwl_lixg/project/17.大豆疫霉基因变异_李磊/01.data/GCF_000149755.1_P.sojae_V3.0_genomic_modified.fna

biopytools filter-snp-indel \
    -i 03.mapping/vcf/combined.g.vcf \
    -o 04.filtered_snp_indel \
    -t 88