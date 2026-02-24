~/software/scripts/32.run_bwa_gatk.sh \
    -i 01.data/clean \
    -r 01.data/genome/genome.fa \
    -o 02.mapping

biopytools gatk-joint \
    -i 02.mapping/02.gvcf \
    -r 01.data/genome/genome.fa \
    -o 03.gatk_joint

biopytools filter-snp-indel \
    -i 03.gatk_joint/joint_genotyping_merged_filtered.vcf.gz \
    -o 04.filtered_snp_indel \
    -t 12

python3 ~/software/scripts/23.bam_stats_reporter.py \
    -i 02.mapping/01.bam \
    -o 02.mapping/04.bam_stat/bam_stat.xlsx