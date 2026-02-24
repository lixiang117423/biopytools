# ~/software/scripts/59.ALLHiC_Juicer_3DDNA流程.sh \
~/software/scripts/65.ALLHiC_Asmkit流程挂载染色体.sh \
    -r /share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/21.ALLHiC_hap/test_pipeline/OV53_1.primary.fa \
    -1 /share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/21.ALLHiC_hap/test_pipeline/fastq/OV53_1-hic_R1.fastq.gz \
    -2 /share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/21.ALLHiC_hap/test_pipeline/fastq/OV53_1-hic_R2.fastq.gz \
    -w /share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/21.ALLHiC_hap/test_pipeline/output \
    -k 12 \
    -t 88 \
    --memory 500G \
    --skip-mapping --skip-allele \
    --skip-prune --skip-partition \
    --skip-extract --skip-rescue \
    --skip-optimize --skip-build \
    --skip-plot