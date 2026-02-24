biopytools fastp \
    -i /share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/01.data/2nd_rnaseq/raw \
    -o /share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/01.data/2nd_rnaseq/clean \
    -t 88

biopytools rnaseq \
    -g /share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/69.二代转录组比对到hifi基因组上/Orychophragmus_violaceus_OV53_1_HiFi.fa \
    -f /share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/69.二代转录组比对到hifi基因组上/complete.genomic.gtf \
    -i /share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/01.data/2nd_rnaseq/clean \
    -o /share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/69.二代转录组比对到hifi基因组上 \
    -p "*_1.clean.fq.gz" -t 88