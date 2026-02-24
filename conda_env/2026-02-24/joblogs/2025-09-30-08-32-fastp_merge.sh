#  biopytools fastp \
#     -i /share/org/YZWL/yzwl_lixg/project/95.gdl_metagenome/01.data/raw/fastq_data \
#     -o /share/org/YZWL/yzwl_lixg/project/95.gdl_metagenome/01.data/clean \
#     --read1-suffix "_1.fastq.gz"
#     --read2-suffix "_2.fastq.gz"

# zcat /share/org/YZWL/yzwl_lixg/project/95.gdl_metagenome/01.data/clean/*_1.clean.fq.gz >> /share/org/YZWL/yzwl_lixg/project/95.gdl_metagenome/01.data/all_1.fq
# zcat /share/org/YZWL/yzwl_lixg/project/95.gdl_metagenome/01.data/clean/*_2.clean.fq.gz >> /share/org/YZWL/yzwl_lixg/project/95.gdl_metagenome/01.data/all_2.fq

metawrap assembly \
    -1 /share/org/YZWL/yzwl_lixg/project/95.gdl_metagenome/01.data/all_1.fq \
    -2 /share/org/YZWL/yzwl_lixg/project/95.gdl_metagenome/01.data/all_2.fq \
    -o /share/org/YZWL/yzwl_lixg/project/95.gdl_metagenome/02.assembly \
    -m 300 -t 80 --megahit 