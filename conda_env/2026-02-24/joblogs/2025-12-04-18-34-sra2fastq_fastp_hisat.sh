# biopytools sra2fastq \
#     -i /share/org/YZWL/yzwl_lixg/project/94.rice_gas/18.根部转录组/01.data/sra \
#     -o /share/org/YZWL/yzwl_lixg/project/94.rice_gas/18.根部转录组/01.data/fastq \
#     -t 88

# biopytools fastp \
#     -i /share/org/YZWL/yzwl_lixg/project/94.rice_gas/18.根部转录组/01.data/fastq \
#     -o /share/org/YZWL/yzwl_lixg/project/94.rice_gas/18.根部转录组/01.data/clean \
#     -t 88 --read1-suffix "_1.fastq.gz" --read2-suffix "_2.fastq.gz"

biopytools rnaseq \
    -g /share/org/YZWL/yzwl_lixg/project/94.rice_gas/18.根部转录组/01.data/genome/MSU.fa \
    -f /share/org/YZWL/yzwl_lixg/project/94.rice_gas/18.根部转录组/01.data/genome/MSU.gtf \
    -i /share/org/YZWL/yzwl_lixg/project/94.rice_gas/18.根部转录组/01.data/clean \
    -o /share/org/YZWL/yzwl_lixg/project/94.rice_gas/18.根部转录组/02.mapping \
    -t 88 -p "*_1.clean.fq.gz"