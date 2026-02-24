biopytools sra2fastq \
    -i /share/org/YZWL/yzwl_lixg/database/yimeijun/resequence/sra \
    -o /share/org/YZWL/yzwl_lixg/database/yimeijun/resequence/raw

biopytools fastp \
    -i /share/org/YZWL/yzwl_lixg/database/yimeijun/resequence/raw \
    -o /share/org/YZWL/yzwl_lixg/database/yimeijun/resequence/clean \
    --read1-suffix "_1.fastq.gz" \
    --read2-suffix "_2.fastq.gz"