biopytools kmc count \
    -d fastq \
    -o database \
    -k 51 -t 64 \
    --read1-suffix _1.clean.fq \
    --read2-suffix _2.clean.fq

biopytools kmc matrix -i database -o database -t 24
