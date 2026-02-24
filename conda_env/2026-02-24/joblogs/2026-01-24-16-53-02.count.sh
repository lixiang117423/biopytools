biopytools kmc count \
    -d ../../01.data/clean \
    -o database \
    -k 51 -t 64 \
    --read1-suffix _1.clean.fq.gz \
    --read2-suffix _2.clean.fq.gz

biopytools kmc matrix -i database -o database -t 64
