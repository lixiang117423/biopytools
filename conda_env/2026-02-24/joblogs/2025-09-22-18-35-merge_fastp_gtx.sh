# bash ./merge.fastq.sh

biopytools fastp \
    -i ~/project/96.yabian/13.bsa/01.data/raw \
    -o ~/project/96.yabian/13.bsa/01.data/clean \
    --read1-suffix "_1.fq" --read2-suffix "_2.fq"

biopytools gtx \
    -i ~/project/96.yabian/13.bsa/01.data/clean \
    -o ~/project/96.yabian/13.bsa/02.gtx \
    -r ~/project/96.yabian/13.bsa/01.data/genome/genome.fa \
    -t 64 --joint-threads 64 --enable-joint