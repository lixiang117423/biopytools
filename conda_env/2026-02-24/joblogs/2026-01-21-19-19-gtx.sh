biopytools fastq2vcf-gtx \
    -i 01.data/raw \
    -g 01.data/genome/EcA_RagTag_scaffolded.fa \
    -o ./ \
    --read1-pattern-fastp _1.fq.gz \
    --read2-pattern-fastp _2.fq.gz \
    -t 64
