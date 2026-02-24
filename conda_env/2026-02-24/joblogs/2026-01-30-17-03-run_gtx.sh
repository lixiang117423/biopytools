biopytools fastp -i 01.data/raw -o 01.data/clean --read1-suffix "_r1.fq.gz" --read2-suffix "_r2.fq.gz"

biopytools fastq2vcf-gtx \
	-i 01.data/clean \
	-g 01.data/genome/OV53_Chr.fa \
    -t 64