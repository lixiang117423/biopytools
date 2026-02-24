biopytools fastq2vcf-gtx \
	-i 01.data/raw \
	-g 01.data/genome/OV53_Chr.fa \
    -t 64 \
	--read1-pattern-fastp "_r1.fq.gz" \
	--read2-pattern-fastp "_r2.fq.gz"