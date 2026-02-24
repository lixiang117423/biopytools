biopytools fastq2vcf-gtx \
	-i ../../01.data/clean \
	-r OV53_1.primary_corrected.chr.fa \
	-p ./ \
	--read1-pattern-fastp "_1.clean.fq.gz" \
	--read2-pattern-fastp "_2.clean.fq.gz"
