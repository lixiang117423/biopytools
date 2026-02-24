biopytools fastq2vcf-gtx \
	-i ./raw/balanced \
	-r /share/org/YZWL/yzwl_lixg/project/21.野生大豆群体/01.data/genome/genome.fa \
	-p ./ \
	--read1-pattern-fastp "_1.fq.gz" \
	--read2-pattern-fastp "_2.fq.gz"
