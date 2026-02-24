biopytools fastq2vcf-gtx \
	-i ./01.data/raw/merge \
	-r /share/org/YZWL/yzwl_lixg/project/18.小花糖芥/06.BSA/01.data/genome/genome.fa \
	-p ./ \
	--read1-pattern-fastp "_1.fq.gz" \
	--read2-pattern-fastp "_2.fq.gz"
