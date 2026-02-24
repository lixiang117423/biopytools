filter_bam sample.clean.bam --nm 3 --threads 32 | samtools view - -b -@ 32 -o hic.filtered.bam
