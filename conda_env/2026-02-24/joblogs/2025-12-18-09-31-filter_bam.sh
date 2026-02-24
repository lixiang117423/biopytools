filter_bam sample.clean.bam 1  --nm 3 --threads 32 | samtools view - -b -@ 32 -o hic.filtered.bam
