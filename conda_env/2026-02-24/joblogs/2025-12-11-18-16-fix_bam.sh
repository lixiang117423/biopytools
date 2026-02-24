samtools view -@ 88 -u -f 1 -F 2316 sample.clean.bam | samtools sort -@ 88 -n -o fixed.namesorted.bam
