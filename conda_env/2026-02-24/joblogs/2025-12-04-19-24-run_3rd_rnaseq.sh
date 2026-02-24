minimap2 -ax splice -uf -k 14 --secondary=no -t 64 ../01.data/genome/MSU.fa SRR25203456_1.fastq.gz > SRR25203456.sam
samtools view -bS  SRR25203456.sam > SRR25203456.bam
samtools sort -@ 88 -o  SRR25203456.sord.bam SRR25203456.bam
