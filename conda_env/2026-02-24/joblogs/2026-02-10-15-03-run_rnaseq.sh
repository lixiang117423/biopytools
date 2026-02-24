biopytools fastp -i 01.data/raw -o 01.data/clean --read1-suffix "_1.fq.gz"  --read2-suffix "_2.fq.gz" 
biopytools rnaseq -i 01.data/clean -o 02.rnaseq -g 01.data/genome/OV53_Chr.fa -f 01.data/genome/OV53_Chr.gtf -t 64
