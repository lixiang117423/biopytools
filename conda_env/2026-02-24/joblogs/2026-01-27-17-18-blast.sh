seqkit fq2fa -j 12 41.unused.fq.gz > 41.unused.fa
blastn -query 41.unused.fa -db ~/database/nt/nt -out ./blast/41.unused.nt.blast.txt -outfmt 6 -evalue 1e-5 -num_threads 8
