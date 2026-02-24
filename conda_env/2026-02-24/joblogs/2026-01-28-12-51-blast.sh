seqkit sample -j 12 ../41.unused.fa -p 0.0001 > sampled_reads.fa

blastn -query sampled_reads.fa -db ~/database/nt/nt -out blast_result.txt -outfmt 6 -evalue 1e-5 -num_threads 2