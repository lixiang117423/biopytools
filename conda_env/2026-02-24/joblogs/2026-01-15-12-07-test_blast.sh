 blastn -num_threads 16 \
    -taxids 4762 \
    -query 02.hifiasm/41/02.fasta/41.primary.fa \
    -db ~/database/nt/nt \
    -out 04.blast/41.hifiasm.primary.nt.blast.txt \
    -outfmt 6 -evalue 1e-5
