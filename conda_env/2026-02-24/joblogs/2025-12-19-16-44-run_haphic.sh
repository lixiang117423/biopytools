biopytools haphic \
    -a OV53_1.primary.fasta \
    -1 OV53_1-hic_R1.fastq.gz \
    -2 OV53_1-hic_R1.fastq.gz \
    -c 12 \
    --threads 64 --processes 64 \
    --memory-limit 500G -o ./ \
    --continue-from reassign
