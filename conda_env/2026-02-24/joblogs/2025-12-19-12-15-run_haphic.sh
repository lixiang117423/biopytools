biopytools haphic \
    -a genome.fa \
    -1 FZY4201_1.fq.gz \
    -2 FZY4201_2.fq.gz \
    -c 8 \
    --threads 64 --processes 64 \
    --memory-limit 500G -o ./
