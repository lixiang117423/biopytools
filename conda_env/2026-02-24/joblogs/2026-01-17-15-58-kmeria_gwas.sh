biopytools kmeria pipeline \
    -i ../../01.data/clean/lab/ \
    --sample 115.kmer-gwas的样品名称.txt \
    -d 117.kmer-gwas的测序深度.txt \
    -p 116.kmer-gwas的分级表型.txt \
    -o output -t 64 \
    --genome-file genome.fa \
    --gff-file genome.gff3
