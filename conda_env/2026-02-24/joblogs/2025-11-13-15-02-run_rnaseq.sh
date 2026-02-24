biopytools rnaseq \
    -g 01.data/genome/MSU.fa \
    -f 01.data/genome/MSU.IGDBv1.Allset.gtf \
    -i 01.data/clean \
    -o 03.output_rice \
    -p "*_1.clean.fq.gz" -t 80
