biopytools hifi-hic \
    -i 01.data/clean/hifi/EcA.clean.fq.gz \
    -1 01.data/hic/clean//eca_hic_1.clean.fq.gz \
    -2 01.data/hic/clean//eca_hic_2.clean.fq.gz \
    -o 02.hifiasm -t 64 -g 180m -p EcA
