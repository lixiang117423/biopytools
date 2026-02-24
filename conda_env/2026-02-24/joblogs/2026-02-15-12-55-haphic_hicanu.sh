biopytools haphic \
    -i 03.hicanu/EcA/02.fasta/EcA.contigs.fasta \
    -1 01.data/hic/clean/eca_hic_1.clean.fq.gz \
    -2 01.data/hic/clean/eca_hic_2.clean.fq.gz \
    -c 8 -o 05.haphic_hicanu -t 64 --prefix EcA
