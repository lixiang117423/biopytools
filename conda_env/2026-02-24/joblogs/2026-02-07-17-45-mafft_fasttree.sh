cat output/*/*.plastome.fasta > mafft_tree/all.plastome.fa

biopytools mafft-fasttree -i mafft_tree/all.plastome.fa -o mafft_tree -t 64
