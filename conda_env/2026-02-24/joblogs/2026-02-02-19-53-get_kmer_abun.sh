biopytools kmertools intersect -m kmer_matrix_with_header.txt -k hap1_group2_chr12_kmer_51.fa -t 64 -o hap1_group2_chr12_matrix.txt
biopytools kmertools intersect -m kmer_matrix_with_header.txt -k hap2_group9_chr12_kmer_51.fa -t 64 -o hap2_group9_chr12_matrix.txt
biopytools kmertools compare -f1 hap1_group2_chr12_matrix.txt -f2 hap2_group9_chr12_matrix.txt  -o unique
