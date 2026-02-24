biopytools kmertools build -i 01.data -o 02.kmer_db -t 64
biopytools kmertools kmer2vcf -i 02.kmer_db/kmer_matrix_with_header.txt -o 02.kmer_db/kmer.vcf.gz -t 64