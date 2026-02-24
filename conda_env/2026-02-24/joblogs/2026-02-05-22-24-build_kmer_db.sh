biopytools kmertools build -i ngs -o kmer_db/20260206 -t 64
biopytools kmertools kmer2vcf -i kmer_db/20260206/kmer_matrix_with_header.txt -o kmer_db/20260206/kmer.vcf.gz -t 64