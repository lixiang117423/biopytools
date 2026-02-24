mkdir -p 04.ragtag/41/hifiasm/共线性
mkdir -p 04.ragtag/41/hicanu/共线性
biopytools hifi-hic -i 01.data/clean/hifi/41.clean.fq.gz -o 02.hifiasm -t 64 -g 85m -p 41
biopytools ragtag -r 01.data/genome/Psoja_T2T.fa -q ~/project/28.大豆疫霉菌三代基因组/02.hifiasm/41/02.fasta/41.primary.fa -o 04.ragtag/41/hifiasm -s 41 -p hifiasm
biopytools genomesyn -i 04.ragtag/41/hifiasm/41_RagTag_scaffolded.fa -I 01.data/genome/Psoja_T2T.fa -o 04.ragtag/41/hifiasm/共线性 --min-length 100000
