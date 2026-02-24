mkdir -p 04.ragtag/N3-3/hifiasm/共线性
mkdir -p 04.ragtag/N3-3/hicanu/共线性
biopytools hifi-hic -i 01.data/clean/hifi/N3-3.clean.fq.gz -o 02.hifiasm -t 64 -g 85m -p N3-3
biopytools ragtag -r 01.data/genome/Psoja_T2T.fa -q ~/project/28.大豆疫霉菌三代基因组/02.hifiasm/N3-3/02.fasta/N3-3.primary.fa -o 04.ragtag/N3-3/hifiasm -s N3-3 -p hifiasm
biopytools genomesyn -i 04.ragtag/N3-3/hifiasm/N3-3_RagTag_scaffolded.fa -I 01.data/genome/Psoja_T2T.fa -o 04.ragtag/N3-3/hifiasm/共线性 --min-length 100000
biopytools hicanu -i 01.data/clean/hifi/N3-3.clean.fq.gz -g 85m -p N3-3 -o 03.hicanu
biopytools ragtag -r 01.data/genome/Psoja_T2T.fa -q ~/project/28.大豆疫霉菌三代基因组/03.hicanu/N3-3/02.fasta/N3-3.contigs.fasta -o 04.ragtag/N3-3/hicanu -s N3-3 -p hicanu
biopytools genomesyn -i 04.ragtag/N3-3/hicanu/N3-3_RagTag_scaffolded.fa -I 01.data/genome/Psoja_T2T.fa -o 04.ragtag/N3-3/hicanu/共线性 --min-length 100000
