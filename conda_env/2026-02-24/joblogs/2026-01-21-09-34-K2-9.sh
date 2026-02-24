mkdir -p 04.ragtag/K2-9/hifiasm/共线性
mkdir -p 04.ragtag/K2-9/hicanu/共线性
biopytools hifi-hic -i 01.data/clean/hifi/K2-9.clean.fq.gz -o 02.hifiasm -t 64 -g 85m -p K2-9
biopytools ragtag -r 01.data/genome/Psoja_T2T.fa -q ~/project/28.大豆疫霉菌三代基因组/02.hifiasm/K2-9/02.fasta/K2-9.primary.fa -o 04.ragtag/K2-9/hifiasm -s K2-9 -p hifiasm
biopytools genomesyn -i 04.ragtag/K2-9/hifiasm/K2-9_RagTag_scaffolded.fa -I 01.data/genome/Psoja_T2T.fa -o 04.ragtag/K2-9/hifiasm/共线性 --min-length 100000
biopytools hicanu -i 01.data/clean/hifi/K2-9.clean.fq.gz -g 85m -p K2-9 -o 03.hicanu
biopytools ragtag -r 01.data/genome/Psoja_T2T.fa -q ~/project/28.大豆疫霉菌三代基因组/03.hicanu/K2-9/02.fasta/K2-9.contigs.fasta -o 04.ragtag/K2-9/hicanu -s K2-9 -p hicanu
biopytools genomesyn -i 04.ragtag/K2-9/hicanu/K2-9_RagTag_scaffolded.fa -I 01.data/genome/Psoja_T2T.fa -o 04.ragtag/K2-9/hicanu/共线性 --min-length 100000
