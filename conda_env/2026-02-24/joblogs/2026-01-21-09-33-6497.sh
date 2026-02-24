mkdir -p 04.ragtag/6497/hifiasm/共线性
mkdir -p 04.ragtag/6497/hicanu/共线性
biopytools hifi-hic -i 01.data/clean/hifi/6497.clean.fq.gz -o 02.hifiasm -t 64 -g 85m -p 6497
biopytools ragtag -r 01.data/genome/Psoja_T2T.fa -q ~/project/28.大豆疫霉菌三代基因组/02.hifiasm/6497/02.fasta/6497.primary.fa -o 04.ragtag/6497/hifiasm -s 6497 -p hifiasm
biopytools genomesyn -i 04.ragtag/6497/hifiasm/6497_RagTag_scaffolded.fa -I 01.data/genome/Psoja_T2T.fa -o 04.ragtag/6497/hifiasm/共线性 --min-length 100000
biopytools hicanu -i 01.data/clean/hifi/6497.clean.fq.gz -g 85m -p 6497 -o 03.hicanu
biopytools ragtag -r 01.data/genome/Psoja_T2T.fa -q ~/project/28.大豆疫霉菌三代基因组/03.hicanu/6497/02.fasta/6497.contigs.fasta -o 04.ragtag/6497/hicanu -s 6497 -p hicanu
biopytools genomesyn -i 04.ragtag/6497/hicanu/6497_RagTag_scaffolded.fa -I 01.data/genome/Psoja_T2T.fa -o 04.ragtag/6497/hicanu/共线性 --min-length 100000
