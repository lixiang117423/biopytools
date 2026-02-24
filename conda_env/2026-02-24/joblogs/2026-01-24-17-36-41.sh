mkdir -p 04.ragtag/41/hifiasm/共线性
mkdir -p 04.ragtag/41/hicanu/共线性
seqkit seq -j 64 -m 5000 01.data/clean/hifi/41.clean.fq.gz > 01.data/clean/hifi/41.clean.fq
biopytools hifi-hic -i 01.data/clean/hifi/41.clean.fq -o 02.hifiasm -t 64 -g 85m -p 41
biopytools ragtag -r 01.data/genome/Psoja_T2T.fa -q ~/project/28.大豆疫霉菌三代基因组/02.hifiasm/41/02.fasta/41.primary.fa -o 04.ragtag/41/hifiasm -s 41 -p hifiasm
biopytools genomesyn -i 04.ragtag/41/hifiasm/41_RagTag_scaffolded.fa -I 01.data/genome/Psoja_T2T.fa -o 04.ragtag/41/hifiasm/共线性 --min-length 100000
biopytools ragtag -r 01.data/genome/Psoja_T2T.fa -q ~/project/28.大豆疫霉菌三代基因组/02.hifiasm/41/02.fasta/41.hap1.fa -o 04.ragtag/41/hifiasm -s 41_hap1 -p hifiasm
biopytools genomesyn -i 04.ragtag/41/hifiasm/41_hap1_RagTag_scaffolded.fa -I 01.data/genome/Psoja_T2T.fa -o 04.ragtag/41/hifiasm/共线性 --min-length 100000
biopytools ragtag -r 01.data/genome/Psoja_T2T.fa -q ~/project/28.大豆疫霉菌三代基因组/02.hifiasm/41/02.fasta/41.hap2.fa -o 04.ragtag/41/hifiasm -s 41_hap2 -p hifiasm
biopytools genomesyn -i 04.ragtag/41/hifiasm/41_hap2_RagTag_scaffolded.fa -I 01.data/genome/Psoja_T2T.fa -o 04.ragtag/41/hifiasm/共线性 --min-length 100000
biopytools hicanu -i 01.data/clean/hifi/41.clean.fq -g 85m -p 41 -o 03.hicanu
biopytools ragtag -r 01.data/genome/Psoja_T2T.fa -q ~/project/28.大豆疫霉菌三代基因组/03.hicanu/41/02.fasta/41.contigs.fasta -o 04.ragtag/41/hicanu -s 41 -p hicanu
biopytools genomesyn -i 04.ragtag/41/hicanu/41_RagTag_scaffolded.fa -I 01.data/genome/Psoja_T2T.fa -o 04.ragtag/41/hicanu/共线性 --min-length 100000
