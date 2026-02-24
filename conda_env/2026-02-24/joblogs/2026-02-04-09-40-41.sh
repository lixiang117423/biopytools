mkdir 01.data/ngs/41
cp ~/project/19.大豆疫霉菌/01.data/clean/41_* 01.data/ngs/41
biopytools hifi-hic -i 01.data/clean/hifi/41.clean.fq.gz --ngs 01.data/ngs/41 -o 02.hifiasm -t 64 -g 85m -p 41
mkdir -p 03.ragtag/41/共线性
biopytools ragtag -r 01.data/genome/Psoja_T2T.fa -q 02.hifiasm/41/03.ngs_polish/04.reassembly/02.fasta/41.primary.fa -o 03.ragtag/41 -s 41 -p hifiasm
biopytools genomesyn -i 03.ragtag/41/41_RagTag_scaffolded.fa -I 01.data/genome/Psoja_T2T.fa -o 03.ragtag/41/共线性 --min-length 100000
biopytools ragtag -r 01.data/genome/Psoja_T2T.fa -q 02.hifiasm/41/03.ngs_polish/04.reassembly/02.fasta/41.hap1.fa -o 03.ragtag/41 -s 41_hap1 -p hifiasm
biopytools genomesyn -i 03.ragtag/41/41_hap1_RagTag_scaffolded.fa -I 01.data/genome/Psoja_T2T.fa -o 03.ragtag/41/共线性 --min-length 100000
biopytools ragtag -r 01.data/genome/Psoja_T2T.fa -q 02.hifiasm/41/03.ngs_polish/04.reassembly/02.fasta/41.hap2.fa -o 03.ragtag/41 -s 41_hap2 -p hifiasm
biopytools genomesyn -i 03.ragtag/41/41_hap2_RagTag_scaffolded.fa -I 01.data/genome/Psoja_T2T.fa -o 03.ragtag/41/共线性 --min-length 100000
mkdir 04.busco/41/primary
biopytools busco -i 03.ragtag/41/41_RagTag_scaffolded.fa -l ~/database/busco/stramenopiles_odb12 -o 04.busco/41/primary -t 64 --offline
mkdir 04.busco/41/hap1
biopytools busco -i 03.ragtag/41/41_hap1_RagTag_scaffolded.fa -l ~/database/busco/stramenopiles_odb12 -o 04.busco/41/hap1 -t 64 --offline
mkdir 04.busco/41/hap2
biopytools busco -i 03.ragtag/41/41_hap2_RagTag_scaffolded.fa -l ~/database/busco/stramenopiles_odb12 -o 04.busco/41/hap2 -t 64 --offline
