mkdir 01.data/ngs/6497
cp ~/project/19.大豆疫霉菌/01.data/clean/6497_* 01.data/ngs/6497
biopytools hifi-hic -i 01.data/raw/hifi/fastq/6497.hifi_reads.fastq.gz --ngs 01.data/ngs/6497 -o 02.hifiasm -t 64 -g 85m -p 6497
mkdir -p 03.ragtag/6497/共线性
biopytools ragtag -r 01.data/genome/Psoja_T2T.fa -q 02.hifiasm/6497/03.ngs_polish/04.reassembly/02.fasta/6497.primary.fa -o 03.ragtag/6497 -s 6497 -p hifiasm
biopytools genomesyn -i 03.ragtag/6497/6497_RagTag_scaffolded.fa -I 01.data/genome/Psoja_T2T.fa -o 03.ragtag/6497/共线性 --min-length 100000
biopytools ragtag -r 01.data/genome/Psoja_T2T.fa -q 02.hifiasm/6497/03.ngs_polish/04.reassembly/02.fasta/6497.hap1.fa -o 03.ragtag/6497 -s 6497_hap1 -p hifiasm
biopytools genomesyn -i 03.ragtag/6497/6497_hap1_RagTag_scaffolded.fa -I 01.data/genome/Psoja_T2T.fa -o 03.ragtag/6497/共线性 --min-length 100000
biopytools ragtag -r 01.data/genome/Psoja_T2T.fa -q 02.hifiasm/6497/03.ngs_polish/04.reassembly/02.fasta/6497.hap2.fa -o 03.ragtag/6497 -s 6497_hap2 -p hifiasm
biopytools genomesyn -i 03.ragtag/6497/6497_hap2_RagTag_scaffolded.fa -I 01.data/genome/Psoja_T2T.fa -o 03.ragtag/6497/共线性 --min-length 100000
mkdir 04.busco/6497/primary
biopytools busco -i 03.ragtag/6497/6497_RagTag_scaffolded.fa -l ~/database/busco/stramenopiles_odb12 -o 04.busco/6497/primary -t 64 --offline
mkdir 04.busco/6497/hap1
biopytools busco -i 03.ragtag/6497/6497_hap1_RagTag_scaffolded.fa -l ~/database/busco/stramenopiles_odb12 -o 04.busco/6497/hap1 -t 64 --offline
mkdir 04.busco/6497/hap2
biopytools busco -i 03.ragtag/6497/6497_hap2_RagTag_scaffolded.fa -l ~/database/busco/stramenopiles_odb12 -o 04.busco/6497/hap2 -t 64 --offline
