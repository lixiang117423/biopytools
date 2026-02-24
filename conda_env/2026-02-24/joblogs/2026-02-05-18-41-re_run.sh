mkdir 01.data/ngs/6497
cp ~/project/19.大豆疫霉菌/01.data/clean/6497_* 01.data/ngs/6497
biopytools hifi-hic -i 01.data/clean/hifi/6497.clean.fq.gz --ngs 01.data/ngs/6497 -o 02.hifiasm -t 64 -g 85m -p 6497
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
mkdir 01.data/ngs/K2-9
cp ~/project/19.大豆疫霉菌/01.data/clean/K2-9_* 01.data/ngs/K2-9
biopytools hifi-hic -i 01.data/clean/hifi/K2-9.clean.fq.gz --ngs 01.data/ngs/K2-9 -o 02.hifiasm -t 64 -g 85m -p K2-9
mkdir -p 03.ragtag/K2-9/共线性
biopytools ragtag -r 01.data/genome/Psoja_T2T.fa -q 02.hifiasm/K2-9/03.ngs_polish/04.reassembly/02.fasta/K2-9.primary.fa -o 03.ragtag/K2-9 -s K2-9 -p hifiasm
biopytools genomesyn -i 03.ragtag/K2-9/K2-9_RagTag_scaffolded.fa -I 01.data/genome/Psoja_T2T.fa -o 03.ragtag/K2-9/共线性 --min-length 100000
biopytools ragtag -r 01.data/genome/Psoja_T2T.fa -q 02.hifiasm/K2-9/03.ngs_polish/04.reassembly/02.fasta/K2-9.hap1.fa -o 03.ragtag/K2-9 -s K2-9_hap1 -p hifiasm
biopytools genomesyn -i 03.ragtag/K2-9/K2-9_hap1_RagTag_scaffolded.fa -I 01.data/genome/Psoja_T2T.fa -o 03.ragtag/K2-9/共线性 --min-length 100000
biopytools ragtag -r 01.data/genome/Psoja_T2T.fa -q 02.hifiasm/K2-9/03.ngs_polish/04.reassembly/02.fasta/K2-9.hap2.fa -o 03.ragtag/K2-9 -s K2-9_hap2 -p hifiasm
biopytools genomesyn -i 03.ragtag/K2-9/K2-9_hap2_RagTag_scaffolded.fa -I 01.data/genome/Psoja_T2T.fa -o 03.ragtag/K2-9/共线性 --min-length 100000
mkdir 04.busco/K2-9/primary
biopytools busco -i 03.ragtag/K2-9/K2-9_RagTag_scaffolded.fa -l ~/database/busco/stramenopiles_odb12 -o 04.busco/K2-9/primary -t 64 --offline
mkdir 04.busco/K2-9/hap1
biopytools busco -i 03.ragtag/K2-9/K2-9_hap1_RagTag_scaffolded.fa -l ~/database/busco/stramenopiles_odb12 -o 04.busco/K2-9/hap1 -t 64 --offline
mkdir 04.busco/K2-9/hap2
biopytools busco -i 03.ragtag/K2-9/K2-9_hap2_RagTag_scaffolded.fa -l ~/database/busco/stramenopiles_odb12 -o 04.busco/K2-9/hap2 -t 64 --offline
mkdir 01.data/ngs/N3-3
cp ~/project/19.大豆疫霉菌/01.data/clean/N3-3_* 01.data/ngs/N3-3
biopytools hifi-hic -i 01.data/clean/hifi/N3-3.clean.fq.gz --ngs 01.data/ngs/N3-3 -o 02.hifiasm -t 64 -g 85m -p N3-3
mkdir -p 03.ragtag/N3-3/共线性
biopytools ragtag -r 01.data/genome/Psoja_T2T.fa -q 02.hifiasm/N3-3/03.ngs_polish/04.reassembly/02.fasta/N3-3.primary.fa -o 03.ragtag/N3-3 -s N3-3 -p hifiasm
biopytools genomesyn -i 03.ragtag/N3-3/N3-3_RagTag_scaffolded.fa -I 01.data/genome/Psoja_T2T.fa -o 03.ragtag/N3-3/共线性 --min-length 100000
biopytools ragtag -r 01.data/genome/Psoja_T2T.fa -q 02.hifiasm/N3-3/03.ngs_polish/04.reassembly/02.fasta/N3-3.hap1.fa -o 03.ragtag/N3-3 -s N3-3_hap1 -p hifiasm
biopytools genomesyn -i 03.ragtag/N3-3/N3-3_hap1_RagTag_scaffolded.fa -I 01.data/genome/Psoja_T2T.fa -o 03.ragtag/N3-3/共线性 --min-length 100000
biopytools ragtag -r 01.data/genome/Psoja_T2T.fa -q 02.hifiasm/N3-3/03.ngs_polish/04.reassembly/02.fasta/N3-3.hap2.fa -o 03.ragtag/N3-3 -s N3-3_hap2 -p hifiasm
biopytools genomesyn -i 03.ragtag/N3-3/N3-3_hap2_RagTag_scaffolded.fa -I 01.data/genome/Psoja_T2T.fa -o 03.ragtag/N3-3/共线性 --min-length 100000
mkdir 04.busco/N3-3/primary
biopytools busco -i 03.ragtag/N3-3/N3-3_RagTag_scaffolded.fa -l ~/database/busco/stramenopiles_odb12 -o 04.busco/N3-3/primary -t 64 --offline
mkdir 04.busco/N3-3/hap1
biopytools busco -i 03.ragtag/N3-3/N3-3_hap1_RagTag_scaffolded.fa -l ~/database/busco/stramenopiles_odb12 -o 04.busco/N3-3/hap1 -t 64 --offline
mkdir 04.busco/N3-3/hap2
biopytools busco -i 03.ragtag/N3-3/N3-3_hap2_RagTag_scaffolded.fa -l ~/database/busco/stramenopiles_odb12 -o 04.busco/N3-3/hap2 -t 64 --offline
