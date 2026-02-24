mkdir -p 06.二代比对的bam/41
mkdir -p 07.筛选的reads重新组装/01.data
mkdir -p 07.筛选的reads重新组装/02.hifiasm

biopytools bwa -i ../01.data/second/41/ -o 06.二代比对的bam/41 -g 02.hifiasm/41/02.fasta/41.primary.fa  -p _1.clean.fq.gz -t 64
biopytools coverage-filter -i 06.二代比对的bam/41/bam/41.bam -f 02.hifiasm/41/02.fasta/41.primary.fa -o 41 -d 06.二代比对的bam/41 --high-cov 95
grep -Fwf 06.二代比对的bam/41/41_high_quality.list 02.hifiasm/41/02.fasta/41.p_ctg.contig_reads.tsv | cut -f2 > 06.二代比对的bam/41/high_quality_read_names.txt
seqkit grep -n -f 06.二代比对的bam/41/high_quality_read_names.txt ../01.data/clean/hifi/41.clean.fq.gz -o 07.筛选的reads重新组装/01.data/41_high_quality_reads.fq.gz
biopytools hifi-hic -i 07.筛选的reads重新组装/01.data/41_high_quality_reads.fq.gz -o 07.筛选的reads重新组装/02.hifiasm -t 64 -g 85m -p 41

mkdir -p 07.筛选的reads重新组装/04.ragtag/41/hifiasm/共线性

biopytools ragtag -r ../01.data/genome/Psoja_T2T.fa -q 07.筛选的reads重新组装/02.hifiasm/41/02.fasta/41.primary.fa -o 07.筛选的reads重新组装/04.ragtag/41/hifiasm -s 41 -p hifiasm
biopytools genomesyn -i 07.筛选的reads重新组装/04.ragtag/41/hifiasm/41_RagTag_scaffolded.fa -I ../01.data/genome/Psoja_T2T.fa -o 07.筛选的reads重新组装/04.ragtag/41/hifiasm/共线性 --min-length 100000
biopytools ragtag -r ../01.data/genome/Psoja_T2T.fa -q 07.筛选的reads重新组装/02.hifiasm/41/02.fasta/41.hap1.fa -o 07.筛选的reads重新组装/04.ragtag/41/hifiasm -s 41_hap1 -p hifiasm
biopytools genomesyn -i 07.筛选的reads重新组装/04.ragtag/41/hifiasm/41_hap1_RagTag_scaffolded.fa -I ../01.data/genome/Psoja_T2T.fa -o 07.筛选的reads重新组装/04.ragtag/41/hifiasm/共线性 --min-length 100000
biopytools ragtag -r ../01.data/genome/Psoja_T2T.fa -q 07.筛选的reads重新组装/02.hifiasm/41/02.fasta/41.hap2.fa -o 07.筛选的reads重新组装/04.ragtag/41/hifiasm -s 41_hap2 -p hifiasm
biopytools genomesyn -i 07.筛选的reads重新组装/04.ragtag/41/hifiasm/41_hap2_RagTag_scaffolded.fa -I ../01.data/genome/Psoja_T2T.fa -o 07.筛选的reads重新组装/04.ragtag/41/hifiasm/共线性 --min-length 100000
