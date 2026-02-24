biopytools hicanu -i 01.data/clean/hifi/41.clean.fq.gz -g 85m -p 41 -o 07.hicanu组装的reads给hifiasm组装/01.hicanu
awk '$1 !~ /^#/ {print $2}' 07.hicanu组装的reads给hifiasm组装/01.hicanu/41/02.fasta/41.contig_reads.tsv > 07.hicanu组装的reads给hifiasm组装/01.hicanu/41/02.fasta/reads.id.txt
mkdir 07.hicanu组装的reads给hifiasm组装/02.hifiasm/41
seqkit grep -f 07.hicanu组装的reads给hifiasm组装/01.hicanu/41/02.fasta/reads.id.txt 01.data/clean/hifi/41.clean.fq.gz -j 64 -o 07.hicanu组装的reads给hifiasm组装/02.hifiasm/41/$s_hicanu_reads.fq.gz
biopytools hifi-hic -i 07.hicanu组装的reads给hifiasm组装/02.hifiasm/41/$s_hicanu_reads.fq.gz -o 07.hicanu组装的reads给hifiasm组装/02.hifiasm -t 64 -g 85m -p 41
