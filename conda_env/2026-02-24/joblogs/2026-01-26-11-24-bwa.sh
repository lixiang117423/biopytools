biopytools bwa -i 01.data/clean -o 02.挂载的基因组 -g 01.data/genome/41_RagTag_scaffolded.fa  -p _1.clean.gz -t 64
biopytools bwa -i 01.data/clean -o 03.未挂载的基因组 -g 01.data/genome/41_RagTag_unscaffolded.fa  -p _1.clean.gz -t 64
