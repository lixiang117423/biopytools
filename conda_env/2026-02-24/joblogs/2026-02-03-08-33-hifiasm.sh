biopytools fastp -i 01.data/raw/hifi/68-1_hifi.fq.gz -o 01.data/clean/hifi/68-1_hifi.clean.fq.gz --single-end

biopytools hifi-hic -i 01.data/clean/hifi/68-1_hifi.clean.fq.gz -o 02.纯hifi组装 -t 64 -g 175m --ngs 01.data/clean/ngs --high-cov 90
