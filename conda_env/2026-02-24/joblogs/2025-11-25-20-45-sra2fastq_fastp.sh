biopytools sra2fastq \
    -i raw/metagenome/sra \
    -o raw/metagenome/fastq 

biopytools fastp \
    -i raw/metagenome/fastq \
    -o clean/metagenome \
    --read1-suffix "_1.fastq.gz" \
    --read2-suffix "_2.fastq.gz"

python3 ~/software/scripts/28.run_myyc.py /share/org/YZWL/yzwl_lixg/project/94.rice_gas/04.mcyc/sample.txt

/share/org/YZWL/yzwl_lixg/software/scripts/31.run_kraken2.sh \
    -d /share/org/YZWL/yzwl_lixg/database/kraken2 \
    -i /share/org/YZWL/yzwl_lixg/project/94.rice_gas/01.data/clean/metagenome \
    -o /share/org/YZWL/yzwl_lixg/project/94.rice_gas/05.kraken2 \
    -t 88