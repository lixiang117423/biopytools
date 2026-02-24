biopytools fastq2vcf-parabricks \
    -i /share/org/YZWL/yzwl_lixg/project/19.大豆疫霉菌/97.测试全自动流程-gtx/01.data/raw \
    -r /share/org/YZWL/yzwl_lixg/project/19.大豆疫霉菌/01.data/genome/Phytophthora_sojae_JS2_genome.fasta \
    -p ./ \
    --read1-pattern-fastp "_1.clean.fq.gz" \
    --read2-pattern-fastp "_2.clean.fq.gz" \
    --skip-qc