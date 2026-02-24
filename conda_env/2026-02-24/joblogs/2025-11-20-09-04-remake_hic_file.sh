# 1. 对最终的 Scaffolds 建立索引 (这步很快)
samtools faidx yahs_out_scaffolds_final.fa

# 2. 提取新的、正确的名字和长度，生成 chrom.sizes.final
# 这一步生成的文件里，第一列就会变成 scaffold_1, scaffold_2 ...
cut -f1,2 yahs_out_scaffolds_final.fa.fai > chrom.sizes.final

# 3. 重新生成 .hic 文件
# 请确保 JUICER_JAR 的路径是对的
export JUICER_JAR="/share/org/YZWL/yzwl_lixg/software/juicer/scripts/common/juicer_tools.jar"

# 这一步可能需要几十分钟，取决于服务器性能
java -Xmx128G -jar ${JUICER_JAR} pre \
    alignments_sorted.txt \
    yahs_out_final_fixed.hic \
    chrom.sizes.final
