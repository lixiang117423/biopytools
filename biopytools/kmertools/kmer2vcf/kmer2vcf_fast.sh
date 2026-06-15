#!/bin/bash
# kmer matrix to VCF fast pipeline (awk-based, single-pass streaming)
# Usage: bash kmer2vcf_fast.sh input.zst output.vcf.gz [chr_length] [kmer_length]
# Example: bash kmer2vcf_fast.sh input.txt.zst output.vcf.gz 100000000 51

set -euo pipefail

INPUT="${1:?Usage: $0 <input.zst> <output.vcf.gz> [chr_length] [kmer_length]}"
OUTPUT="${2:?Usage: $0 <input.zst> <output.vcf.gz> [chr_length] [kmer_length]}"
CHR_LEN="${3:-100000000}"
KMER_LEN="${4:-51}"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] 开始转换 | Starting conversion"
echo "  输入 | Input:    $INPUT"
echo "  输出 | Output:   $OUTPUT"
echo "  染色体长度 | Chr length: $CHR_LEN"
echo "  Kmer长度 | Kmer length: $KMER_LEN"

# 根据文件后缀选择解压方式 | Choose decompression based on file extension
decompress() {
    case "$INPUT" in
        *.zst)  zstd -d -c "$INPUT" ;;
        *.gz)   zcat "$INPUT" ;;
        *)      cat "$INPUT" ;;
    esac
}

# 解压 + awk转换 + bgzip压缩 | Decompress + awk convert + bgzip compress
# 所有操作通过管道流式处理，不产生中间文件 | All operations stream through pipe, no temp files
decompress | awk -v chr_len="$CHR_LEN" -v kmer_len="$KMER_LEN" '
BEGIN {
    chr = 1
    pos = 0
    n_samples = 0
}
NR == 1 {
    # 写VCF header | Write VCF header
    printf("##fileformat=VCFv4.1\n")
    printf("##INFO=<ID=KL,Number=1,Type=Integer,Description=\"Kmer length (%d)\">\n", kmer_len)
    printf("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")

    # 解析样本名 | Parse sample names
    # 如果第一列是 KMER/ID 等标识符，从第2列开始取样本名
    # If first column is KMER/ID identifier, start sample names from column 2
    start_col = 1
    if ($1 == "KMER" || $1 == "ID" || $1 == "KMERID" || $1 == "KMER_ID") {
        start_col = 2
    }

    printf("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
    for (i = start_col; i <= NF; i++) {
        printf("\t%s", $i)
        n_samples++
    }
    printf("\n")
    next
}
{
    # 跳过不足2列的行 | Skip lines with less than 2 columns
    if (NF < 2) next

    # 分配染色体位置 | Assign chromosome position
    pos++
    if (pos > chr_len) {
        chr++
        pos = 1
    }

    # 输出VCF行 | Output VCF line
    printf("chr%d\t%d\t%s\tA\tG\t100\tPASS\tKL=%d\tGT", chr, pos, $1, kmer_len)
    for (i = 2; i <= NF; i++) {
        printf("\t%s/%s", $i, $i)
    }
    printf("\n")

    # 进度显示 | Progress display (every 1M lines)
    if (NR % 1000000 == 0) {
        printf("[%dM lines processed, chr%d, pos%d]\n", NR/1000000, chr, pos) > "/dev/stderr"
    }
}
END {
    printf("\n转换完成 | Conversion completed: %d lines, %d chromosomes\n", NR-1, chr) > "/dev/stderr"
}
' | bgzip -c > "$OUTPUT"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] 完成 | Done: $OUTPUT"
