def generate_reverse_complement(seq):
    """生成反向互补序列"""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    return ''.join(complement[base] for base in reversed(seq))

def expand_kmer_library(input_fasta, input_bed, output_fasta, output_bed):
    """扩展k-mer库，添加反向互补序列"""
    # 读取原始库，生成扩展版本
    # 处理去重逻辑
    # 输出扩展后的文件

def is_self_complement(seq):
    """检查序列是否与自身反向互补相同"""
    return seq == generate_reverse_complement(seq)