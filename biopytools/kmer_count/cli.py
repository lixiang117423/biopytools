"""
🖥️ 命令行接口 | Command Line Interface
"""

import sys
import os
import argparse
from .config import KmerCountConfig
from .main import KmerCountAnalyzer

def create_parser():
    """🔧 创建命令行参数解析器 | Create command line argument parser"""
    parser = argparse.ArgumentParser(
        description='🧬 K-mer丰度分析工具 | K-mer abundance analysis tool',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
🌟 示例用法 | Example usage:
  %(prog)s -i /data/fastq -p "*_1.fq.gz" -k kmers.fasta -o results/
  %(prog)s -i ./samples -p "*_R1.fastq" -k kmers.fasta -b kmers.bed -w 500000 -o results/
        """
    )
    
    # 必需参数 | Required arguments
    required = parser.add_argument_group('📋 必需参数 | Required arguments')
    # required.add_argument('-i', '--input', required=True,
    #                     help='📁 FASTQ文件输入目录 | FASTQ files input directory')
    # required.add_argument('-p', '--pattern', required=True,
    #                     help='📁 FASTQ文件模式，如*_1.fq.gz | FASTQ file pattern, e.g. *_1.fq.gz')
    
    required.add_argument('-i', '--input', required=True,
                    help='📁 输入文件目录 | Input files directory')
    required.add_argument('-p', '--pattern', required=True,
                        help='📁 文件模式，支持FASTQ和FASTA格式，如*_1.fq.gz、*.fasta、*.fa | File pattern, support FASTQ and FASTA formats, e.g. *_1.fq.gz, *.fasta, *.fa')

    required.add_argument('-k', '--kmer-lib', required=True,
                        help='🧬 K-mer库文件(FASTA格式) | K-mer library file (FASTA format)')
    required.add_argument('-o', '--output', required=True,
                        help='📂 输出目录 | Output directory')
    
    # 可选参数 | Optional arguments
    optional = parser.add_argument_group('⚙️ 可选参数 | Optional arguments')
    optional.add_argument('-b', '--bed-file',
                        help='📋 BED文件路径 | BED file path')
    optional.add_argument('-m', '--kmer-size', type=int, default=51,
                        help='📏 K-mer长度 (默认: %(default)s) | K-mer size (default: %(default)s)')
    optional.add_argument('-s', '--hash-size', default='1000M',
                        help='🗂️ 哈希表大小 (默认: %(default)s) | Hash table size (default: %(default)s)')
    optional.add_argument('-t', '--threads', type=int, default=8,
                        help='🧵 线程数 (默认: %(default)s) | Number of threads (default: %(default)s)')
    optional.add_argument('-w', '--window-size', type=int, default=500000,
                        help='🪟 滑动窗口大小bp (默认: %(default)s) | Sliding window size in bp (default: %(default)s)')
    optional.add_argument('--step-size', type=int,
                        help='👣 滑动窗口步长bp (默认: window-size/5) | Sliding window step size in bp (default: window-size/5)')
    optional.add_argument('-C', '--canonical', action='store_true',
                        help='🔄 统计正向和反向互补链 | Count both forward and reverse complement')
    optional.add_argument('--keep-temp', action='store_true',
                        help='💾 保留临时文件 | Keep temporary files')
    optional.add_argument('--keep-binary', action='store_true',
                        help='🔢 保留0/1存在缺失矩阵 | Keep 0/1 presence/absence matrix')
    optional.add_argument('--jellyfish-path', default='jellyfish',
                        help='🐙 Jellyfish程序路径 (默认: %(default)s) | Jellyfish program path (default: %(default)s)')
    optional.add_argument('-v', '--verbose', action='store_true',
                        help='📝 详细输出 | Verbose output')
    
    return parser


def main():
    """🚀 主函数 - 这是run_kmer_count脚本的入口点 | Main function - Entry point for run_kmer_count script"""
    parser = create_parser()
    args = parser.parse_args()
    
    try:
        # 创建配置 | Create configuration
        config = KmerCountConfig(args)
        
        # 创建分析器并运行 | Create analyzer and run
        analyzer = KmerCountAnalyzer(config)
        analyzer.run_analysis()
        
    except KeyboardInterrupt:
        print("\n⏹️ 分析被用户中断 | Analysis interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"\n💥 错误 | Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
