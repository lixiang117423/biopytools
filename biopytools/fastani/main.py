"""fastANI命令行入口|fastANI CLI Entry

示例|Examples: biopytools fastani -i genome_dir/ -o output_dir/
"""

import argparse
import os
import sys

from .calculator import FastaniCalculator
from .config import FastaniConfig
from .utils import FastaniLogger


def parse_arguments(argv=None):
    """解析命令行参数|Parse command-line arguments"""
    parser = argparse.ArgumentParser(
        prog='biopytools fastani',
        description='fastANI全基因组ANI计算:两基因组"拼写相似度平均分",'
                    '95%+通常同种|fastANI whole-genome ANI: 95%+ usually same species',
        epilog='示例|Examples: biopytools fastani -i genome_dir/ -o output_dir/')
    parser.add_argument('-i', '--input', default=None,
                        help='all-vs-all输入(目录/列表文件/单fasta,与-q/-r互斥)|'
                             'all-vs-all input (dir/list/single fasta)')
    parser.add_argument('-q', '--query', default=None,
                        help='query侧输入(定向模式)|Query side (directional mode)')
    parser.add_argument('-r', '--reference', default=None,
                        help='reference侧输入(定向模式)|Reference side (directional mode)')
    parser.add_argument('-o', '--output-dir', default='./fastani_output',
                        help='输出目录|Output directory')
    parser.add_argument('-t', '--threads', type=int, default=12,
                        help='线程数|Thread count')
    parser.add_argument('-k', '--kmer', type=int, default=16,
                        help='k-mer大小(<=16)|K-mer size (<=16)')
    parser.add_argument('--frag-len', type=int, default=3000,
                        help='片段长度|Fragment length')
    parser.add_argument('--min-fraction', type=float, default=0.2,
                        help='信任ANI的最小共享比例|Min shared fraction to trust ANI')
    parser.add_argument('--iterated', dest='iterated', action='store_true',
                        default=True, help=argparse.SUPPRESS)
    parser.add_argument('--no-iterated', dest='iterated', action='store_false',
                        help='关闭大数据集自动遍历(强制all-vs-all)|Disable auto '
                             'iterated mode (force all-vs-all)')
    parser.add_argument('--iterated-threshold', type=int, default=100,
                        help='触发遍历的基因组数阈值(默认100)|Genome count '
                             'threshold for iterated mode (default 100)')
    parser.add_argument('--ref-batch-size', type=int, default=50,
                        help='遍历模式 reference 分批大小(默认50,越小内存越低)'
                             '|Reference batch size in iterated mode (default '
                             '50; smaller = lower memory)')
    parser.add_argument('--log-level', default='INFO',
                        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
                        help='日志级别|Log level')
    return parser.parse_args(argv)


def main():
    """主函数|Main function"""
    args = parse_arguments()
    try:
        config = FastaniConfig(
            input=args.input, query=args.query, reference=args.reference,
            output_dir=args.output_dir, threads=args.threads, kmer=args.kmer,
            frag_len=args.frag_len, min_fraction=args.min_fraction,
            iterated=args.iterated, iterated_threshold=args.iterated_threshold,
            ref_batch_size=args.ref_batch_size, log_level=args.log_level)
        config.validate()

        logs_dir = os.path.join(config.output_dir, '99_logs')
        logger = FastaniLogger(logs_dir, config.log_level).get_logger()
        logger.info(
            f"fastANI模块启动|fastANI module started "
            f"(mode={'all-vs-all' if config.all_vs_all else 'query-vs-ref'}, "
            f"threads={config.threads})")

        calculator = FastaniCalculator(config, logger)
        ok = calculator.run()
        sys.exit(0 if ok else 1)
    except ValueError as e:
        print(f"错误|Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
