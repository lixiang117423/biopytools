"""seq-len 命令行入口|seq-len command-line entry"""

import argparse
import sys

from .config import SeqLenConfig
from .calculator import SeqLenCalculator
from .utils import SeqLenLogger


def parse_arguments(argv=None):
    """解析命令行参数|Parse CLI args"""
    parser = argparse.ArgumentParser(
        prog='seq-len',
        description='FASTA序列长度统计(文件/文件夹,合并+汇总)|'
                    'FASTA sequence length statistics (file/folder, merged + summary)')
    parser.add_argument('-i', '--input', required=True,
                        help='FASTA 文件或文件夹|FASTA file or folder')
    parser.add_argument('-o', '--output', required=True,
                        help='输出 TSV 路径或目录|Output TSV path or directory')
    parser.add_argument('--prefix', default=None,
                        help='输出前缀(默认取输入名)|Output prefix (default: input name)')
    parser.add_argument('--min-length', type=int, default=0,
                        help='最小长度过滤(0=不过滤)|Min length filter (0=no filter)')
    parser.add_argument('--sort', action='store_true',
                        help='按长度降序(默认保持输入顺序)|Sort by length descending')
    parser.add_argument('--no-summary', action='store_true',
                        help='不输出汇总表|Skip summary table')
    parser.add_argument('--log-file', default=None, help='日志文件|Log file')
    parser.add_argument('--log-level', default='INFO',
                        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
                        help='日志级别|Log level')
    parser.add_argument('-v', '--verbose', action='store_true', help='详细日志|Verbose')
    return parser.parse_args(argv)


def main():
    """主函数|Main entry"""
    args = parse_arguments()
    try:
        config = SeqLenConfig(
            input_path=args.input, output=args.output, prefix=args.prefix,
            min_length=args.min_length, sort=args.sort,
            summary=not args.no_summary,
            log_file=args.log_file, log_level=args.log_level, verbose=args.verbose)
        config.validate()

        logger = SeqLenLogger(config.log_file, config.log_level, config.verbose).get_logger()
        logger.info(f"输入|Input: {config.input_path} "
                    f"({'文件夹|folder' if config.is_folder else '文件|file'}, "
                    f"{len(config.input_files)} 个 FASTA|fastas)")
        logger.info(f"输出|Output: {config.tsv_path}")
        SeqLenCalculator(config, logger).run()
        sys.exit(0)
    except Exception as e:
        print(f"错误|Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
