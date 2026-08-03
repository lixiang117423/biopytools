"""
NeedleIdentity主程序模块|Needle Identity Main Module
"""

import argparse
import sys
from pathlib import Path

from .calculator import NeedleIdentityCalculator
from .config import NeedleIdentityConfig
from .utils import NeedleIdentityLogger, check_needle


def parse_arguments():
    """解析命令行参数|Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='序列两两identity计算(EMBOSS needle)|Pairwise sequence identity (EMBOSS needle)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='示例|Examples: %(prog)s -i sequences.fa -o output_dir/',
    )

    parser.add_argument('-i', '--input', required=True,
                        help='输入FASTA序列文件|Input FASTA file')
    parser.add_argument('-o', '--output-dir', default='./output',
                        help='输出目录|Output directory (default: ./output)')
    parser.add_argument('--needle-path', default=None,
                        help='needle可执行文件路径|needle executable path (default: ~/miniforge3/envs/needle_v.1.0.3/bin/needle)')
    parser.add_argument('--threads', type=int, default=12,
                        help='并行线程数|Threads (default: 12)')
    parser.add_argument('--gapopen', type=float, default=10.0,
                        help='gap开放罚分|Gap open penalty (default: 10.0)')
    parser.add_argument('--gapextend', type=float, default=0.5,
                        help='gap延伸罚分|Gap extend penalty (default: 0.5)')

    return parser.parse_args()


def main():
    """主函数|Main function"""
    args = parse_arguments()

    try:
        config = NeedleIdentityConfig(
            input_file=args.input,
            output_dir=args.output_dir,
            needle_path=args.needle_path,
            threads=args.threads,
            gapopen=args.gapopen,
            gapextend=args.gapextend,
        )
        config.validate()

        output_path = Path(config.output_dir)
        log_file = output_path / "99_logs" / f"{config.output_prefix}.needle_identity.log"

        logger_manager = NeedleIdentityLogger(log_file=str(log_file))
        logger = logger_manager.get_logger()

        logger.info("=" * 80)
        logger.info(" 序列两两identity计算|Pairwise Sequence Identity")
        logger.info("=" * 80)

        if not check_needle(config, logger):
            logger.error("EMBOSS不可用，请检查needle路径|EMBOSS unavailable, check needle path")
            sys.exit(1)

        calculator = NeedleIdentityCalculator(config, logger)
        result_file = calculator.run()

        logger.info(f"完成|Completed: {result_file}")
        logger.info("=" * 80)
        sys.exit(0)

    except Exception as e:
        print(f"错误|Error: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()
