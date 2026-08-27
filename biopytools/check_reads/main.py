"""check_reads 命令行入口|check_reads CLI entry"""
import argparse
import sys

from .calculator import CheckReadsCalculator
from .config import CheckReadsConfig
from .utils import ModuleLogger


def parse_arguments(argv=None):
    """解析命令行参数|Parse arguments"""
    parser = argparse.ArgumentParser(
        prog="check-reads",
        description="fastq 完整性检查(gz 校验/0字节/配对)|FASTQ integrity check "
                    "(gzip/empty/pairing)")
    parser.add_argument("-i", "--input", required=True,
                        help="fastq 目录(逗号分隔多个,递归扫描)|FASTQ dir(s), "
                             "comma-separated, recursive")
    parser.add_argument("-o", "--output-dir", default="./check_reads_output",
                        help="输出目录(默认./check_reads_output)|Output directory")
    parser.add_argument("-t", "--threads", type=int, default=12,
                        help="并行线程数(默认12)|Parallel threads (default 12)")
    parser.add_argument("--log-level", default="INFO",
                        help="日志级别(默认INFO)|Log level (default INFO)")
    return parser.parse_args(argv)


def main(argv=None):
    """主函数|Main function. Returns exit code."""
    args = parse_arguments(argv)
    try:
        config = CheckReadsConfig(
            input_dir=args.input,
            output_dir=args.output_dir,
            threads=args.threads,
            log_level=args.log_level,
        )
        config.validate()
    except ValueError as e:
        print(f"参数错误|Parameter error: {e}", file=sys.stderr)
        return 1

    logger = ModuleLogger(log_file=str(config.log_file),
                          log_level=config.log_level).get_logger()
    return CheckReadsCalculator(config, logger).run()


if __name__ == "__main__":
    sys.exit(main())
