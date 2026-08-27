"""phy2fa 命令行入口|phy2fa CLI entry"""
import argparse
import sys

from .calculator import Phy2FaCalculator
from .config import Phy2FaConfig
from .utils import ModuleLogger


def parse_arguments(argv=None):
    """解析命令行参数|Parse arguments"""
    parser = argparse.ArgumentParser(
        prog="phy2fa",
        description="Phylip 序列矩阵 → FASTA(支持 sequential/interleaved)"
                    "|Phylip sequence matrix → FASTA (sequential/interleaved)")
    parser.add_argument("-i", "--input", required=True,
                        help="Phylip 文件(.phy/.phylip/.gz)|Phylip file")
    parser.add_argument("-o", "--output-dir", default="./phy2fa_output",
                        help="输出目录(默认./phy2fa_output)|Output directory")
    parser.add_argument("--line-width", type=int, default=60,
                        help="FASTA 换行宽度,0=不换行(默认60)|FASTA line wrap "
                             "width, 0=no wrap (default 60)")
    parser.add_argument("--log-level", default="INFO",
                        help="日志级别(默认INFO)|Log level (default INFO)")
    return parser.parse_args(argv)


def main(argv=None):
    """主函数|Main function. Returns exit code."""
    args = parse_arguments(argv)
    try:
        config = Phy2FaConfig(
            input=args.input,
            output_dir=args.output_dir,
            line_width=args.line_width,
            log_level=args.log_level,
        )
        config.validate()
    except ValueError as e:
        print(f"参数错误|Parameter error: {e}", file=sys.stderr)
        return 1

    logger = ModuleLogger(log_file=str(config.log_file),
                          log_level=config.log_level).get_logger()
    return Phy2FaCalculator(config, logger).run()


if __name__ == "__main__":
    sys.exit(main())
