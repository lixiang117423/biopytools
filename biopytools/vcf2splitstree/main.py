"""vcf2splitstree 命令行入口|vcf2splitstree CLI entry"""
import argparse
import sys

from .calculator import Vcf2SplitstreeCalculator
from .config import Vcf2SplitstreeConfig
from .utils import ModuleLogger


def parse_arguments(argv=None):
    """解析命令行参数|Parse arguments"""
    parser = argparse.ArgumentParser(
        prog="vcf2splitstree",
        description="VCF → SplitsTree6 距离矩阵 CSV|VCF to SplitsTree6 distance "
                    "matrix CSV")
    parser.add_argument("-i", "--input", required=True,
                        help="VCF 变异文件(.vcf/.vcf.gz)|VCF variants file")
    parser.add_argument("-o", "--output-dir", default="./vcf2splitstree_output",
                        help="输出目录(默认./vcf2splitstree_output)|Output directory")
    parser.add_argument("--log-level", default="INFO",
                        help="日志级别(默认INFO)|Log level (default INFO)")
    return parser.parse_args(argv)


def main(argv=None):
    """主函数|Main function. Returns exit code."""
    args = parse_arguments(argv)
    try:
        config = Vcf2SplitstreeConfig(
            input=args.input,
            output_dir=args.output_dir,
            log_level=args.log_level,
        )
        config.validate()
    except ValueError as e:
        print(f"参数错误|Parameter error: {e}", file=sys.stderr)
        return 1

    logger = ModuleLogger(log_file=str(config.log_file),
                          log_level=config.log_level).get_logger()
    return Vcf2SplitstreeCalculator(config, logger).run()


if __name__ == "__main__":
    sys.exit(main())
