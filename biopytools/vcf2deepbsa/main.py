"""vcf2deepbsa 命令行入口|vcf2deepbsa CLI entry"""
import argparse
import sys

from .config import Vcf2DeepBsaConfig
from .converter import Vcf2DeepBsaConverter
from .utils import ModuleLogger


def parse_arguments(argv=None):
    """解析命令行参数|Parse arguments"""
    parser = argparse.ArgumentParser(
        prog="vcf2deepbsa",
        description="VCF转DeepBSA输入CSV(提取FORMAT中的AD)|"
                    "VCF to DeepBSA input CSV (extract AD from FORMAT)")
    parser.add_argument("-i", "--input-vcf", required=True,
                        help="输入 VCF 文件(支持 .gz)|Input VCF file (.gz supported)")
    parser.add_argument("-o", "--output-dir", default="./vcf2deepbsa_output",
                        help="输出目录(默认: ./vcf2deepbsa_output)|"
                             "Output directory (default: ./vcf2deepbsa_output)")
    parser.add_argument("--log-level", default="INFO",
                        help="日志级别(默认: INFO)|Log level (default: INFO)")
    return parser.parse_args(argv)


def main(argv=None) -> int:
    """主函数,返回退出码|Main function, returns exit code"""
    args = parse_arguments(argv)
    try:
        config = Vcf2DeepBsaConfig(
            input_vcf=args.input_vcf,
            output_dir=args.output_dir,
            log_level=args.log_level,
        )
        config.validate()
    except ValueError as e:
        print(f"参数错误|Parameter error: {e}", file=sys.stderr)
        return 1

    log_file = str(config.output_path / "vcf2deepbsa.log")
    logger = ModuleLogger(log_file=log_file, log_level=config.log_level).get_logger()
    logger.info("=" * 60)
    logger.info("vcf2deepbsa 启动|vcf2deepbsa started")

    converter = Vcf2DeepBsaConverter(config, logger)
    csv_path = converter.run()

    logger.info(f"输出CSV|Output CSV: {csv_path}")
    logger.info("vcf2deepbsa 完成|vcf2deepbsa finished")
    return 0


if __name__ == "__main__":
    sys.exit(main())
