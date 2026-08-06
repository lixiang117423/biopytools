"""vcf2pav 命令行入口|vcf2pav CLI entry"""
import argparse
import sys

from .config import Vcf2PavConfig
from .converter import Vcf2PavConverter
from .utils import ModuleLogger


def parse_arguments(argv=None):
    """解析命令行参数|Parse arguments"""
    parser = argparse.ArgumentParser(
        prog="vcf2pav",
        description="VCF转PAV(Presence/Absence Variation)矩阵|"
                    "VCF to PAV (Presence/Absence Variation) matrix converter")
    parser.add_argument("-i", "--input", required=True,
                        help="输入 VCF 文件|Input VCF file")
    parser.add_argument("-o", "--output-dir", default="./vcf2pav_output",
                        help="输出目录|Output directory (default ./vcf2pav_output)")
    parser.add_argument("--log-level", default="INFO",
                        help="日志级别|Log level (default INFO)")
    return parser.parse_args(argv)


def main(argv=None):
    """主函数|Main function. Returns exit code."""
    args = parse_arguments(argv)
    try:
        config = Vcf2PavConfig(
            input_vcf=args.input,
            output_dir=args.output_dir,
            log_level=args.log_level,
        )
        config.validate()
    except ValueError as e:
        print(f"参数错误|Parameter error: {e}", file=sys.stderr)
        return 1

    log_file = str(config.output_path / "vcf2pav.log")
    logger = ModuleLogger(log_file=log_file, log_level=config.log_level).get_logger()
    logger.info("=" * 60)
    logger.info("vcf2pav 启动|vcf2pav started")
    logger.info(f"输入文件|Input: {config.input_vcf}")
    logger.info(f"输出目录|Output: {config.output_dir}")

    converter = Vcf2PavConverter(config, logger)
    output_path = converter.run()

    logger.info(f"PAV矩阵|PAV matrix: {output_path}")
    logger.info("vcf2pav 完成|vcf2pav finished")
    return 0


if __name__ == "__main__":
    sys.exit(main())
