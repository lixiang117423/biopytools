"""
序列提取命令行入口|Sequence extraction CLI entry point
"""

import argparse
import sys

from .config import SeqExtractConfig
from .utils import SeqExtractLogger, SeqExtractRunner


def parse_arguments():
    """解析命令行参数|Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="序列提取工具(seqkit封装,自动识别ID/ID文件/BED)|Sequence extraction tool (seqkit wrapper, auto-detect ID/ID file/BED)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="示例|Examples: biopytools seq-extract -i gene.id.txt -s gene.fa -o gene.genomic.fa",
    )

    parser.add_argument(
        "-i", "--input",
        required=True,
        help="查询:单个ID、ID文件(一列)或BED文件(>=2列)|Query: single ID, ID file (1 column), or BED file (>=2 columns)",
    )
    parser.add_argument(
        "-s", "--sequence",
        required=True,
        help="目标序列FASTA文件|Target sequence FASTA file",
    )
    parser.add_argument(
        "-o", "--output",
        default=None,
        help="输出文件(默认自动推导:{query}.{subject}.fa)|Output file (default: auto-derived {query}.{subject}.fa)",
    )
    parser.add_argument(
        "--bed",
        action="store_true",
        dest="force_bed",
        help="强制BED模式(跳过自动检测)|Force BED mode (skip auto-detection)",
    )

    return parser.parse_args()


def main():
    """主函数|Main function"""
    args = parse_arguments()

    try:
        config = SeqExtractConfig(
            input_query=args.input,
            sequence_file=args.sequence,
            output_file=args.output,
            force_bed=args.force_bed,
        )
        config.validate()

        logger_manager = SeqExtractLogger()
        logger = logger_manager.get_logger()

        runner = SeqExtractRunner(config, logger)
        success = runner.run()

        sys.exit(0 if success else 1)

    except ValueError as e:
        print(f"错误|Error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"未知错误|Unknown error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
