"""
FOF生成主程序模块|FOF Generation Main Module
"""

import argparse
import sys

from .config import FofConfig
from .utils import FofLogger
from .generator import FofGenerator


def parse_arguments():
    """解析命令行参数|Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='生成样品名到文件路径的FOF映射表|Generate sample-to-file-path FOF mapping table',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例|Examples:
  %(prog)s -i ./data/ -o samples.fof
  %(prog)s -i ./data/ -o samples.fof -s .fastq.gz -s .fq.gz
  %(prog)s -i sample.bam -o sample.fof
  %(prog)s -i ./data/ -o samples.fof -r
        '''
    )

    parser.add_argument('-i', '--input',
                        required=True,
                        help='输入文件或目录|Input file or directory')

    parser.add_argument('-o', '--output',
                        required=True,
                        help='输出FOF文件路径|Output FOF file path')

    parser.add_argument('-s', '--suffix',
                        action='append',
                        default=None,
                        help='文件后缀过滤（可多次指定，如 -s .fastq.gz -s .fq.gz）|'
                             'File suffix filter (can be specified multiple times)')

    parser.add_argument('-r', '--recursive',
                        action='store_true',
                        default=False,
                        help='递归扫描子目录|Recursively scan subdirectories')

    return parser.parse_args()


def main():
    """主函数|Main function"""
    args = parse_arguments()

    try:
        # 创建配置|Create configuration
        config = FofConfig(
            input_path=args.input,
            output_file=args.output,
            suffixes=args.suffix if args.suffix else [],
            recursive=args.recursive
        )
        config.validate()

        # 初始化日志|Initialize logging
        logger_manager = FofLogger(config.output_file_obj.parent)
        logger = logger_manager.get_logger()

        logger.info(f"输入路径|Input path: {config.input_path}")
        logger.info(f"输出文件|Output file: {config.output_file}")
        if config.suffixes:
            logger.info(f"后缀过滤|Suffix filter: {', '.join(config.suffixes)}")
        if config.recursive:
            logger.info(f"递归扫描|Recursive scan: 是|Yes")

        # 创建生成器并执行|Create generator and execute
        generator = FofGenerator(config, logger)
        success = generator.generate()

        if success:
            sys.exit(0)
        else:
            sys.exit(1)

    except ValueError as e:
        print(f"配置错误|Configuration error: {e}")
        sys.exit(1)

    except Exception as e:
        print(f"程序执行出错|Program execution error: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    main()
