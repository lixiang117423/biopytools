"""
FASTA序列ID重命名主程序|FASTA ID Renamer Main Program
"""

import argparse
import sys
from pathlib import Path

from .config import FastaIDRenamerConfig
from .processor import FastaIDRenamer
from .utils import RenamerLogger


def parse_arguments():
    """解析命令行参数|Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="FASTA序列ID重命名工具 - 按顺序重命名所有序列|FASTA ID Renamer - Rename all sequences in order",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例|Examples:
  # 基本用法 - 重命名所有序列|Basic usage - rename all sequences
  python fasta_id_renamer.py -i input.fa -o output.fa

  # 指定染色体数量，自动提取|Specify chromosome count for auto-extraction
  python fasta_id_renamer.py -i input.fa -o output.fa -n 20

  # 自定义前缀和格式|Custom prefix and format
  python fasta_id_renamer.py -i input.fa -o output.fa -p chromosome --no-zero-padding

  # 三位数字填充|Three-digit padding
  python fasta_id_renamer.py -i input.fa -o output.fa -w 3
        """
    )

    # 必需参数|Required parameters
    parser.add_argument('-i', '--input',
                       required=True,
                       help='输入FASTA文件|Input FASTA file')

    parser.add_argument('-o', '--output',
                       required=True,
                       help='输出FASTA文件|Output FASTA file')

    # 重命名规则|Renaming rules
    parser.add_argument('-p', '--prefix',
                       default='Chr',
                       help='序列前缀(默认: Chr)|Sequence prefix (default: Chr)')

    parser.add_argument('--no-zero-padding',
                       action='store_true',
                       help='不使用零填充(如Chr1而非Chr01)|Do not use zero padding (Chr1 instead of Chr01)')

    parser.add_argument('-w', '--padding-width',
                       type=int,
                       default=2,
                       help='填充宽度(默认: 2)|Padding width (default: 2)')

    # 染色体提取|Chromosome extraction
    parser.add_argument('-n', '--chr-count',
                       type=int,
                       default=0,
                       help='染色体数量，提取前N条作为染色体(默认: 0)|Chromosome count, extract first N as chromosomes (default: 0)')

    # ID映射|ID mapping
    parser.add_argument('--no-mapping',
                       action='store_true',
                       help='不保存ID映射文件|Do not save ID mapping file')

    parser.add_argument('--mapping-file',
                       default=None,
                       help='ID映射文件路径|ID mapping file path')

    # 其他选项|Other options
    parser.add_argument('--log-level',
                       choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
                       default='INFO',
                       help='日志级别(默认: INFO)|Log level (default: INFO)')

    return parser.parse_args()


def main():
    """主函数|Main function"""
    args = parse_arguments()

    try:
        # 创建配置|Create configuration
        config = FastaIDRenamerConfig(
            input_file=args.input,
            output_file=args.output,
            prefix=args.prefix,
            use_zero_padding=not args.no_zero_padding,
            padding_width=args.padding_width,
            chr_count=args.chr_count,
            save_mapping=not args.no_mapping,
            mapping_file=args.mapping_file
        )

        # 验证配置|Validate configuration
        config.validate()

        # 创建日志器|Create logger
        log_file = Path(args.output).parent / "fasta_id_renamer.log"
        logger_manager = RenamerLogger(log_file=str(log_file), log_level=args.log_level)
        logger = logger_manager.get_logger()

        logger.info("FASTA序列ID重命名工具启动|FASTA ID Renamer started")
        logger.info(f"版本|Version: 1.0.0")

        # 创建处理器并执行|Create processor and run
        renamer = FastaIDRenamer(config, logger)
        renamer.process()

        logger.info("程序正常退出|Program exited successfully")
        sys.exit(0)

    except ValueError as e:
        print(f"参数错误|Parameter error: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"错误|Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
