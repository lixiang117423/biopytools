"""
序列子集提取主程序 | Sequence Subsequence Extraction Main Program
"""

import argparse
import sys
from pathlib import Path
from .core import SequenceExtractor
from .utils import SubseqLogger


class Config:
    """配置类 | Configuration class"""
    def __init__(self):
        self.input_fasta = None
        self.output_fasta = None
        self.id_list_file = None
        self.keep_order = True

        # 模式匹配相关 | Pattern matching related
        self.pattern = None
        self.pattern_type = 'contains'  # 'contains', 'startswith', 'endswith', 'regex'
        self.case_sensitive = True

        # 长度筛选相关 | Length filtering related
        self.min_length = 0
        self.max_length = float('inf')


def parse_arguments():
    """解析命令行参数 | Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='🧬 序列子集提取工具 | Sequence Subsequence Extraction Tool',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例 | Examples:

  # 📋 根据ID列表提取序列 | Extract sequences by ID list
  python -m biopytools.subseq -i input.fasta -l id_list.txt -o output.fasta

  # 🔍 根据模式匹配提取序列 | Extract sequences by pattern matching
  python -m biopytools.subseq -i input.fasta -p "chr1_" -o output.fasta --pattern-type startswith

  # 📏 根据长度范围提取序列 | Extract sequences by length range
  python -m biopytools.subseq -i input.fasta -o output.fasta --min-length 1000 --max-length 5000

  # 🔄 忽略大小写匹配 | Case-insensitive matching
  python -m biopytools.subseq -i input.fasta -p "GENE" -o output.fasta --ignore-case
        """
    )

    # 必需参数 | Required parameters
    parser.add_argument('-i', '--input',
                        required=True,
                        help='📄 输入FASTA文件路径 | Input FASTA file path')

    parser.add_argument('-o', '--output',
                        required=True,
                        help='📤 输出FASTA文件路径 | Output FASTA file path')

    # 提取方式选择 | Extraction method selection (mutually exclusive)
    extraction_group = parser.add_mutually_exclusive_group(required=True)

    extraction_group.add_argument('-l', '--id-list',
                                 help='📋 ID列表文件路径 | ID list file path')

    extraction_group.add_argument('-p', '--pattern',
                                 help='🔍 模式匹配字符串 | Pattern matching string')

    extraction_group.add_argument('--length-only',
                                 action='store_true',
                                 help='📏 仅使用长度筛选 | Use only length filtering')

    # 模式匹配选项 | Pattern matching options
    pattern_group = parser.add_argument_group('模式匹配选项 | Pattern Matching Options')

    pattern_group.add_argument('--pattern-type',
                               choices=['contains', 'startswith', 'endswith', 'regex'],
                               default='contains',
                               help='📝 模式类型 | Pattern type (default: contains)')

    pattern_group.add_argument('--ignore-case',
                               action='store_true',
                               help='🔤 忽略大小写 | Case insensitive')

    # 长度筛选选项 | Length filtering options
    length_group = parser.add_argument_group('长度筛选选项 | Length Filtering Options')

    length_group.add_argument('--min-length',
                              type=int,
                              default=0,
                              help='📏 最小序列长度 | Minimum sequence length (default: 0)')

    length_group.add_argument('--max-length',
                              type=int,
                              help='📏 最大序列长度 | Maximum sequence length (default: unlimited)')

    # 其他选项 | Other options
    parser.add_argument('--no-order',
                       action='store_true',
                       help='🔄 不保持ID列表顺序 | Do not keep ID list order')

    parser.add_argument('--log-dir',
                       default='.',
                       help='📁 日志输出目录 | Log output directory (default: current directory)')

    return parser.parse_args()


def main():
    """主函数 | Main function"""
    try:
        # 解析参数 | Parse arguments
        args = parse_arguments()

        # 配置对象 | Configuration object
        config = Config()
        config.input_fasta = args.input
        config.output_fasta = args.output
        config.keep_order = not args.no_order

        # 根据提取方式设置配置 | Set configuration based on extraction method
        if args.id_list:
            config.id_list_file = args.id_list
            extraction_method = 'id_list'
        elif args.pattern:
            config.pattern = args.pattern
            config.pattern_type = args.pattern_type
            config.case_sensitive = not args.ignore_case
            extraction_method = 'pattern'
        elif args.length_only:
            config.min_length = args.min_length
            config.max_length = args.max_length if args.max_length else float('inf')
            extraction_method = 'length'
        else:
            print("❌ 必须指定一种提取方式 | Must specify an extraction method", file=sys.stderr)
            sys.exit(1)

        # 设置日志 | Setup logging
        log_dir = Path(args.log_dir)
        logger_manager = SubseqLogger(log_dir)
        logger = logger_manager.get_logger()

        logger.info("🧬 序列子集提取工具启动 | Sequence Subsequence Extraction Tool Started")
        logger.info(f"📝 提取方式 | Extraction method: {extraction_method}")

        # 创建序列提取器 | Create sequence extractor
        extractor = SequenceExtractor(config, logger)

        # 执行提取 | Perform extraction
        success = False
        if extraction_method == 'id_list':
            success = extractor.extract_sequences_by_id_list()
        elif extraction_method == 'pattern':
            success = extractor.extract_sequences_by_pattern()
        elif extraction_method == 'length':
            success = extractor.extract_sequences_by_length()

        if success:
            logger.info("✅ 序列提取完成 | Sequence extraction completed successfully")
            sys.exit(0)
        else:
            logger.error("❌ 序列提取失败 | Sequence extraction failed")
            sys.exit(1)

    except KeyboardInterrupt:
        print("\n🛑 用户中断操作 | User interrupted operation", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"❌ 程序错误 | Program error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()