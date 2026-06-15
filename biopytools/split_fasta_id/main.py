"""
FASTA ID分割主程序模块|FASTA ID Splitting Main Module
"""

import argparse
import sys
from .config import SplitConfig
from .utils import SplitLogger
from .processor import FastaProcessor

class FastaIDSplitter:
    """FASTA ID分割主类|Main FASTA ID Splitter Class"""

    def __init__(self, **kwargs):
        # 初始化配置|Initialize configuration
        self.config = SplitConfig(**kwargs)
        self.config.validate()

        # 初始化日志|Initialize logging
        self.logger_manager = SplitLogger(self.config.output_file)
        self.logger = self.logger_manager.get_logger()

        # 初始化处理器|Initialize processor
        self.processor = FastaProcessor(self.config, self.logger)

    def run_splitting(self):
        """运行完整的FASTA ID分割流程|Run complete FASTA ID splitting pipeline"""
        try:
            self.logger.info("开始FASTA ID分割流程|Starting FASTA ID splitting pipeline")
            self.logger.info("=" * 60)
            self.logger.info(f"输入文件|Input file: {self.config.input_file}")
            self.logger.info(f"输出文件|Output file: {self.config.output_file}")
            self.logger.info(f"提取位置|Extract position: {self.config.position} (第{self.config.position + 1}个元素|element {self.config.position + 1})")
            self.logger.info(f"分隔符|Delimiter: {self.config.delimiter}")
            self.logger.info("=" * 60)

            # 执行分割|Execute splitting
            self.processor.split_fasta_ids()

            self.logger.info("=" * 60)
            self.logger.info("FASTA ID分割流程完成|FASTA ID splitting pipeline completed")
            self.logger.info(f"结果保存在|Results saved in: {self.config.output_file}")

        except Exception as e:
            self.logger.error(f"分割流程在执行过程中意外终止|Splitting pipeline terminated unexpectedly: {e}")
            sys.exit(1)

def main():
    """主函数|Main function"""
    parser = argparse.ArgumentParser(
        description='FASTA序列ID分割脚本(模块化版本)|FASTA Sequence ID Splitting Script (Modular Version)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例|Examples:
  %(prog)s -i input.fasta -o output.fasta -p 0
        """
    )

    # 必需参数|Required arguments
    parser.add_argument('-i', '--input', required=True,
                       help='输入FASTA文件路径|Input FASTA file path')

    # 可选参数|Optional arguments
    parser.add_argument('-o', '--output', default='output.fasta',
                       help='输出FASTA文件路径|Output FASTA file path')
    parser.add_argument('-p', '--position', type=int, default=0,
                       help='提取位置(0表示第一个元素)|Extract position (0 means first element)')
    parser.add_argument('-d', '--delimiter', default='auto',
                       help='分隔符类型|Delimiter type: "auto"(auto detect), "space", "tab", "both"(space and tab), or any character like "," or "|"')

    # 处理选项|Processing options
    parser.add_argument('--keep-original', action='store_true',
                       help='保留原始文件作为备份|Keep original file as backup')
    parser.add_argument('--no-skip-empty', action='store_true',
                       help='不跳过空的序列名称行|Do not skip empty sequence name lines')
    parser.add_argument('--preserve-comments', action='store_true',
                       help='保留序列名称行中的注释|Preserve comments in sequence name lines')

    args = parser.parse_args()

    # 创建分割器并运行|Create splitter and run
    splitter = FastaIDSplitter(
        input_file=args.input,
        output_file=args.output,
        position=args.position,
        delimiter=args.delimiter,
        keep_original=args.keep_original,
        skip_empty=not args.no_skip_empty,
        preserve_comments=args.preserve_comments
    )

    splitter.run_splitting()

if __name__ == "__main__":
    main()
