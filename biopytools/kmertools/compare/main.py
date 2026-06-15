"""
主程序模块|Main Module
"""

import sys
from .config import KmerCompareConfig
from .utils import KmerCompareLogger, format_number
from .calculator import KmerCompareCalculator


class KmerMatrixComparator:
    """Kmer矩阵比较主类|Main Kmer Matrix Comparator Class"""

    def __init__(self, **kwargs):
        """
        初始化比较器|Initialize comparator

        Args:
            **kwargs: 配置参数|Configuration parameters
        """
        # 初始化配置|Initialize configuration
        self.config = KmerCompareConfig(**kwargs)
        self.config.validate()

        # 初始化日志|Initialize logging
        self.logger_manager = KmerCompareLogger(
            log_file=self.config.log_file,
            log_level="INFO"
        )
        self.logger = self.logger_manager.get_logger()

        # 初始化计算器|Initialize calculator
        self.calculator = KmerCompareCalculator(self.config, self.logger)

    def run_comparison(self):
        """运行完整的比较流程|Run complete comparison pipeline"""
        try:
            self.logger.info("=" * 100)
            self.logger.info(" Kmer矩阵比较工具|Kmer Matrix Comparison Tool")
            self.logger.info("=" * 100)
            self.logger.info(f"文件1|File1: {self.config.file1}")
            self.logger.info(f"文件2|File2: {self.config.file2}")
            self.logger.info(f"输出前缀|Output prefix: {self.config.output_prefix}")
            self.logger.info(f"窗口大小|Window size: {format_number(self.config.window_size)}")
            self.logger.info("")

            # Step 1: 提取kmer序列集合|Step 1: Extract kmer sequence sets
            self.logger.info("=" * 100)
            self.logger.info(" 步骤 1/4: 提取kmer序列集合|Step 1/4: Extract kmer sequence sets")
            self.logger.info("=" * 100)
            seqs1 = self.calculator.extract_sequences(self.config.file1)
            seqs2 = self.calculator.extract_sequences(self.config.file2)
            self.logger.info("")

            # Step 2: 找出特有序列|Step 2: Find unique sequences
            self.logger.info("=" * 100)
            self.logger.info(" 步骤 2/4: 找出特有序列|Step 2/4: Find unique sequences")
            self.logger.info("=" * 100)
            unique1, unique2 = self.calculator.find_unique_sequences(seqs1, seqs2)
            self.logger.info("")

            # Step 3: 处理文件1|Step 3: Process file 1
            self.logger.info("=" * 100)
            self.logger.info(" 步骤 3/4: 处理文件1并统计|Step 3/4: Process file 1 and calculate statistics")
            self.logger.info("=" * 100)
            sample_names1 = self.calculator.process_file_with_stats(
                self.config.file1, unique1, self.config.output_file1,
                self.config.file1.replace('_matrix.txt', '').replace('.txt', '').split('/')[-1]
            )
            self.logger.info("")

            # Step 4: 处理文件2|Step 4: Process file 2
            self.logger.info("=" * 100)
            self.logger.info(" 步骤 4/4: 处理文件2并统计|Step 4/4: Process file 2 and calculate statistics")
            self.logger.info("=" * 100)
            sample_names2 = self.calculator.process_file_with_stats(
                self.config.file2, unique2, self.config.output_file2,
                self.config.file2.replace('_matrix.txt', '').replace('.txt', '').split('/')[-1]
            )
            self.logger.info("")

            # 输出汇总信息|Output summary information
            self.logger.info("=" * 100)
            self.logger.info(f" 分析完成|Analysis Completed Successfully")
            self.logger.info("=" * 100)
            self.logger.info(f"文件1特有kmer|Unique kmers in file1: {format_number(len(unique1))}")
            self.logger.info(f"文件2特有kmer|Unique kmers in file2: {format_number(len(unique2))}")
            self.logger.info(f"文件1统计结果|File1 statistics: {self.config.output_file1}")
            self.logger.info(f"文件2统计结果|File2 statistics: {self.config.output_file2}")
            self.logger.info(f"文件1样品数|File1 samples: {len(sample_names1)}")
            self.logger.info(f"文件2样品数|File2 samples: {len(sample_names2)}")
            self.logger.info(f"日志文件|Log file: {self.config.log_file}")
            self.logger.info("=" * 100)

        except Exception as e:
            self.logger.error(f"分析过程中发生错误|Error during analysis: {e}")
            import traceback
            self.logger.error(f"错误详情|Error details:\n{traceback.format_exc()}")
            sys.exit(1)


def parse_arguments():
    """解析命令行参数|Parse command line arguments"""
    import argparse

    parser = argparse.ArgumentParser(
        description='Kmer矩阵比较工具|Compare two kmer matrix files and find unique kmers',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例|Examples:
  # 基本用法|Basic usage
  %(prog)s -f1 matrix1.txt -f2 matrix2.txt -o comparison

  # 指定窗口大小|Specify window size
  %(prog)s -f1 matrix1.txt -f2 matrix2.txt -o comparison -w 50000

说明|Notes:
  - 使用流式处理，不会将整个文件载入内存
  - Uses streaming processing, does not load entire files into memory
  - 输出两个文件：{prefix}_file1_stats.txt 和 {prefix}_file2_stats.txt
  - Outputs two files: {prefix}_file1_stats.txt and {prefix}_file2_stats.txt
        """
    )

    # 必需参数|Required arguments
    required = parser.add_argument_group('必需参数|Required arguments')
    required.add_argument('-f1', '--file1', required=True,
                        help='第一个kmer矩阵文件|First kmer matrix file')
    required.add_argument('-f2', '--file2', required=True,
                        help='第二个kmer矩阵文件|Second kmer matrix file')
    required.add_argument('-o', '--output-prefix', required=True,
                        help='输出文件前缀|Output file prefix')

    # 可选参数|Optional arguments
    optional = parser.add_argument_group('可选参数|Optional arguments')
    optional.add_argument('-w', '--window-size', type=int, default=100000,
                        help='窗口大小（行数）|Window size in lines (default: %(default)s)')

    return parser.parse_args()


def main():
    """主函数|Main function"""
    args = parse_arguments()

    try:
        # 创建比较器并运行|Create comparator and run
        comparator = KmerMatrixComparator(
            file1=args.file1,
            file2=args.file2,
            output_prefix=args.output_prefix,
            window_size=args.window_size
        )

        comparator.run_comparison()

        sys.exit(0)

    except KeyboardInterrupt:
        print("\n分析被用户中断|Analysis interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"错误|Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
