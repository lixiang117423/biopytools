"""
FASTQ GC过滤主程序模块|FASTQ GC Filter Main Module
"""

import argparse
import sys
from .config import FastqGcFilterConfig
from .utils import FastqGcFilterLogger
from .filter import FastqFilter


class FastqGcFilter:
    """FASTQ GC过滤主类|Main FASTQ GC Filter Class"""

    def __init__(self, **kwargs):
        """
        初始化过滤器|Initialize filter

        Args:
            **kwargs: 配置参数|Configuration parameters
        """
        # 初始化配置|Initialize configuration
        self.config = FastqGcFilterConfig(**kwargs)
        self.config.validate()

        # 初始化日志|Initialize logging
        self.logger_manager = FastqGcFilterLogger()
        self.logger = self.logger_manager.get_logger()

        # 初始化过滤器|Initialize filter
        self.filter = FastqFilter(self.config, self.logger)

    def run(self):
        """运行过滤|Run filtering"""
        try:
            total_reads, passed_reads = self.filter.filter_fastq()

            # 输出统计信息|Output statistics
            self.logger.info("-" * 50)
            self.logger.info(f"处理完成|Processing completed")
            self.logger.info(f"总reads数|Total reads: {self.filter.format_number(total_reads)}")
            self.logger.info(f"通过筛选|Passed filters: {self.filter.format_number(passed_reads)}")
            self.logger.info(f"过滤掉|Filtered out: {self.filter.format_number(total_reads - passed_reads)}")
            if total_reads > 0:
                pass_rate = (passed_reads / total_reads) * 100
                self.logger.info(f"通过率|Pass rate: {pass_rate:.2f}%")
            self.logger.info(f"输出文件|Output file: {self.config.output_file}")

            return True

        except Exception as e:
            self.logger.error(f"程序执行出错|Program execution error: {str(e)}")
            return False


def main():
    """主函数|Main function"""
    parser = argparse.ArgumentParser(
        description='筛选FASTQ文件中GC含量和序列长度在指定范围内的reads|Filter FASTQ reads by GC content and sequence length',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # 必需参数|Required arguments
    parser.add_argument('-i', '--input', required=True,
                       help='输入FASTQ文件路径|Input FASTQ file path (支持.gz压缩|supports .gz compression)')
    parser.add_argument('-o', '--output', required=True,
                       help='输出FASTQ文件路径|Output FASTQ file path (支持.gz压缩|supports .gz compression)')

    # GC含量过滤参数|GC content filtering parameters
    parser.add_argument('--min-gc', type=float, default=25.0,
                       help='最小GC含量百分比|Minimum GC content percentage')
    parser.add_argument('--max-gc', type=float, default=100.0,
                       help='最大GC含量百分比|Maximum GC content percentage')

    # 序列长度过滤参数|Sequence length filtering parameters
    parser.add_argument('--min-length', type=int, default=50,
                       help='最短序列长度|Minimum sequence length')
    parser.add_argument('--max-length', type=int, default=None,
                       help='最长序列长度|Maximum sequence length')

    args = parser.parse_args()

    # 创建过滤器并运行|Create filter and run
    filter_tool = FastqGcFilter(
        input_file=args.input,
        output_file=args.output,
        min_gc=args.min_gc,
        max_gc=args.max_gc,
        min_length=args.min_length,
        max_length=args.max_length
    )

    success = filter_tool.run()
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
