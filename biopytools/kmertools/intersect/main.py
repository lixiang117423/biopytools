"""
主程序模块|Main Module
"""

import argparse
import sys
import os
from .config import KmerIntersectConfig
from .utils import KmerIntersectLogger
from .calculator import KmerIntersectCalculator


class KmerIntersectAnalyzer:
    """Kmer交集分析主类|Main Kmer Intersection Analyzer Class"""

    def __init__(self, **kwargs):
        """
        初始化分析器|Initialize analyzer

        Args:
            **kwargs: 配置参数|Configuration parameters
        """
        # 初始化配置|Initialize configuration
        self.config = KmerIntersectConfig(**kwargs)
        self.config.validate()

        # 初始化日志|Initialize logging
        log_file = self.config.output_file + '.log'
        self.logger_manager = KmerIntersectLogger(
            log_file=log_file,
            log_level="INFO"
        )
        self.logger = self.logger_manager.get_logger()

        # 初始化计算器|Initialize calculator
        self.calculator = KmerIntersectCalculator(self.config, self.logger)

    def run_analysis(self):
        """运行完整的分析流程|Run complete analysis pipeline"""
        try:
            self.logger.info("=" * 100)
            self.logger.info(" Kmer交集提取工具|Kmer Intersection Extraction Tool")
            self.logger.info("=" * 100)
            self.logger.info(f"kmer矩阵文件|Kmer matrix file: {self.config.kmer_matrix}")
            self.logger.info(f"目标kmer fasta|Target kmer fasta: {self.config.kmer_fasta}")
            self.logger.info(f"输出文件|Output file: {self.config.output_file}")
            self.logger.info(f"线程数|Threads: {self.config.threads}")
            self.logger.info(f"使用反向互补查询|Use reverse complement: {self.config.use_reverse_complement}")
            self.logger.info("")

            # Step 1: 加载目标kmer|Step 1: Load target kmers
            target_kmers, kmer_id_list, kmer_id_to_seq = self.calculator.load_target_kmers()
            self.logger.info("")

            # Step 2: 提取丰度信息|Step 2: Extract abundance information
            results, found_count, sample_names = self.calculator.extract_kmer_abundance(
                target_kmers=target_kmers,
                sample_names=None  # 从文件header读取|Read from file header
            )
            self.logger.info("")

            # Step 3: 写入结果|Step 3: Write results
            self.calculator.write_results(
                results=results,
                sample_names=sample_names,
                target_kmers=target_kmers,
                kmer_id_list=kmer_id_list,
                kmer_id_to_seq=kmer_id_to_seq
            )
            self.logger.info("")

            # 输出汇总信息|Output summary information
            self.logger.info("=" * 100)
            self.logger.info(f" 分析完成|Analysis Completed Successfully")
            self.logger.info("=" * 100)
            self.logger.info(f"目标kmer总数|Total target kmers: {len(target_kmers)}")
            self.logger.info(f"找到kmer数|Kmers found: {found_count}")
            self.logger.info(f"未找到kmer数|Kmers not found: {len(target_kmers) - found_count}")
            if len(target_kmers) > 0:
                found_rate = (found_count / len(target_kmers)) * 100
                self.logger.info(f"找到率|Found rate: {found_rate:.2f}%")
            self.logger.info(f"输出文件|Output file: {self.config.output_file}")
            self.logger.info(f"日志文件|Log file: {self.config.output_file}.log")
            self.logger.info("=" * 100)

        except Exception as e:
            self.logger.error(f"分析过程中发生错误|Error during analysis: {e}")
            import traceback
            self.logger.error(f"错误详情|Error details:\n{traceback.format_exc()}")
            sys.exit(1)


def parse_arguments():
    """解析命令行参数|Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='Kmer交集提取工具|Extract kmer abundance from matrix based on target kmer list',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例|Examples:
  # 基本用法|Basic usage
  %(prog)s -m kmer_matrix.txt -k target_kmers.fa -o results.txt

  # 使用反向互补查询|Use reverse complement query
  %(prog)s -m kmer_matrix.txt -k target_kmers.fa -o results.txt -r

  # 输出CSV格式|Output CSV format
  %(prog)s -m kmer_matrix.txt -k target_kmers.fa -o results.csv -f csv

  # 不保留未找到的kmer|Do not keep kmers not found
  %(prog)s -m kmer_matrix.txt -k target_kmers.fa -o results.txt --no-keep-not-found
        """
    )

    # 必需参数|Required arguments
    required = parser.add_argument_group('必需参数|Required arguments')
    required.add_argument('-m', '--kmer-matrix', required=True,
                        help='kmer矩阵文件|Kmer matrix file (TSV format)')
    required.add_argument('-k', '--kmer-fasta', required=True,
                        help='目标kmer的fasta文件|Target kmer fasta file')
    required.add_argument('-o', '--output-file', required=True,
                        help='输出文件路径|Output file path')

    # 可选参数|Optional arguments
    optional = parser.add_argument_group('可选参数|Optional arguments')
    optional.add_argument('-t', '--threads', type=int, default=12,
                        help='线程数|Number of threads (default: %(default)s)')
    optional.add_argument('-r', '--use-reverse-complement', action='store_true', default=True,
                        help='使用反向互补查询（默认启用）|Use reverse complement query (default: enabled)')
    optional.add_argument('--no-reverse-complement', dest='use_reverse_complement',
                        action='store_false',
                        help='不使用反向互补查询|Do not use reverse complement query')
    optional.add_argument('--keep-not-found', action='store_true', default=True,
                        help='保留未找到的kmer（默认保留）|Keep kmers not found (default: enabled)')
    optional.add_argument('--no-keep-not-found', dest='keep_not_found',
                        action='store_false',
                        help='不保留未找到的kmer|Do not keep kmers not found')
    optional.add_argument('-f', '--output-format', choices=['tsv', 'csv'], default='tsv',
                        help='输出格式|Output format (default: %(default)s)')
    optional.add_argument('-w', '--window-size', type=int, default=100000,
                        help='窗口大小：每N个kmer统计一次总数|Window size for statistics per N kmers (default: %(default)s)')

    return parser.parse_args()


def main():
    """主函数|Main function"""
    args = parse_arguments()

    try:
        # 创建分析器并运行|Create analyzer and run
        analyzer = KmerIntersectAnalyzer(
            kmer_matrix=args.kmer_matrix,
            kmer_fasta=args.kmer_fasta,
            output_file=args.output_file,
            threads=args.threads,
            use_reverse_complement=args.use_reverse_complement,
            keep_not_found=args.keep_not_found,
            output_format=args.output_format,
            window_size=args.window_size
        )

        analyzer.run_analysis()

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
