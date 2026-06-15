"""
主程序模块|Main Module
"""

import argparse
import sys
from .config import Kmer2VcfConfig
from .utils import Kmer2VcfLogger, sort_file_by_kmer
from .parser import KmerAbundanceParser
from .calculator import KmerToVcfConverter


class Kmer2VcfAnalyzer:
    """Kmer转VCF分析主类|Main Kmer to VCF Analyzer Class"""

    def __init__(self, **kwargs):
        # 初始化配置|Initialize configuration
        self.config = Kmer2VcfConfig(**kwargs)
        self.config.validate()

        # 初始化日志|Initialize logging
        self.logger_manager = Kmer2VcfLogger(
            log_file=self.config.get_log_path(),
            log_level="INFO"
        )
        self.logger = self.logger_manager.get_logger()

        # 初始化解析器和转换器|Initialize parser and converter
        self.parser = KmerAbundanceParser(self.config, self.logger)
        self.converter = KmerToVcfConverter(self.config, self.logger)

    def run_analysis(self):
        """运行完整的分析流程|Run complete analysis pipeline"""
        try:
            self.logger.info("=" * 100)
            self.logger.info(" Kmer丰度转VCF工具|Kmer Abundance to VCF Converter")
            self.logger.info("=" * 100)
            self.logger.info(f"输入文件|Input file: {self.config.input_matrix}")
            self.logger.info(f"输出VCF|Output VCF: {self.config.output_vcf}")
            self.logger.info(f"样本数|Number of samples: {self.config.num_samples}")

            # 判断使用快速模式还是标准模式|Determine whether to use fast mode or standard mode
            use_fast_mode = self.config.fast_mode and not self.config.standard_mode

            if use_fast_mode:
                self.logger.info(f"运行模式|Running mode: 快速模式（单次遍历）|Fast mode (single pass)")
                if self.config.chr_number > 0:
                    self.logger.info(f"染色体数量|Chromosome count: {self.config.chr_number:,}")
                else:
                    self.logger.info(f"染色体长度|Chromosome length: {self.config.chr_length:,}")
                if self.config.min_freq > 0:
                    self.logger.info(f"最小频次过滤|Minimum frequency filter: {self.config.min_freq}")
                self.logger.info("")

                # 使用快速模式|Use fast mode
                from .fast_mode import FastModeConverter
                converter = FastModeConverter(self.config, self.logger)
                kept_lines = converter.convert()

                # 输出汇总信息（快速模式）|Output summary (fast mode)
                self.logger.info("=" * 100)
                self.logger.info(f" 分析完成|Analysis completed successfully")
                self.logger.info("=" * 100)
                self.logger.info(f"输出kmer数量|Kmers output: {format_number(kept_lines)}")
                self.logger.info(f"输出文件|Output file:")
                self.logger.info(f"  VCF: {self.config.output_vcf}")
                self.logger.info(f"  LOG: {self.config.get_log_path()}")
                self.logger.info("=" * 100)

            else:
                self.logger.info(f"运行模式|Running mode: 标准模式（3遍处理）|Standard mode (3-pass processing)")
                self.logger.info(f"最小聚合频次|Minimum aggregated count: {self.config.min_agg_count}")
                self.logger.info("")

                # 使用标准模式（原有逻辑）|Use standard mode (original logic)
                # Pass 1: 转换为0/1格式|Pass 1: Convert to binary format
                total_lines = self.parser.convert_to_binary_format()
                self.logger.info("")

                # Pass 2: 排序（仅标准模式需要）|Pass 2: Sort (only needed for standard mode)
                # 分桶模式下，排序在每个桶内完成|In bucket mode, sorting is done within each bucket
                use_bucket_mode = self.config.should_use_bucket_mode()

                if not use_bucket_mode:
                    self.logger.info("=" * 100)
                    self.logger.info(" Step 2: 排序kmer矩阵|Step 2: Sorting kmer matrix")
                    self.logger.info("=" * 100)

                    binary_file = self.config.get_temp_file_path('binary_matrix.txt')
                    sorted_file = self.config.get_temp_file_path('sorted_binary_matrix.txt')

                    sort_success = sort_file_by_kmer(
                        binary_file,
                        sorted_file,
                        self.config.temp_dir,
                        sort_memory='10G',
                        logger=self.logger
                    )

                    if not sort_success:
                        self.logger.error("排序失败|Sorting failed")
                        self.config.cleanup()
                        sys.exit(1)

                    self.logger.info("")
                else:
                    self.logger.info("=" * 100)
                    self.logger.info(" Step 2: 分桶模式下跳过独立排序步骤|Step 2: Skipping independent sorting in bucket mode")
                    self.logger.info("  (排序将在每个桶内并行完成)|(Sorting will be done in parallel within each bucket)")
                    self.logger.info("=" * 100)
                    self.logger.info("")

                # Pass 3: 统计聚合频次并转换VCF|Pass 3: Calculate aggregated frequencies and convert to VCF
                kept_lines = self.converter.process_sorted_file_and_convert()
                self.logger.info("")

                # 清理临时文件|Clean up temporary files
                self.logger.info("清理临时文件|Cleaning up temporary files...")
                self.config.cleanup()

                # 输出汇总信息|Output summary
                self.logger.info("=" * 100)
                self.logger.info(f" 分析完成|Analysis completed successfully")
                self.logger.info("=" * 100)
                self.logger.info(f"原始行数|Original lines: {format_number(total_lines)}")
                self.logger.info(f"保留行数|Kept lines: {format_number(kept_lines)}")
                if total_lines > 0:
                    percentage = (kept_lines / total_lines) * 100
                    self.logger.info(f"保留比例|Retention rate: {percentage:.2f}%")
                self.logger.info(f"输出文件|Output files:")
                self.logger.info(f"  VCF: {self.config.output_vcf}")
                self.logger.info(f"  TXT: {self.config.get_output_txt_path()}")
                self.logger.info(f"  LOG: {self.config.get_log_path()}")
                self.logger.info("=" * 100)

        except Exception as e:
            self.logger.error(f"分析过程中发生错误|Error during analysis: {e}")
            self.config.cleanup()
            sys.exit(1)


def format_number(num: int) -> str:
    """格式化数字|Format number"""
    if num >= 1_000_000:
        return f"{num / 1_000_000:.2f}M"
    elif num >= 1_000:
        return f"{num / 1_000:.2f}K"
    else:
        return str(num)


def main():
    """主函数|Main function"""
    parser = argparse.ArgumentParser(
        description='Kmer丰度转VCF工具|Kmer Abundance to VCF Converter',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例|Examples:
  # 快速模式（默认，推荐用于大数据集）|Fast mode (default, recommended for large datasets)
  %(prog)s -i abundance_matrix.tsv -o output.vcf

  # 快速模式 + 频次过滤|Fast mode + frequency filter
  %(prog)s -i abundance_matrix.tsv -o output.vcf --min-freq 5

  # 标准模式（3遍处理，用于需要排序的场景）|Standard mode (3-pass processing, for scenarios requiring sorting)
  %(prog)s -i abundance_matrix.tsv -o output.vcf --standard-mode

  # 标准模式 + 聚合频次过滤|Standard mode + aggregated count filter
  %(prog)s -i abundance_matrix.tsv -o output.vcf --standard-mode -m 5

  # 输出压缩VCF|Output compressed VCF
  %(prog)s -i abundance_matrix.tsv -o output.vcf.gz

  # 指定临时文件目录（仅标准模式需要）|Specify temporary directory (only needed for standard mode)
  %(prog)s -i abundance_matrix.tsv -o output.vcf -T /custom/temp/path

模式选择|Mode Selection:
  --fast-mode          快速模式（默认）：1遍遍历，适合10亿+kmer|Fast mode (default): single pass, suitable for 1B+ kmers
  --standard-mode      标准模式：3遍处理+排序，适合需要有序输出的场景|Standard mode: 3-pass + sorting, for ordered output
        """
    )

    # 必需参数|Required arguments
    parser.add_argument('-i', '--input-matrix', required=True,
                       help=' 输入kmer丰度矩阵文件|Input kmer abundance matrix file (TSV format)')
    parser.add_argument('-o', '--output-vcf', required=True,
                       help=' 输出VCF文件路径|Output VCF file path (.vcf or .vcf.gz)')

    # 模式选择参数|Mode selection arguments
    mode_group = parser.add_mutually_exclusive_group()
    mode_group.add_argument('--fast-mode', action='store_true', default=True,
                           help=' 快速模式（单次遍历，默认）|Fast mode (single pass, default)')
    mode_group.add_argument('--standard-mode', action='store_true',
                           help=' 标准模式（3遍处理+排序）|Standard mode (3-pass processing + sorting)')

    # 快速模式参数|Fast mode parameters
    parser.add_argument('--chr-length', type=int, default=100000000,
                       help=' 每条染色体长度（快速模式，默认100M）|Chromosome length for fast mode (default: 100M)')
    parser.add_argument('--chr-number', type=int, default=0,
                       help=' 染色体数量（如果设置则优先使用）|Number of chromosomes (if set, takes priority over --chr-length)')
    parser.add_argument('--min-freq', type=int, default=0,
                       help=' 最小出现频次过滤（快速模式）|Minimum frequency filter for fast mode (default: 0=no filter)')
    parser.add_argument('--kmer-length', type=int, default=51,
                       help=' Kmer长度，用于VCF INFO字段|Kmer length for VCF INFO field (default: 51)')
    parser.add_argument('--no-header', action='store_true',
                       help=' 输入文件没有header行（第一行就是样本名）|Input file has no header line (first line is sample names)')

    # 标准模式参数|Standard mode parameters
    parser.add_argument('-m', '--min-agg-count', type=int, default=3,
                       help=' 最小聚合频次阈值（标准模式）|Minimum aggregated count threshold for standard mode (default: 3)')

    # 通用参数|Common parameters
    parser.add_argument('-t', '--threads', type=int, default=12,
                       help=' 线程数（标准模式）|Number of threads for standard mode (default: 12)')
    parser.add_argument('-T', '--temp-dir', default=None,
                       help=' 临时文件目录（标准模式）|Temporary directory for standard mode (default: ./temp)')

    args = parser.parse_args()

    # 创建分析器并运行|Create analyzer and run
    analyzer = Kmer2VcfAnalyzer(
        input_matrix=args.input_matrix,
        output_vcf=args.output_vcf,
        fast_mode=args.fast_mode,
        standard_mode=args.standard_mode,
        chr_length=args.chr_length,
        chr_number=args.chr_number,
        min_freq=args.min_freq,
        kmer_length=args.kmer_length,
        no_header=args.no_header,
        min_agg_count=args.min_agg_count,
        threads=args.threads,
        temp_dir=args.temp_dir
    )

    analyzer.run_analysis()


if __name__ == "__main__":
    main()
