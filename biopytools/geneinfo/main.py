"""
GFF3工具主程序模块|GFF3 Tools Main Module
"""

import argparse
import sys
from .config import GFFConfig
from .utils import GFFLogger
from .gff_extractor import GFFExtractor
from .results import SummaryGenerator

class GFFAnalyzer:
    """GFF3分析主类|Main GFF3 Analyzer Class"""

    def __init__(self, **kwargs):
        # 初始化配置|Initialize configuration
        self.config = GFFConfig(**kwargs)
        self.config.validate()

        # 初始化日志|Initialize logging
        self.logger_manager = GFFLogger(
            self.config.output_file,
            verbose=self.config.verbose,
            quiet=self.config.quiet
        )
        self.logger = self.logger_manager.get_logger()

        # 初始化处理器|Initialize processors
        self.gff_extractor = GFFExtractor(self.config, self.logger)
        self.summary_generator = SummaryGenerator(self.config, self.logger)

    def run_extraction(self):
        """
        运行完整的提取流程|Run complete extraction pipeline
        """
        try:
            self.logger.info("开始GFF3基因转录本提取分析|Starting GFF3 gene transcript extraction analysis")
            self.logger.info(f"输入路径|Input path: {self.config.gff3_file}")
            self.logger.info(f"输出文件|Output file: {self.config.output_file}")

            # 提取基因转录本信息|Extract gene transcript information
            if len(self.config.gff3_files) > 1:
                self.logger.info(f"批量模式|Batch mode: 共 {len(self.config.gff3_files)} 个GFF3文件")
                transcript_data = self.gff_extractor.extract_all_files()
            else:
                transcript_data = self.gff_extractor.extract_gene_transcript_info()

            # 写入结果|Write results
            self.gff_extractor.write_results(transcript_data)

            # 生成总结报告|Generate summary report
            self.summary_generator.generate_summary_report(transcript_data)

            # 打印摘要|Print summary
            self.gff_extractor.print_summary(transcript_data)

            self.logger.info("提取完成|Extraction completed successfully")

        except Exception as e:
            self.logger.error(f"提取失败|Extraction failed: {e}")
            sys.exit(1)

def main():
    """主函数|Main function"""
    parser = argparse.ArgumentParser(
        description="从GFF3文件中为每个转录本提取整合的基因和转录本信息|Extract integrated gene and transcript information for each transcript from GFF3 files",
        epilog="示例|Example: %(prog)s -i input.gff3 -o gene_transcript_info.tsv",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    # 必需参数|Required arguments
    required = parser.add_argument_group('必需参数|Required arguments')
    required.add_argument(
        '--input', '-i',
        required=True,
        help="输入GFF3文件或目录路径|Input GFF3 file or directory path"
    )
    required.add_argument(
        '--output', '-o',
        required=True,
        help="输出的TSV文件路径|Output TSV file path"
    )

    # 可选参数|Optional arguments
    optional = parser.add_argument_group('可选参数|Optional arguments')
    optional.add_argument(
        '--gene-type',
        default='gene',
        help="基因特征类型|Gene feature type"
    )
    optional.add_argument(
        '--transcript-types',
        nargs='+',
        default=['mRNA', 'transcript'],
        help="转录本特征类型列表|Transcript feature types list"
    )

    # 日志选项|Logging options
    logging = parser.add_argument_group('日志选项|Logging options')
    logging.add_argument('-v', '--verbose',
                        action='count',
                        default=0,
                        help='增加输出详细程度|Increase output verbosity')
    logging.add_argument('--quiet',
                        action='store_true',
                        help='静默模式，仅输出错误信息|Quiet mode, only output errors')
    logging.add_argument('--log-file',
                        type=str,
                        help='日志文件路径|Log file path')
    logging.add_argument('--log-level',
                        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
                        default='INFO',
                        help='日志级别|Log level')

    # 高级选项|Advanced options
    advanced = parser.add_argument_group('高级选项|Advanced options')
    advanced.add_argument('--dry-run',
                         action='store_true',
                         help='试运行模式，不实际执行|Dry run mode, no actual execution')
    advanced.add_argument('-t', '--threads',
                         type=int,
                         default=1,
                         help='线程数|Number of threads')

    args = parser.parse_args()

    # 创建分析器并运行|Create analyzer and run
    try:
        analyzer = GFFAnalyzer(
            gff3_file=args.input,
            output_file=args.output,
            gene_type=args.gene_type,
            transcript_types=set(args.transcript_types),
            verbose=(args.verbose > 0),
            quiet=args.quiet,
            log_file=args.log_file,
            log_level=args.log_level,
            dry_run=args.dry_run,
            threads=args.threads
        )

        analyzer.run_extraction()

    except ValueError as e:
        print(f"参数错误|Parameter error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"错误|Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
