"""
VCF基因型统计主程序模块|VCF Genotype Statistics Main Module
"""

import argparse
import sys
from .config import VCFStatsConfig
from .utils import VCFStatsLogger
from .data_processing import VCFProcessor, StatisticsCalculator
from .results import ResultsExporter, SummaryGenerator


class VCFStatsAnalyzer:
    """VCF基因型统计主类|Main VCF Genotype Statistics Class"""

    def __init__(self, **kwargs):
        # 初始化配置|Initialize configuration
        self.config = VCFStatsConfig(**kwargs)
        self.config.validate()

        # 初始化日志|Initialize logging
        self.logger_manager = VCFStatsLogger(
            self.config.output_path,
            verbose=self.config.verbose,
            quiet=self.config.quiet
        )
        self.logger = self.logger_manager.get_logger()

        # 初始化各个处理器|Initialize processors
        self.vcf_processor = VCFProcessor(self.config, self.logger)
        self.stats_calculator = StatisticsCalculator(self.config, self.logger)
        self.results_exporter = ResultsExporter(self.config, self.logger)
        self.summary_generator = SummaryGenerator(self.config, self.logger)

    def run_analysis(self):
        """运行完整的VCF基因型统计分析流程|Run complete VCF genotype statistics analysis pipeline"""
        try:
            self.logger.info("=" * 80)
            self.logger.info("开始VCF基因型统计分析|Starting VCF Genotype Statistics Analysis")
            self.logger.info("=" * 80)
            self.logger.info(f"VCF文件|VCF file: {self.config.vcf_file}")
            self.logger.info(f"输出目录|Output directory: {self.config.output_dir}")

            if self.config.min_depth > 0:
                self.logger.info(f"最小深度过滤|Minimum depth filter: {self.config.min_depth}")
            if self.config.min_qual > 0:
                self.logger.info(f"最小质量过滤|Minimum quality filter: {self.config.min_qual}")
            if self.config.exclude_missing:
                self.logger.info("将排除缺失基因型|Will exclude missing genotypes")

            # 步骤1: 处理VCF文件|Step 1: Process VCF file
            self.logger.info("步骤1/3: 处理VCF文件|Step 1/3: Processing VCF file")
            if not self.vcf_processor.process_vcf():
                self.logger.error("VCF文件处理失败|VCF file processing failed")
                sys.exit(1)

            # 步骤2: 计算统计结果|Step 2: Calculate statistics
            self.logger.info("步骤2/3: 计算统计结果|Step 2/3: Calculating statistics")
            stats_results = self.stats_calculator.calculate_rates(self.vcf_processor.sample_stats)

            if not stats_results:
                self.logger.error("统计计算失败|Statistics calculation failed")
                sys.exit(1)

            self.logger.info(f"成功计算了 {len(stats_results)} 个样本的统计结果|Successfully calculated statistics for {len(stats_results)} samples")

            # 步骤3: 导出结果|Step 3: Export results
            self.logger.info("步骤3/3: 导出结果|Step 3/3: Exporting results")
            self.results_exporter.export_summary_statistics(stats_results, self.vcf_processor.detailed_stats)
            self.results_exporter.export_detailed_statistics(self.vcf_processor.detailed_stats)
            self.results_exporter.export_per_sample_files(stats_results, self.vcf_processor.detailed_stats)

            # 生成分析总结|Generate analysis summary
            self.summary_generator.generate_analysis_summary(stats_results)

            self.logger.info("=" * 20 + " VCF基因型统计分析完成|VCF Genotype Statistics Analysis Completed " + "=" * 20)
            self.logger.info(f"结果文件已保存到|Results saved to: {self.config.output_dir}")

        except Exception as e:
            self.logger.error(f"分析流程在执行过程中意外终止|Analysis pipeline terminated unexpectedly: {e}")
            sys.exit(1)


def main():
    """主函数|Main function"""
    parser = argparse.ArgumentParser(
        description="VCF基因型统计分析脚本|VCF Genotype Statistics Analysis Script",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例|Examples:
  %(prog)s -i variants.vcf -o vcf_stats_output
        """
    )

    # 必需参数|Required arguments
    required = parser.add_argument_group('必需参数|required arguments')
    required.add_argument("-i", "--input", required=True,
                         help="输入VCF文件路径|Input VCF file path")

    # 输出配置|Output configuration
    output = parser.add_argument_group('输出配置|Output configuration')
    output.add_argument("-o", "--output", default="vcf_stats_output",
                       help="输出目录|Output directory")

    # 过滤参数|Filtering parameters
    filtering = parser.add_argument_group('过滤参数|Filtering parameters')
    filtering.add_argument("-d", "--min-depth", type=int, default=0,
                          help="最小深度过滤阈值|Minimum depth filter threshold")
    filtering.add_argument("-q", "--min-qual", type=float, default=0.0,
                          help="最小质量分数过滤阈值|Minimum quality score filter threshold")
    filtering.add_argument("-e", "--exclude-missing", action="store_true",
                          help="排除缺失基因型|Exclude missing genotypes")

    # 输出控制|Output control
    output_ctrl = parser.add_argument_group('输出控制|Output control')
    output_ctrl.add_argument("-D", "--no-detailed", action="store_true",
                            help="不输出详细统计结果|Do not output detailed statistics")
    output_ctrl.add_argument("-S", "--no-summary", action="store_true",
                            help="不输出汇总统计结果|Do not output summary statistics")

    # 日志选项|Logging options
    logging = parser.add_argument_group('日志选项|Logging options')
    logging.add_argument('--verbose', '-v', action='count', default=0,
                        help='增加输出详细程度|Increase output verbosity')
    logging.add_argument('--quiet', action='store_true',
                        help='静默模式|Quiet mode')
    logging.add_argument('--log-file', type=str,
                        help='日志文件路径|Log file path')
    logging.add_argument('--log-level', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
                        default='INFO',
                        help='日志级别|Log level')

    # 高级选项|Advanced options
    advanced = parser.add_argument_group('高级选项|Advanced options')
    advanced.add_argument('--dry-run', action='store_true',
                         help='试运行模式|Dry run mode')

    args = parser.parse_args()

    # 创建分析器并运行|Create analyzer and run
    analyzer = VCFStatsAnalyzer(
        vcf_file=args.input,
        output_dir=args.output,
        min_depth=args.min_depth,
        min_qual=args.min_qual,
        exclude_missing=args.exclude_missing,
        output_detailed=not args.no_detailed,
        output_summary=not args.no_summary,
        verbose=(args.verbose > 0),
        quiet=args.quiet,
        log_file=args.log_file,
        log_level=args.log_level,
        dry_run=args.dry_run
    )

    analyzer.run_analysis()


if __name__ == "__main__":
    main()
