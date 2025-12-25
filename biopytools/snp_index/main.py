"""
SNP Index主程序模块 | SNP Index Main Program Module
"""

import argparse
import logging
import sys
import time
from typing import Optional, List, Dict

from .config import SNPIndexConfig
from .utils import setup_logger
from .calculator import SNPIndexCalculator
from .analyzer import SNPIndexAnalyzer
from .visualizer import SNPIndexVisualizer


class SNPIndexProcessor:
    """SNP Index处理器主类 | SNP Index Processor Main Class"""

    def __init__(self, config: SNPIndexConfig):
        """
        初始化处理器 | Initialize processor

        Args:
            config: 配置对象 | Configuration object
        """
        self.config = config
        self.logger = setup_logger(
            self.__class__.__name__,
            config.log_file,
            config.get_log_level()
        )

        # 验证配置 | Validate configuration
        self.config.validate()

        # 初始化组件 | Initialize components
        self.calculator = None
        self.analyzer = None
        self.visualizer = None

    def calculate_snp_index(self) -> bool:
        """
        计算SNP index | Calculate SNP index

        Returns:
            bool: 计算是否成功 | Whether calculation succeeded
        """
        if not self.config.input_vcf:
            self.logger.error("未指定输入VCF文件 | Input VCF file not specified")
            return False

        self.logger.info("=" * 60)
        self.logger.info("开始SNP index计算 | Starting SNP index calculation")
        self.logger.info("=" * 60)

        # 设置默认输出文件路径 | Set default output file path
        if not self.config.output_file:
            self.config.output_file = self.config.get_output_path(f"{self.config.prefix}_results.tsv")

        # 创建计算器 | Create calculator
        self.calculator = SNPIndexCalculator(self.config)

        # 执行计算 | Perform calculation
        return self.calculator.calculate()

    def analyze_results(self, result_file: Optional[str] = None) -> bool:
        """
        分析SNP index结果 | Analyze SNP index results

        Args:
            result_file: 结果文件路径 | Result file path

        Returns:
            bool: 分析是否成功 | Whether analysis succeeded
        """
        if not result_file:
            result_file = self.config.output_file

        if not result_file:
            self.logger.error("未指定结果文件 | Result file not specified")
            return False

        self.logger.info("=" * 60)
        self.logger.info("开始结果分析 | Starting result analysis")
        self.logger.info("=" * 60)

        # 创建分析器 | Create analyzer
        self.analyzer = SNPIndexAnalyzer(result_file, self.config)

        # 执行分析 | Perform analysis
        success = self.analyzer.analyze()

        if success:
            # 导出结果 | Export results
            output_prefix = self.config.get_output_path(self.config.prefix)
            self.analyzer.export_results(output_prefix)

        return success

    def _save_sliding_window_results(self) -> bool:
        """
        计算和保存滑动窗口结果 | Calculate and save sliding window results

        Returns:
            bool: 保存是否成功 | Whether save succeeded
        """
        try:
            self.logger.info("=" * 60)
            self.logger.info("计算滑动窗口分析 | Calculating sliding window analysis")
            self.logger.info("=" * 60)

            # 检查是否有数据 | Check if data is available
            if not self.analyzer or not hasattr(self.analyzer, 'data') or not self.analyzer.data:
                self.logger.warning("没有可用的数据进行滑动窗口分析 | No data available for sliding window analysis")
                return False

            # 导入滑动窗口分析器 | Import sliding window analyzer
            from .sliding_window import SlidingWindowAnalyzer

            # 创建滑动窗口分析器 | Create sliding window analyzer
            analyzer = SlidingWindowAnalyzer(self.analyzer.data, self.config)

            # 创建滑动窗口 | Create sliding windows
            windows = analyzer.create_sliding_windows()

            if not windows:
                self.logger.warning("没有创建任何滑动窗口 | No sliding windows created")
                return False

            # 保存窗口数据 | Save window data
            window_output_file = self.config.get_output_path(f"{self.config.prefix}_sliding_windows.tsv")
            if not analyzer.export_windows(window_output_file):
                return False

            # 计算置信区间 | Calculate confidence intervals
            try:
                confidence_intervals = analyzer.calculate_confidence_intervals(self.config.confidence_level)
                if confidence_intervals:
                    ci_output_file = self.config.get_output_path(f"{self.config.prefix}_confidence_intervals.txt")
                    self._save_confidence_intervals(confidence_intervals, ci_output_file)
            except Exception as e:
                self.logger.warning(f"置信区间计算失败 | Confidence interval calculation failed: {str(e)}")

            # 识别候选区域 | Identify candidate regions
            candidate_regions = analyzer.identify_candidate_regions()
            if candidate_regions:
                regions_output_file = self.config.get_output_path(f"{self.config.prefix}_candidate_regions.tsv")
                self._save_candidate_regions(candidate_regions, regions_output_file)

            # 显示摘要统计 | Show summary statistics
            stats = analyzer.get_summary_statistics()
            self._log_sliding_window_stats(stats)

            self.logger.info("滑动窗口分析完成 | Sliding window analysis completed")
            return True

        except Exception as e:
            self.logger.error(f"滑动窗口分析时出错 | Error in sliding window analysis: {str(e)}", exc_info=True)
            return False

    def _save_confidence_intervals(self, confidence_intervals: dict, output_file: str) -> None:
        """
        保存置信区间结果 | Save confidence interval results

        Args:
            confidence_intervals: 置信区间数据 | Confidence interval data
            output_file: 输出文件 | Output file
        """
        try:
            with open(output_file, 'w', encoding='utf-8') as f:
                f.write("滑动窗口置信区间分析结果 | Sliding Window Confidence Interval Analysis Results\n")
                f.write("=" * 80 + "\n\n")
                f.write(f"置信水平 | Confidence Level: {confidence_intervals.get('confidence_level', 'N/A')}\n")
                f.write(f"样本数 | Sample Count: {len(confidence_intervals.get('all_deltas', []))}\n")
                f.write("\n")

                f.write("统计量 | Statistics:\n")
                f.write(f"  平均值 | Mean: {confidence_intervals.get('mean', 'N/A'):.4f}\n")
                f.write(f"  标准差 | Standard Deviation: {confidence_intervals.get('std', 'N/A'):.4f}\n")
                f.write(f"  误差范围 | Margin of Error: {confidence_intervals.get('margin_error', 'N/A'):.4f}\n")
                f.write("\n")

                f.write("置信区间 | Confidence Intervals:\n")
                f.write(f"  下限 | Lower Bound: {confidence_intervals.get('ci_lower', 'N/A'):.4f}\n")
                f.write(f"  上限 | Upper Bound: {confidence_intervals.get('ci_upper', 'N/A'):.4f}\n")
                f.write("\n")

                f.write("基于百分位的置信区间 | Percentile-based Confidence Intervals:\n")
                f.write(f"  下限 | Lower Bound: {confidence_intervals.get('ci_lower_percentile', 'N/A'):.4f}\n")
                f.write(f"  上限 | Upper Bound: {confidence_intervals.get('ci_upper_percentile', 'N/A'):.4f}\n")

        except Exception as e:
            self.logger.error(f"保存置信区间结果时出错 | Error saving confidence intervals: {str(e)}")

    def _save_candidate_regions(self, candidate_regions: List[Dict], output_file: str) -> None:
        """
        保存候选区域结果 | Save candidate region results

        Args:
            candidate_regions: 候选区域列表 | Candidate region list
            output_file: 输出文件 | Output file
        """
        try:
            import csv

            with open(output_file, 'w', newline='', encoding='utf-8') as f:
                writer = csv.writer(f, delimiter='\t')

                # 写入表头 | Write header
                writer.writerow([
                    'Chromosome', 'Start', 'End', 'Length', 'Window_Count',
                    'Mean_Delta_SNP_Index', 'Max_Delta_SNP_Index'
                ])

                # 写入数据 | Write data
                for region in candidate_regions:
                    length = region['end'] - region['start']
                    writer.writerow([
                        region['chromosome'],
                        region['start'],
                        region['end'],
                        length,
                        region['window_count'],
                        f"{region['mean_delta']:.4f}",
                        f"{region['max_delta']:.4f}"
                    ])

        except Exception as e:
            self.logger.error(f"保存候选区域结果时出错 | Error saving candidate regions: {str(e)}")

    def _log_sliding_window_stats(self, stats: Dict) -> None:
        """
        记录滑动窗口统计信息 | Log sliding window statistics

        Args:
            stats: 统计数据 | Statistics data
        """
        self.logger.info(f"总窗口数 | Total windows: {stats.get('total_windows', 0):,}")
        self.logger.info(f"染色体数 | Chromosome count: {stats.get('chromosome_count', 0)}")
        self.logger.info(f"平均ΔSNP index | Mean ΔSNP index: {stats.get('mean_delta_mean', 0):.4f}")
        self.logger.info(f"平均每窗口SNP数 | Average SNPs per window: {stats.get('mean_snp_count', 0):.1f}")

    def create_visualizations(self, result_file: Optional[str] = None) -> bool:
        """
        创建可视化图表 | Create visualization plots

        Args:
            result_file: 结果文件路径 | Result file path

        Returns:
            bool: 创建是否成功 | Whether creation succeeded
        """
        if not result_file:
            result_file = self.config.output_file

        if not result_file:
            self.logger.error("未指定结果文件 | Result file not specified")
            return False

        self.logger.info("=" * 60)
        self.logger.info("开始创建可视化图表 | Starting visualization creation")
        self.logger.info("=" * 60)

        try:
            # 检查matplotlib是否可用 | Check if matplotlib is available
            import matplotlib
            matplotlib.use('Agg')  # 使用非交互式后端 | Use non-interactive backend
        except ImportError:
            self.logger.warning("matplotlib未安装，跳过可视化 | matplotlib not installed, skipping visualization")
            return False

        # 加载数据 | Load data
        if not self.analyzer:
            analyzer = SNPIndexAnalyzer(result_file, self.config)
            if not analyzer.load_data():
                return False
            data = analyzer.data
            sample_names = analyzer.sample_names
        else:
            data = self.analyzer.data
            sample_names = self.analyzer.sample_names

        if not data:
            self.logger.error("没有可用于可视化的数据 | No data available for visualization")
            return False

        # 创建可视化器 | Create visualizer
        self.visualizer = SNPIndexVisualizer(data, sample_names, self.config)

        output_prefix = self.config.get_output_path(self.config.prefix)
        success_count = 0
        total_plots = 4

        # 创建综合分析图 | Create comprehensive analysis plot
        if self.visualizer.create_comprehensive_plot(f"{output_prefix}_comprehensive.png"):
            success_count += 1

        # 创建曼哈顿图 | Create Manhattan plot
        if self.visualizer.create_manhattan_plot(f"{output_prefix}_manhattan.png"):
            success_count += 1

        # 创建分布图 | Create distribution plots
        if self.visualizer.create_distribution_plots(output_prefix):
            success_count += 1

        # 创建相关性图 | Create correlation plot
        if self.visualizer.create_correlation_plot(f"{output_prefix}_correlation.png"):
            success_count += 1

        # 创建滑动窗口折线图 | Create sliding window line plot
        if self.config.enable_sliding_window_plot:
            total_plots += 1
            if self.visualizer.create_sliding_window_plot(f"{output_prefix}_sliding_window.png"):
                success_count += 1

        # 创建多染色体分离图 | Create multi-chromosome separated plot
        if self.config.create_multi_chrom_plot:
            total_plots += 1
            if self.visualizer.create_multi_chromosome_sliding_plot(f"{output_prefix}_multi_chrom_sliding.png"):
                success_count += 1

        self.logger.info(f"可视化创建完成 | Visualization creation completed: {success_count}/{total_plots} 个图表成功")
        return success_count > 0

    def run_full_pipeline(self) -> bool:
        """
        运行完整流程 | Run full pipeline

        Returns:
            bool: 流程是否成功 | Whether pipeline succeeded
        """
        start_time = time.time()

        self.logger.info("=" * 60)
        self.logger.info("SNP Index分析流程开始 | SNP Index Analysis Pipeline Started")
        self.logger.info("=" * 60)
        self.logger.info(f"输入VCF文件 | Input VCF file: {self.config.input_vcf}")
        self.logger.info(f"输出目录 | Output directory: {self.config.output_dir}")
        self.logger.info(f"输出前缀 | Output prefix: {self.config.prefix}")

        try:
            # 步骤1: 计算SNP index | Step 1: Calculate SNP index
            if self.config.input_vcf:
                if not self.calculate_snp_index():
                    self.logger.error("SNP index计算失败 | SNP index calculation failed")
                    return False

            # 步骤2: 分析结果 | Step 2: Analyze results
            if not self.analyze_results():
                self.logger.error("结果分析失败 | Result analysis failed")
                return False

            # 步骤3: 计算和保存滑动窗口结果 | Step 3: Calculate and save sliding window results
            if not self._save_sliding_window_results():
                self.logger.warning("滑动窗口结果保存失败，但不影响主要结果 | Sliding window results save failed, but doesn't affect main results")

            # 步骤4: 创建可视化 | Step 4: Create visualizations
            if not self.create_visualizations():
                self.logger.warning("可视化创建失败，但不影响主要结果 | Visualization creation failed, but doesn't affect main results")

            # 输出总结信息 | Output summary
            elapsed_time = time.time() - start_time
            self.logger.info("=" * 60)
            self.logger.info("流程总结 | Pipeline Summary")
            self.logger.info("=" * 60)
            self.logger.info(f"总运行时间 | Total runtime: {elapsed_time:.2f} 秒 | seconds")
            self.logger.info(f"输出目录 | Output directory: {self.config.output_dir}")

            if self.calculator:
                stats = self.calculator.get_statistics()
                self.logger.info(f"处理的变异位点 | Processed variants: {stats.get('processed_variants', 0):,}")

            self.logger.info("SNP Index分析流程成功完成 | SNP Index Analysis Pipeline Completed Successfully")
            return True

        except KeyboardInterrupt:
            self.logger.warning("流程被用户中断 | Pipeline interrupted by user")
            return False
        except Exception as e:
            self.logger.error(f"流程执行出错 | Pipeline execution error: {str(e)}", exc_info=True)
            return False


def parse_arguments():
    """
    解析命令行参数 | Parse command line arguments
    """
    parser = argparse.ArgumentParser(
        description='SNP Index计算和分析工具 | SNP Index Calculation and Analysis Tool',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
使用示例 | Usage Examples:

  # 基本计算 | Basic calculation
  %(prog)s -i input.vcf.gz -o results/

  # 完整分析流程 | Complete analysis pipeline
  %(prog)s -i input.vcf.gz -o results/ -p sample1 -v

  # 自定义参数 | Custom parameters
  %(prog)s -i input.vcf.gz -o results/ --min-depth 20 --extreme-threshold 0.9

  # 只分析已有结果（包含滑动窗口）| Analyze existing results (including sliding window)
  %(prog)s --analyze-only -r results.tsv -o output/

  # 禁用滑动窗口图 | Disable sliding window plot
  %(prog)s -i input.vcf.gz -o results/ --disable-sliding-window-plot
        '''
    )

    # 必需参数 | Required arguments
    required = parser.add_argument_group('required arguments')
    required.add_argument('-i', '--input', metavar='VCF_FILE',
                         help='输入VCF文件路径 | Input VCF file path (required for calculation)')

    # 通用参数 | General arguments
    general = parser.add_argument_group('general arguments')
    general.add_argument('-o', '--output', metavar='DIR', default='./snp_index_output',
                        help='输出目录路径 | Output directory path (default: ./snp_index_output)')
    general.add_argument('-p', '--prefix', default='snp_index',
                        help='输出文件前缀 | Output file prefix (default: snp_index)')
    general.add_argument('-r', '--result-file', metavar='TSV_FILE',
                        help='已有结果文件路径 | Existing result file path (for analysis only)')

    # 过滤参数 | Filtering parameters
    filtering = parser.add_argument_group('filtering parameters')
    filtering.add_argument('--min-depth', type=int, default=10,
                          help='最小测序深度 | Minimum sequencing depth (default: 10)')
    filtering.add_argument('--min-quality', type=int, default=20,
                          help='最小质量值 | Minimum quality value (default: 20)')
    filtering.add_argument('--min-mapping-quality', type=int, default=20,
                          help='最小mapping质量 | Minimum mapping quality (default: 20)')
    filtering.add_argument('--sample-names', nargs=2, metavar=('SAMPLE1', 'SAMPLE2'),
                          help='指定要分析的样本名称 | Specify sample names to analyze')

    # 分析参数 | Analysis parameters
    analysis = parser.add_argument_group('analysis parameters')
    analysis.add_argument('--extreme-threshold', type=float, default=0.8,
                         help='极端ΔSNP index阈值 | Extreme ΔSNP index threshold (default: 0.8)')
    analysis.add_argument('--region-threshold', type=float, default=0.5,
                         help='区域检测阈值 | Region detection threshold (default: 0.5)')
    analysis.add_argument('--min-region-snps', type=int, default=5,
                         help='区域最少SNP数量 | Minimum SNPs for region (default: 5)')
    analysis.add_argument('--max-region-gap', type=int, default=10000,
                         help='区域最大gap(bp) | Maximum gap in region (default: 10000)')

    # 滑动窗口参数 | Sliding window parameters
    sliding_window = parser.add_argument_group('sliding window parameters')
    sliding_window.add_argument('--window-size', type=int, default=1000000,
                              help='滑动窗口大小(bp) | Sliding window size in bp (default: 1000000)')
    sliding_window.add_argument('--step-size', type=int, default=100000,
                              help='滑动步长(bp) | Sliding step size in bp (default: 100000)')
    sliding_window.add_argument('--min-window-snps', type=int, default=5,
                              help='窗口最少SNP数 | Minimum SNPs per window (default: 5)')
    sliding_window.add_argument('--confidence-level', type=float, default=0.95,
                              help='置信水平 | Confidence level (default: 0.95)')

    # 可视化参数 | Visualization parameters
    visualization = parser.add_argument_group('visualization parameters')
    visualization.add_argument('--disable-sliding-window-plot', action='store_true',
                              help='禁用滑动窗口折线图 | Disable sliding window line plot')
    visualization.add_argument('--create-multi-chrom-plot', action='store_true',
                              help='创建多染色体分离图 | Create multi-chromosome separated plot')

    # 模式选择 | Mode selection
    mode = parser.add_argument_group('mode selection')
    mode.add_argument('--calculate-only', action='store_true',
                     help='只计算SNP index，不分析 | Calculate SNP index only, no analysis')
    mode.add_argument('--analyze-only', action='store_true',
                     help='只分析已有结果，不计算 | Analyze existing results only, no calculation')
    mode.add_argument('--skip-visualization', action='store_true',
                     help='跳过可视化 | Skip visualization')

    # 日志参数 | Logging parameters
    log_group = parser.add_argument_group('logging options')
    log_group.add_argument('-v', '--verbose', action='count', default=0,
                         help='详细输出模式 | Verbose mode (-v: INFO, -vv: DEBUG)')
    log_group.add_argument('--quiet', action='store_true',
                         help='静默模式 | Quiet mode (ERROR only)')
    log_group.add_argument('--log-file', metavar='FILE',
                         help='日志文件路径 | Log file path')

    # 执行控制 | Execution control
    exec_group = parser.add_argument_group('execution options')
    exec_group.add_argument('-f', '--force', action='store_true',
                           help='强制覆盖已存在文件 | Force overwrite existing files')
    exec_group.add_argument('-t', '--threads', type=int, default=1,
                           help='线程数 | Number of threads (default: 1)')

    # 版本信息 | Version information
    parser.add_argument('-V', '--version', action='version',
                       version='%(prog)s 1.0.0')

    return parser.parse_args()


def main():
    """主函数 | Main function"""
    args = parse_arguments()

    # 检查参数一致性 | Check parameter consistency
    if args.calculate_only and args.analyze_only:
        print("错误: 不能同时指定 --calculate-only 和 --analyze-only", file=sys.stderr)
        sys.exit(1)

    if not args.input and not args.result_file:
        print("错误: 必须指定输入VCF文件 (-i) 或结果文件 (-r)", file=sys.stderr)
        sys.exit(1)

    if args.analyze_only and not args.result_file:
        print("错误: --analyze-only 模式需要指定结果文件 (-r)", file=sys.stderr)
        sys.exit(1)

    # 创建配置 | Create configuration
    config = SNPIndexConfig(
        input_vcf=args.input,
        output_dir=args.output,
        prefix=args.prefix,
        min_depth=args.min_depth,
        min_quality=args.min_quality,
        min_mapping_quality=args.min_mapping_quality,
        extreme_threshold=args.extreme_threshold,
        region_threshold=args.region_threshold,
        min_region_snps=args.min_region_snps,
        max_region_gap=args.max_region_gap,
        sample_names=args.sample_names,
        threads=args.threads,
        force=args.force,
        quiet=args.quiet,
        verbose=args.verbose,
        log_file=args.log_file,
        # 滑动窗口参数 | Sliding window parameters
        window_size=getattr(args, 'window_size', 1000000),
        step_size=getattr(args, 'step_size', 100000),
        min_window_snps=getattr(args, 'min_window_snps', 5),
        confidence_level=getattr(args, 'confidence_level', 0.95),
        enable_sliding_window_plot=not getattr(args, 'disable_sliding_window_plot', False),
        create_multi_chrom_plot=getattr(args, 'create_multi_chrom_plot', False)
    )

    # 创建处理器 | Create processor
    processor = SNPIndexProcessor(config)

    # 运行相应的流程 | Run appropriate pipeline
    success = False
    if args.analyze_only:
        success = processor.analyze_results(args.result_file)
        if not args.skip_visualization:
            processor.create_visualizations(args.result_file)
    elif args.calculate_only:
        success = processor.calculate_snp_index()
    else:
        success = processor.run_full_pipeline()

    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()