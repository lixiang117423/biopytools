"""
FASTP质控主程序模块|FASTP Quality Control Main Module
"""

import argparse
import sys
import time
from .config import FastpConfig
from .utils import FastpLogger, CommandRunner
from .data_processing import SampleFinder
from .processing import FastpCore
from .results import SummaryGenerator

# 版本信息|Version information
VERSION = "1.0.0"


class FastpProcessor:
    """FASTP质控主类|Main FASTP Quality Control Class"""

    def __init__(self, **kwargs):
        # 初始化配置|Initialize configuration
        self.config = FastpConfig(**kwargs)
        self.config.validate()

        # 初始化日志|Initialize logging
        self.logger_manager = FastpLogger(
            self.config.output_path,
            log_name="fastp_processing.log",
            log_level=self.config.log_level,
            quiet=self.config.quiet
        )
        self.logger = self.logger_manager.get_logger()

        # 初始化命令执行器|Initialize command runner
        self.cmd_runner = CommandRunner(self.logger)

        # 初始化各个处理器|Initialize processors
        self.sample_finder = SampleFinder(self.config, self.logger)
        self.fastp_core = FastpCore(self.config, self.logger, self.cmd_runner)
        self.summary_generator = SummaryGenerator(self.config, self.logger)

    def run_batch_processing(self):
        """运行批处理|Run batch processing"""

        self.logger.info("=" * 60)
        self.logger.info("开始FASTQ数据质控批处理|Starting FASTQ data quality control batch processing")
        self.logger.info("=" * 60)
        self.logger.info(f"输入|Input: {'文件|File' if self.config.is_single_file else '目录|Directory'}: {self.config.input_dir}")
        self.logger.info(f"输出目录|Output directory: {self.config.output_dir}")

        # 只在目录模式下显示文件模式|Only show file pattern in directory mode
        if not self.config.is_single_file:
            self.logger.info(f"数据模式|Data mode: {'单末端|Single-end' if self.config.single_end else '双末端|Paired-end'}")
            self.logger.info(f"文件模式|File pattern: *{self.config.read1_suffix}" +
                            (f", *{self.config.read2_suffix}" if not self.config.single_end else ""))

        self.logger.info(f"线程数|Threads: {self.config.threads}")
        self.logger.info(f"质量阈值|Quality threshold: {self.config.quality_threshold}")
        self.logger.info(f"最小长度|Minimum length: {self.config.min_length}")
        self.logger.info("=" * 60)

        # 验证fastp可执行性|Validate fastp executable
        if not self.fastp_core.validate_fastp():
            sys.exit(1)

        # 创建输出目录|Create output directories
        self.fastp_core.create_output_directories()

        # 查找样本配对|Find sample pairs
        sample_pairs = self.sample_finder.find_sample_pairs()

        # 验证样本配对|Validate sample pairs
        if not self.sample_finder.validate_sample_pairs(sample_pairs):
            sys.exit(1)

        self.logger.info(f"找到 {len(sample_pairs)} 个有效样本配对|Found {len(sample_pairs)} valid sample pairs:")
        for sample_name, _, _ in sample_pairs:
            self.logger.info(f"  - {sample_name}")

        # 处理所有样本|Process all samples
        successful_count = 0
        failed_count = 0

        for sample_name, read1_file, read2_file in sample_pairs:
            if self.fastp_core.process_sample(sample_name, read1_file, read2_file):
                successful_count += 1
            else:
                failed_count += 1

        # 生成总结报告|Generate summary report
        self.summary_generator.generate_summary_report(
            successful_count, failed_count, len(sample_pairs)
        )

        # 输出最终统计|Output final statistics
        self.logger.info("=" * 60)
        self.logger.info("FASTQ质控批处理完成|FASTQ quality control batch processing completed!")
        self.logger.info(f"总样本数|Total samples: {len(sample_pairs)}")
        self.logger.info(f"成功处理|Successfully processed: {successful_count}")
        self.logger.info(f"失败样本|Failed samples: {failed_count}")
        self.logger.info(f"成功率|Success rate: {(successful_count/len(sample_pairs))*100:.1f}%")
        self.logger.info(f"质控后的清洁数据位于|Clean data location: {self.config.output_dir}")
        self.logger.info(f"质控报告位于|QC reports location: {self.config.report_path}")
        self.logger.info("=" * 60)


def main():
    """主函数|Main function"""
    start_time = time.time()

    parser = argparse.ArgumentParser(
        description="FASTQ数据质控批处理脚本|FASTQ Data Quality Control Batch Processing Script",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例|Examples:
  %(prog)s -i raw_data/ -o clean_data/
        '''
    )

    # 必需参数|Required arguments
    required = parser.add_argument_group('必需参数|Required arguments')
    required.add_argument("-i", "--input", required=True,
                         help="输入原始FASTQ数据目录|Input raw FASTQ data directory")
    required.add_argument("-o", "--output-dir", required=True,
                         help="输出清洁FASTQ数据目录|Output clean FASTQ data directory")

    # 可选参数|Optional arguments
    optional = parser.add_argument_group('可选参数|Optional arguments')
    optional.add_argument("--fastp-path", default="fastp",
                         help="fastp可执行文件路径|fastp executable path")
    optional.add_argument("-t", "--threads", type=int, default=12,
                         help="线程数|Number of threads")
    optional.add_argument("-q", "--quality-threshold", type=int, default=30,
                         help="质量阈值|Quality threshold")
    optional.add_argument("-l", "--min-length", type=int, default=50,
                         help="最小长度|Minimum length")
    optional.add_argument("-u", "--unqualified-percent", type=int, default=40,
                         help="不合格碱基百分比阈值|Unqualified base percentage threshold")
    optional.add_argument("-n", "--n-base-limit", type=int, default=10,
                         help="N碱基数量限制|N base count limit")
    optional.add_argument("--read1-suffix", default="_1.fq.gz",
                         help="Read1文件后缀（单末端模式也使用此参数）|Read1 file suffix (also used for single-end mode)")
    optional.add_argument("--read2-suffix", default="_2.fq.gz",
                         help="Read2文件后缀|Read2 file suffix")
    optional.add_argument("--single-end", action="store_true",
                         help="单末端模式|Single-end mode")

    # 日志参数|Logging parameters
    log_group = parser.add_argument_group('日志选项|Logging options')
    log_group.add_argument("-v", "--verbose", action="count", default=0,
                          help="详细输出模式(-v: INFO, -vv: DEBUG)|Verbose mode (-v: INFO, -vv: DEBUG)")
    log_group.add_argument("--quiet", action="store_true",
                          help="静默模式(只输出ERROR)|Quiet mode (ERROR only)")
    log_group.add_argument("--log-level",
                          help="日志级别(DEBUG/INFO/WARNING/ERROR/CRITICAL)|Log level")
    log_group.add_argument("--log-file",
                          help="日志文件路径|Log file path")

    # 执行控制|Execution control
    exec_group = parser.add_argument_group('执行选项|Execution options')
    exec_group.add_argument("-f", "--force", action="store_true",
                           help="强制覆盖已存在文件|Force overwrite existing files")
    exec_group.add_argument("--dry-run", action="store_true",
                           help="模拟运行(不实际执行)|Dry run without execution")

    # 版本信息|Version information
    parser.add_argument("-V", "--version", action="version",
                       version=f'%(prog)s {VERSION}')

    args = parser.parse_args()

    # 确定日志级别|Determine log level
    if args.log_level:
        log_level = args.log_level
    elif args.verbose >= 2:
        log_level = "DEBUG"
    elif args.verbose == 1:
        log_level = "INFO"
    elif args.quiet:
        log_level = "ERROR"
    else:
        log_level = "INFO"

    # 创建处理器并运行|Create processor and run
    processor = FastpProcessor(
        input_dir=args.input,
        output_dir=args.output_dir,
        fastp_path=args.fastp_path,
        threads=args.threads,
        quality_threshold=args.quality_threshold,
        min_length=args.min_length,
        unqualified_percent=args.unqualified_percent,
        n_base_limit=args.n_base_limit,
        read1_suffix=args.read1_suffix,
        read2_suffix=args.read2_suffix,
        single_end=args.single_end,
        log_level=log_level,
        quiet=args.quiet,
        verbose=args.verbose,
        force=args.force,
        dry_run=args.dry_run
    )

    try:
        # 输出程序信息|Output program information
        processor.logger.info("=" * 60)
        processor.logger.info("Program: FASTQ Quality Control Batch Processing")
        processor.logger.info(f"Version: {VERSION}")
        processor.logger.info("=" * 60)

        if args.dry_run:
            processor.logger.info("模拟运行模式-不会实际执行命令|DRY RUN mode - commands will not be executed")

        # 执行批处理|Run batch processing
        processor.run_batch_processing()

        # 输出总结信息|Output summary
        elapsed_time = time.time() - start_time
        processor.logger.info("=" * 60)
        processor.logger.info("Pipeline Summary")
        processor.logger.info("=" * 60)
        processor.logger.info(f"Total runtime: {elapsed_time:.2f} seconds")
        processor.logger.info(f"Output directory: {args.output_dir}")
        processor.logger.info("Pipeline completed successfully")

    except KeyboardInterrupt:
        processor.logger.warning("用户中断程序执行|Process interrupted by user")
        sys.exit(130)
    except Exception as e:
        processor.logger.critical(f"Pipeline failed: {str(e)}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
