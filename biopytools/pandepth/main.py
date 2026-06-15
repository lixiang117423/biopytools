"""
PanDepth覆盖度计算主程序模块|PanDepth Coverage Calculation Main Module
"""

import argparse
import sys
import os
from .config import PanDepthConfig
from .utils import PanDepthLogger, CommandRunner
from .calculator import PanDepthCalculator
from .results import PanDepthResultsMerger


class PanDepthRunner:
    """PanDepth覆盖度计算运行器|PanDepth Coverage Calculation Runner"""

    def __init__(self, **kwargs):
        # 初始化配置|Initialize configuration
        self.config = PanDepthConfig(**kwargs)
        self.config.validate()

        # 初始化日志|Initialize logging
        self.logger_manager = PanDepthLogger(
            self.config.output_path,
            verbose=self.config.verbose,
            quiet=self.config.quiet
        )
        self.logger = self.logger_manager.get_logger()

        # 初始化命令执行器|Initialize command runner
        self.cmd_runner = CommandRunner(self.logger)

        # 初始化计算器|Initialize calculator
        self.calculator = PanDepthCalculator(self.config, self.logger, self.cmd_runner)

        # 初始化结果合并器|Initialize results merger
        self.results_merger = PanDepthResultsMerger(self.logger)

    def run(self):
        """运行覆盖度计算|Run coverage calculation"""
        try:
            self.logger.info("PanDepth覆盖度计算开始|PanDepth coverage calculation started")
            self.logger.info(f"输入路径|Input path: {self.config.input_path}")
            self.logger.info(f"输出目录|Output directory: {self.config.output_dir}")
            self.logger.info(f"线程数|Threads: {self.config.threads}")

            if self.config.is_batch_mode:
                self.logger.info("批量模式|Batch mode: processing multiple BAM files")
            else:
                self.logger.info("单文件模式|Single file mode: processing one BAM file")

            if self.config.gff_file:
                self.logger.info(f"GFF/GTF文件|GFF/GTF file: {self.config.gff_file}")
                self.logger.info(f"特征类型|Feature type: {self.config.feature_type}")

            if self.config.bed_file:
                self.logger.info(f"BED文件|BED file: {self.config.bed_file}")

            if self.config.window_size:
                self.logger.info(f"窗口大小|Window size: {self.config.window_size} bp")

            # 运行计算|Run calculation
            results = self.calculator.run()

            # 合并结果|Merge results
            self.logger.info("\n" + "=" * 60)
            self.logger.info("合并统计结果|Merging statistics results")
            self.logger.info("=" * 60)

            # 根据是否有GFF文件决定合并方式|Decide merge method based on whether GFF file is used
            has_gff = self.config.gff_file is not None
            merge_success = self.results_merger.merge_results(self.config.output_dir, has_gff)

            if merge_success:
                self.logger.info("结果合并成功|Results merged successfully")
            else:
                self.logger.warning("结果合并失败或跳过|Results merging failed or skipped")

            # 输出总结信息|Output summary
            self.logger.info("\n" + "=" * 60)
            self.logger.info("覆盖度计算完成|Coverage calculation completed")
            self.logger.info("=" * 60)

            return results

        except KeyboardInterrupt:
            self.logger.warning("操作被用户中断|Operation interrupted by user")
            sys.exit(130)
        except Exception as e:
            self.logger.error(f"覆盖度计算失败|Coverage calculation failed: {e}", exc_info=True)
            sys.exit(1)


def main():
    """主函数|Main function"""
    parser = argparse.ArgumentParser(
        description="PanDepth覆盖度计算工具|PanDepth Coverage Calculation Tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例|Examples:
  %(prog)s -i sample.bam -o coverage_results
  %(prog)s -i bam_files/ -o coverage_results -t 24
  %(prog)s -i sample.bam -g genes.gff -o gene_coverage
  %(prog)s -i sample.bam -b regions.bed -o region_coverage
        '''
    )

    # 必需参数|Required parameters
    required = parser.add_argument_group('必需参数|Required parameters')
    required.add_argument("-i", "--input",
                         required=True,
                         help="输入BAM文件或BAM文件目录|Input BAM file or BAM file directory")
    required.add_argument("-o", "--output",
                         required=True,
                         help="输出目录|Output directory")

    # 目标区域选项|Target region options
    target = parser.add_argument_group('目标区域选项|Target region options')
    target.add_argument("-g", "--gff",
                        help="GFF/GTF文件用于基因覆盖度|GFF/GTF file for gene coverage")
    target.add_argument("-f", "--feature",
                        default="CDS",
                        choices=['CDS', 'exon'],
                        help="GFF/GTF特征类型|GFF/GTF feature type (default: CDS)")
    target.add_argument("-b", "--bed",
                        help="BED文件用于特定区域覆盖度|BED file for specific region coverage")
    target.add_argument("-w", "--window",
                        type=int,
                        help="滑动窗口大小(bp)|Sliding window size in bp")

    # 过滤选项|Filter options
    filter_opts = parser.add_argument_group('过滤选项|Filter options')
    filter_opts.add_argument("-q", "--min-mapq",
                            type=int,
                            default=0,
                            help="最小比对质量|Minimum mapping quality (default: 0)")
    filter_opts.add_argument("-d", "--min-depth",
                            type=int,
                            default=1,
                            help="最小深度用于统计|Minimum depth for statistics (default: 1)")
    filter_opts.add_argument("-x", "--exclude-flag",
                            type=int,
                            default=1796,
                            help="排除reads的FLAG标志|FLAG bits to exclude reads (default: 1796)")

    # 其他选项|Other options
    other = parser.add_argument_group('其他选项|Other options')
    other.add_argument("-t", "--threads",
                      type=int,
                      default=12,
                      help="线程数|Number of threads (default: 12)")
    other.add_argument("-r", "--reference",
                      help="参考基因组文件(用于CRAM解码或GC计算)|Reference genome file for CRAM decode or GC calculation")
    other.add_argument("-c", "--enable-gc",
                      action='store_true',
                      help="启用GC含量计算|Enable GC content calculation")
    other.add_argument("-a", "--all-sites",
                      action='store_true',
                      help="输出所有位点深度|Output all site depths")
    other.add_argument("--pandepth-path",
                      default='~/software/PanDepth-2.26-Linux-x86_64/pandepth',
                      help="PanDepth程序路径|PanDepth program path")

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

    args = parser.parse_args()

    # 创建运行器并运行|Create runner and run
    try:
        runner = PanDepthRunner(
            input_path=args.input,
            output_dir=args.output,
            pandepth_path=args.pandepth_path,
            gff_file=args.gff,
            bed_file=args.bed,
            window_size=args.window,
            feature_type=args.feature,
            min_mapq=args.min_mapq,
            min_depth=args.min_depth,
            exclude_flag=args.exclude_flag,
            threads=args.threads,
            reference=args.reference,
            enable_gc=args.enable_gc,
            output_all_sites=args.all_sites,
            verbose=(args.verbose > 0),
            quiet=args.quiet,
            log_file=args.log_file,
            log_level=args.log_level
        )

        runner.run()

    except ValueError as e:
        print(f"参数错误|Parameter error: {e}", file=sys.stderr)
        sys.exit(1)
    except KeyboardInterrupt:
        print("操作被用户中断|Operation interrupted by user", file=sys.stderr)
        sys.exit(130)
    except Exception as e:
        print(f"错误|Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
