"""
BAM覆盖度统计主程序模块|BAM Coverage Statistics Main Module
"""

import sys
import os
from .config import BAMCoverageConfig
from .utils import BAMCoverageLogger, SAMToolsHelper, CoverageDataProcessor


class BAMCoverageAnalyzer:
    """BAM覆盖度统计主类|Main BAM Coverage Statistics Class"""

    def __init__(self, **kwargs):
        # 初始化配置|Initialize configuration
        self.config = BAMCoverageConfig(**kwargs)
        self.config.validate()

        # 初始化日志|Initialize logging
        self.logger_manager = BAMCoverageLogger(
            self.config.output_path,
            verbose=self.config.verbose,
            quiet=self.config.quiet
        )
        self.logger = self.logger_manager.get_logger()

        # 初始化辅助工具|Initialize helper tools
        self.samtools = SAMToolsHelper(self.logger, threads=self.config.threads)
        self.processor = CoverageDataProcessor()

    def run_analysis(self):
        """运行完整的覆盖度分析流程|Run complete coverage analysis pipeline"""
        try:
            self.logger.info("=" * 80)
            self.logger.info("开始BAM覆盖度统计分析|Starting BAM Coverage Statistics Analysis")
            self.logger.info("=" * 80)

            # 检查samtools|Check samtools
            if not self.samtools.check_samtools():
                self.logger.error("samtools不可用，退出分析|samtools not available, exiting")
                sys.exit(1)

            # 获取BAM文件列表|Get BAM file list
            bam_files = self.config.get_bam_files()
            self.logger.info(f"找到|Found {len(bam_files)} 个BAM文件|BAM files")

            # 获取染色体长度（如果需要）|Get chromosome length (if needed)
            if self.config.end is None:
                self.logger.info("未指定终止位置，将从BAM文件获取染色体长度|End position not specified, will get chromosome length from BAM file")
                chrom_length = self.samtools.get_chromosome_length(bam_files[0], self.config.chromosome)
                if chrom_length is None:
                    self.logger.error(f"无法获取染色体长度|Cannot get chromosome length: {self.config.chromosome}")
                    sys.exit(1)
                self.config.end = chrom_length
                self.logger.info(f"染色体长度|Chromosome length: {chrom_length}")

            # 显示分析参数|Display analysis parameters
            self.logger.info(f"目标区域|Target region: {self.config.chromosome}:{self.config.start}-{self.config.end}")
            self.logger.info(f"区域大小|Region size: {self.config.end - self.config.start + 1} bp")
            self.logger.info(f"最小mapping质量|Min MAPQ: {self.config.min_mapq}")
            self.logger.info(f"最小碱基质量|Min base quality: {self.config.min_baseq}")

            # 提取覆盖度|Extract coverage
            self.logger.info("")
            self.logger_manager.step("步骤1/3: 提取覆盖度数据|Step 1/3: Extracting coverage data")

            temp_dir = self.config.output_path / 'temp'
            temp_dir.mkdir(exist_ok=True)

            coverage_files = []
            sample_names = []

            for idx, bam_file in enumerate(bam_files, 1):
                sample_name = os.path.basename(bam_file).replace('.bam', '')
                sample_names.append(sample_name)

                self.logger.info(f"[{idx}/{len(bam_files)}] 处理样本|Processing sample: {sample_name}")

                # 检查索引|Check index
                if not self.samtools.check_bam_index(bam_file):
                    self.logger.info(f"  创建索引|Creating index for {sample_name}")
                    if not self.samtools.create_bam_index(bam_file):
                        self.logger.warning(f"  索引创建失败，跳过该样本|Index creation failed, skipping {sample_name}")
                        continue

                # 提取覆盖度|Extract coverage
                output_file = temp_dir / f"{sample_name}.depth"

                if self.samtools.extract_coverage(
                    bam_file,
                    self.config.chromosome,
                    self.config.start,
                    self.config.end,
                    self.config.min_mapq,
                    self.config.min_baseq,
                    str(output_file)
                ):
                    coverage_files.append(str(output_file))

                    # 显示统计信息|Display statistics
                    line_count = sum(1 for _ in open(output_file))
                    self.logger.info(f"  成功：{line_count} 个位置|Success: {line_count} positions")
                else:
                    self.logger.warning(f"  提取失败，跳过该样本|Extraction failed, skipping {sample_name}")

            if not coverage_files:
                self.logger.error("没有成功提取任何样本的覆盖度|No samples successfully extracted")
                sys.exit(1)

            self.logger.info(f"成功提取|Successfully extracted: {len(coverage_files)}/{len(bam_files)} 个样本|samples")

            # 合并覆盖度数据|Merge coverage data
            if self.config.merge_output and len(coverage_files) > 1:
                self.logger.info("")
                self.logger_manager.step("步骤2/3: 合并覆盖度数据|Step 2/3: Merging coverage data")

                merged_file = self.config.output_path / f"{self.config.output_prefix}_merged.txt"

                if self.processor.merge_coverage_data(coverage_files, sample_names, str(merged_file), self.logger):
                    self.logger.info(f"合并文件已保存|Merged file saved: {merged_file}")
                else:
                    self.logger.error("合并失败|Merge failed")
                    sys.exit(1)
            elif len(coverage_files) == 1:
                # 单个样本，直接重命名|Single sample, just rename
                import shutil
                merged_file = self.config.output_path / f"{self.config.output_prefix}_merged.txt"
                shutil.copy(coverage_files[0], merged_file)
                self.logger.info(f"单个样本，输出文件|Single sample, output file: {merged_file}")

            # 生成统计摘要|Generate summary statistics
            if self.config.generate_summary:
                self.logger.info("")
                self.logger_manager.step("步骤3/3: 生成统计摘要|Step 3/3: Generating summary statistics")

                summary_file = self.config.output_path / f"{self.config.output_prefix}_summary.txt"

                if self.processor.calculate_summary_statistics(str(merged_file), str(summary_file), self.logger):
                    self.logger.info(f"统计摘要已保存|Summary statistics saved: {summary_file}")

                    # 显示摘要|Display summary
                    with open(summary_file, 'r') as f:
                        self.logger.info("")
                        self.logger.info("统计摘要|Summary Statistics:")
                        for line in f:
                            self.logger.info(line.strip())
                else:
                    self.logger.warning("统计摘要生成失败|Failed to generate summary statistics")

            # 完成|Complete
            self.logger.info("")
            self.logger.info("=" * 80)
            self.logger.info("BAM覆盖度统计分析完成|BAM Coverage Statistics Analysis Completed")
            self.logger.info("=" * 80)
            self.logger.info(f"输出目录|Output directory: {self.config.output_dir}")
            self.logger.info(f"  合并覆盖度文件|Merged coverage file: {merged_file.name}")
            if self.config.generate_summary:
                self.logger.info(f"  统计摘要文件|Summary statistics file: {summary_file.name}")
            self.logger.info(f"  临时文件目录|Temp files directory: {temp_dir}")

        except Exception as e:
            self.logger.error(f"分析过程中发生错误|Error during analysis: {e}")
            import traceback
            self.logger.error(traceback.format_exc())
            sys.exit(1)


def main():
    """主函数|Main function"""
    import argparse

    parser = argparse.ArgumentParser(
        description="BAM覆盖度统计工具|BAM Coverage Statistics Tool",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # 必需参数|Required arguments
    required = parser.add_argument_group('必需参数|Required arguments')
    required.add_argument('-i', '--input',
                         required=True,
                         help='输入路径（BAM文件或包含BAM的目录）|Input path (BAM file or directory containing BAM files)')
    required.add_argument('-c', '--chromosome',
                         required=True,
                         help='染色体名称|Chromosome name (e.g., chr1, Chr12)')
    required.add_argument('-s', '--start',
                         required=True,
                         type=int,
                         help='起始位置|Start position (1-based)')

    # 输出配置|Output configuration
    output = parser.add_argument_group('输出配置|Output configuration')
    output.add_argument('-o', '--output-dir',
                       default='./bam_coverage_stats_output',
                       help='输出目录|Output directory')
    output.add_argument('-p', '--output-prefix',
                       default='coverage',
                       help='输出文件前缀|Output file prefix')

    # 位置参数|Position parameters
    position = parser.add_argument_group('位置参数|Position parameters')
    position.add_argument('-e', '--end',
                         type=int,
                         help='终止位置|End position')

    # 过滤参数|Filtering parameters
    filtering = parser.add_argument_group('过滤参数|Filtering parameters')
    filtering.add_argument('--min-mapq',
                          type=int,
                          default=0,
                          help='最小mapping质量|Minimum mapping quality')
    filtering.add_argument('--min-baseq',
                          type=int,
                          default=0,
                          help='最小碱基质量|Minimum base quality')

    # 输出选项|Output options
    output_opts = parser.add_argument_group('输出选项|Output options')
    output_opts.add_argument('--no-merge',
                            action='store_true',
                            help='不合并样本输出|Do not merge sample outputs')
    output_opts.add_argument('--no-summary',
                            action='store_true',
                            help='不生成统计摘要|Do not generate summary statistics')

    # 日志选项|Logging options
    logging = parser.add_argument_group('日志选项|Logging options')
    logging.add_argument('-v', '--verbose',
                        action='count',
                        default=0,
                        help='增加输出详细程度|Increase output verbosity')
    logging.add_argument('--quiet',
                        action='store_true',
                        help='静默模式|Quiet mode')
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
                         help='试运行模式|Dry run mode')
    advanced.add_argument('-t', '--threads',
                         type=int,
                         default=64,
                         help='线程数|Number of threads')

    args = parser.parse_args()

    # 创建分析器并运行|Create analyzer and run
    try:
        analyzer = BAMCoverageAnalyzer(
            chromosome=args.chromosome,
            start=args.start,
            end=args.end,
            input=args.input,
            output_dir=args.output_dir,
            output_prefix=args.output_prefix,
            min_mapq=args.min_mapq,
            min_baseq=args.min_baseq,
            merge_output=not args.no_merge,
            generate_summary=not args.no_summary,
            verbose=(args.verbose > 0),
            quiet=args.quiet,
            log_file=args.log_file,
            log_level=args.log_level,
            dry_run=args.dry_run,
            threads=args.threads
        )

        analyzer.run_analysis()

    except ValueError as e:
        print(f"参数错误|Parameter error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"错误|Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
