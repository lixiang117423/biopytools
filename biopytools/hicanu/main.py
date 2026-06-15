"""
HiCanu组装主程序模块|HiCanu Assembly Main Module
"""

import argparse
import sys
import time
from .config import HiCanuConfig
from .utils import CanuLogger, CommandRunner
from .calculator import HiCanuCalculator


class HiCanuPipeline:
    """HiCanu组装流程主类|HiCanu Assembly Pipeline Main Class"""

    def __init__(self, **kwargs):
        """
        初始化HiCanu组装流程|Initialize HiCanu assembly pipeline

        Args:
            **kwargs: 配置参数|Configuration parameters
        """
        # 初始化配置|Initialize configuration
        self.config = HiCanuConfig(**kwargs)
        self.config.validate()

        # 初始化日志|Initialize logging
        self.logger_manager = CanuLogger(self.config.log_dir)
        self.logger = self.logger_manager.get_logger()

        # 初始化命令执行器|Initialize command runner
        self.cmd_runner = CommandRunner(self.logger, self.config.raw_dir)

        # 初始化计算器|Initialize calculator
        self.calculator = HiCanuCalculator(self.config, self.logger, self.cmd_runner)

    def run(self):
        """运行组装流程|Run assembly pipeline"""
        start_time = time.time()

        try:
            # 输出流程信息|Output pipeline information
            self.logger.info("=" * 60)
            self.logger.info("HiCanu基因组组装流程|HiCanu Genome Assembly Pipeline")
            self.logger.info("=" * 60)
            self.logger.info(f"输入reads|Input reads: {self.config.reads_file}")
            self.logger.info(f"基因组大小|Genome size: {self.config.genome_size}")
            self.logger.info(f"工作目录|Work directory: {self.config.work_dir}")
            self.logger.info(f"输出前缀|Output prefix: {self.config.prefix}")
            self.logger.info("=" * 60)

            # 运行组装|Run assembly
            success = self.calculator.run_assembly()

            if not success:
                self.logger.error("组装流程失败|Assembly pipeline failed")
                sys.exit(1)

            # 输出总结信息|Output summary information
            elapsed_time = time.time() - start_time
            self.logger.info("=" * 60)
            self.logger.info("流程总结|Pipeline Summary")
            self.logger.info("=" * 60)
            self.logger.info(f"总运行时间|Total runtime: {elapsed_time:.2f} seconds ({elapsed_time/60:.2f} minutes)")
            self.logger.info(f"工作目录|Work directory: {self.config.work_dir}")

            output_files = self.config.get_output_files()
            self.logger.info("主要输出文件|Main output files:")
            self.logger.info(f"  Contigs: {output_files['contigs']}")
            self.logger.info(f"  Unitigs: {output_files['unitigs']}")
            self.logger.info(f"  报告|Report: {output_files['report']}")
            self.logger.info("=" * 60)
            self.logger.info("流程成功完成|Pipeline completed successfully")
            self.logger.info("=" * 60)

        except KeyboardInterrupt:
            self.logger.warning("流程被用户中断|Pipeline interrupted by user")
            sys.exit(130)
        except Exception as e:
            self.logger.error(f"流程执行出错|Pipeline execution error: {str(e)}", exc_info=True)
            sys.exit(1)


def main():
    """主函数|Main function"""
    parser = argparse.ArgumentParser(
        description='HiCanu基因组组装工具|HiCanu Genome Assembly Tool',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例|Examples:
  %(prog)s -i reads.fastq -g 120m -p sample1 -o output
        '''
    )

    # 必需参数|Required parameters
    required = parser.add_argument_group('必需参数|Required parameters')
    required.add_argument('-i', '--reads', required=True,
                         dest='reads_file',
                         help='输入reads文件路径(FASTA/FASTQ格式)|Input reads file path (FASTA/FASTQ format)')
    required.add_argument('-g', '--genome-size', required=True,
                         dest='genome_size',
                         help='基因组大小(如120m, 1g)|Genome size (e.g., 120m, 1g)')
    required.add_argument('-p', '--prefix', required=True,
                         help='输出文件前缀|Output file prefix')

    # 可选参数|Optional parameters
    optional = parser.add_argument_group('可选参数|Optional parameters')
    optional.add_argument('-o', '--output-dir',
                         default='./hicanu_output',
                         dest='base_dir',
                         help='输出目录路径|Output directory path')
    optional.add_argument('--canu-path',
                         default='~/miniforge3/envs/canu_v.2.3/bin/canu',
                         help='Canu可执行文件路径|Path to Canu executable')
    optional.add_argument('--min-read-length',
                         type=int, default=1000,
                         help='最小reads长度|Minimum read length')
    optional.add_argument('--min-overlap-length',
                         type=int, default=500,
                         help='最小重叠长度|Minimum overlap length')
    optional.add_argument('--corrected-error-rate',
                         type=float,
                         help='纠错后错误率|Corrected error rate')
    optional.add_argument('--raw-error-rate',
                         type=float,
                         help='原始错误率|Raw error rate')
    optional.add_argument('--max-input-coverage',
                         type=int,
                         help='最大输入覆盖度|Maximum input coverage')
    optional.add_argument('--stage',
                         choices=['haplotype', 'correct', 'trim', 'assemble', 'trim-assemble'],
                         default='assemble',
                         help='组装阶段|Assembly stage')

    # 计算资源|Computing resources
    resources = parser.add_argument_group('计算资源|Computing resources')
    resources.add_argument('-t', '--threads',
                          type=int, default=12,
                          help='线程数|Number of threads')
    resources.add_argument('-m', '--memory',
                          default='80G',
                          help='内存限制|Memory limit')
    resources.add_argument('--use-grid',
                          action='store_true',
                          help='使用网格调度|Use grid engine')
    resources.add_argument('--grid-options',
                          help='网格调度选项|Grid engine options')

    # 执行控制|Execution control
    execution = parser.add_argument_group('执行控制|Execution control')
    execution.add_argument('--dry-run',
                          action='store_true',
                          help='模拟运行(不实际执行)|Dry run (do not execute)')
    execution.add_argument('--keep-intermediate',
                          action='store_true',
                          help='保留中间文件|Keep intermediate files')
    execution.add_argument('--no-resume',
                          action='store_true',
                          help='禁用断点续传（强制重新运行所有步骤）|Disable resume mode (force rerun all steps)')
    execution.add_argument('--resume',
                          action='store_true',
                          help='启用断点续传（默认已启用）|Enable resume mode (enabled by default)')

    # 日志参数|Logging parameters
    log_group = parser.add_argument_group('日志选项|Logging options')
    log_group.add_argument('-v', '--verbose',
                          action='count', default=0,
                          help='详细输出模式(-v: INFO, -vv: DEBUG)|Verbose mode (-v: INFO, -vv: DEBUG)')
    log_group.add_argument('--quiet',
                          action='store_true',
                          help='静默模式(只输出ERROR)|Quiet mode (ERROR only)')

    # 版本信息|Version information
    parser.add_argument('-V', '--version',
                       action='version',
                       version='%(prog)s 1.0.0')

    args = parser.parse_args()

    # 确定resume参数值|Determine resume parameter value
    # 默认启用断点续传，使用--no-resume可禁用|Resume enabled by default, use --no-resume to disable
    resume_value = not args.no_resume

    # 创建并运行流程|Create and run pipeline
    try:
        pipeline = HiCanuPipeline(
            reads_file=args.reads_file,
            genome_size=args.genome_size,
            prefix=args.prefix,
            canu_path=args.canu_path,
            base_dir=args.base_dir,
            min_read_length=args.min_read_length,
            min_overlap_length=args.min_overlap_length,
            corrected_error_rate=args.corrected_error_rate,
            raw_error_rate=args.raw_error_rate,
            max_input_coverage=args.max_input_coverage,
            threads=args.threads,
            memory=args.memory,
            use_grid=args.use_grid,
            grid_options=args.grid_options,
            stage=args.stage,
            dry_run=args.dry_run,
            keep_intermediate=args.keep_intermediate,
            resume=resume_value
        )

        pipeline.run()

    except ValueError as e:
        print(f"参数错误|Parameter error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"程序错误|Program error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
