"""
RNA-Bloom转录组从头组装主程序模块|RNA-Bloom De Novo Transcriptome Assembly Main Module
"""

import argparse
import sys
from .config import RNABloomConfig
from .utils import RNABloomLogger, CommandRunner
from .assembler import TranscriptomeAssembler


class RNABloomAssembler:
    """RNA-Bloom转录组从头组装主类|Main RNA-Bloom Assembly Class"""

    def __init__(self, **kwargs):
        # 初始化配置|Initialize configuration
        self.config = RNABloomConfig(**kwargs)
        self.config.validate()

        # 初始化日志|Initialize logging
        self.logger_manager = RNABloomLogger(self.config.output_path)
        self.logger = self.logger_manager.get_logger()

        # 初始化命令执行器|Initialize command runner
        self.cmd_runner = CommandRunner(self.logger, self.config.output_path)

        # 初始化组装器|Initialize assembler
        self.assembler = TranscriptomeAssembler(
            self.config,
            self.logger,
            self.cmd_runner
        )

    def run_assembly(self):
        """运行完整的转录组组装流程|Run complete transcriptome assembly pipeline"""
        try:
            self.logger.info("RNA-Bloom转录组从头组装流程开始|RNA-Bloom de novo transcriptome assembly started")
            self.logger.info(f"组装模式|Assembly mode: {self.config.get_assembly_mode()}")
            self.logger.info(f"输入文件|Input files: {', '.join(self.config.get_input_files())}")
            self.logger.info(f"输出目录|Output directory: {self.config.output_dir}")
            self.logger.info(f"线程数|Threads: {self.config.threads}")

            # 运行组装|Run assembly
            success = self.assembler.run_assembly()

            if success:
                # 验证输出|Validate output
                if self.assembler.validate_output():
                    self.logger.info("=" * 60)
                    self.logger.info("组装成功完成！|Assembly completed successfully!")
                    self.logger.info("=" * 60)
                    self.logger.info(f"输出文件位于|Output files in: {self.config.output_dir}")

                    # 列出输出文件|List output files
                    output_files = self.assembler.get_output_files()
                    for file_type, file_path in output_files.items():
                        self.logger.info(f"  - {file_path}")

                    self.logger.info("=" * 60)
                    return 0
                else:
                    self.logger.error("输出验证失败|Output validation failed")
                    return 1
            else:
                self.logger.error("组装失败|Assembly failed")
                return 1

        except KeyboardInterrupt:
            self.logger.warning("操作被用户中断|Operation interrupted by user")
            return 130
        except Exception as e:
            self.logger.error(f"组装流程在执行过程中意外终止|Assembly pipeline terminated unexpectedly: {e}", exc_info=True)
            return 1


def main():
    """主函数|Main function"""
    parser = argparse.ArgumentParser(
        description="RNA-Bloom转录组从头组装：使用RNA-Bloom进行无参考转录组组装|RNA-Bloom de novo transcriptome assembly: Assemble transcripts without reference using RNA-Bloom",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例|Examples:
  # 短reads双端组装|Short-read paired-end assembly
  %(prog)s --left reads_1.fq --right reads_2.fq -o ./results

  # 短reads单端组装|Short-read single-end assembly
  %(prog)s --sef reads.fq -o ./results

  # 长reads组装|Long-read assembly
  %(prog)s --long long_reads.fq -o ./results

  # 长reads + 短reads混合组装|Long-read + short-read assembly
  %(prog)s --long long_reads.fq --sef short_reads.fq -o ./results

  # 单细胞混合组装|Single-cell pooled assembly
  %(prog)s --cell-list cells.txt -o ./results
        '''
    )

    # 必需参数|Required parameters
    required = parser.add_argument_group('必需参数（至少指定一种）|Required parameters (specify at least one)')

    required.add_argument("--left", "-1",
                         help="左端reads文件|Left reads file (paired-end)")
    required.add_argument("--right", "-2",
                         help="右端reads文件|Right reads file (paired-end)")
    required.add_argument("--sef",
                         help="单端正向reads文件|Single-end forward reads file")
    required.add_argument("--ser",
                         help="单端反向reads文件|Single-end reverse reads file")
    required.add_argument("--long",
                         help="长reads文件|Long reads file (ONT/PacBio)")
    required.add_argument("--cell-list",
                         help="单细胞列表文件|Single-cell list file (pooled assembly)")

    # 必需输出参数|Required output parameter
    required_output = parser.add_argument_group('输出参数|Output parameters')
    required_output.add_argument("-o", "--output-dir", required=True,
                                help="输出目录|Output directory")

    # 可选参数|Optional parameters
    optional = parser.add_argument_group('可选参数|Optional parameters')

    optional.add_argument("-t", "--threads", type=int, default=12,
                         help="线程数|Number of threads (default: 12)")

    optional.add_argument("--rnabloom-path",
                         default="rnabloom",
                         help="RNA-Bloom工具路径|RNA-Bloom tool path (default: rnabloom)")

    # Bloom filter配置|Bloom filter configuration
    bloom_group = parser.add_argument_group('Bloom filter配置|Bloom filter configuration')

    bloom_group.add_argument("--mem", type=float,
                            help="Bloom filter总大小(GB)|Total Bloom filter size in GB")

    bloom_group.add_argument("--fpr", type=float,
                            help="假阳性率|False positive rate (0-1)")

    bloom_group.add_argument("--nk", "--num-kmers", type=int, dest="num_kmers",
                            help="唯一k-mer数量|Number of unique kmers")

    # 数据类型配置|Data type configuration
    data_type = parser.add_argument_group('数据类型配置|Data type configuration')

    data_type.add_argument("--stranded", action="store_true",
                          help="链特异性数据|Strand-specific data")

    data_type.add_argument("--revcomp-left", action="store_true",
                          help="反向互补左端reads|Reverse-complement left reads")

    data_type.add_argument("--revcomp-right", action="store_true",
                          help="反向互补右端reads|Reverse-complement right reads")

    data_type.add_argument("--pacbio", action="store_true",
                          help="PacBio数据|PacBio data (default: ONT)")

    # 参考引导|Reference-guided
    reference = parser.add_argument_group('参考引导组装|Reference-guided assembly')

    reference.add_argument("--ref", "--reference",
                          help="参考转录本文件|Reference transcript file")

    # 输出选项|Output options
    output = parser.add_argument_group('输出选项|Output options')

    output.add_argument("--min-length", type=int, default=200,
                       help="最小转录本长度|Minimum transcript length (default: 200)")

    output.add_argument("--uracil", action="store_true",
                       help="输出尿嘧啶(U)而非胸腺嘧啶(T)|Write uracil (U) instead of thymine (T)")

    output.add_argument("--no-nr", action="store_false", dest="export_non_redundant",
                       help="不导出去冗余转录本|Do not export non-redundant transcripts")

    # 步骤控制|Step control
    control = parser.add_argument_group('处理控制|Processing control')

    control.add_argument("--stage", type=int, choices=[1, 2, 3],
                        help="停止阶段|Stop at stage (1-3)")

    args = parser.parse_args()

    # 创建组装器并运行|Create assembler and run
    try:
        assembler = RNABloomAssembler(
            # 输入文件|Input files
            left_reads=args.left,
            right_reads=args.right,
            single_end_forward=args.sef,
            single_end_reverse=args.ser,
            long_reads=args.long,
            cell_list=args.cell_list,

            # 输出|Output
            output_dir=args.output_dir,

            # 处理参数|Processing parameters
            threads=args.threads,
            rnabloom_path=args.rnabloom_path,

            # Bloom filter|Bloom filter
            memory_gb=args.mem,
            false_positive_rate=args.fpr,
            num_kmers=args.num_kmers,

            # 数据类型|Data type
            stranded=args.stranded,
            revcomp_left=args.revcomp_left,
            revcomp_right=args.revcomp_right,
            is_pacbio=args.pacbio,

            # 参考|Reference
            reference_transcripts=args.ref,

            # 输出选项|Output options
            min_length=args.min_length,
            write_uracil=args.uracil,
            export_non_redundant=args.export_non_redundant,

            # 步骤控制|Step control
            stage=args.stage
        )

        return assembler.run_assembly()

    except ValueError as e:
        print(f"参数错误|Parameter error: {e}", file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        print("操作被用户中断|Operation interrupted by user", file=sys.stderr)
        return 130
    except Exception as e:
        print(f"错误|Error: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
