"""
三代转录组比对主程序模块|Long RNA-seq Alignment Main Module
"""

import argparse
import sys
import time
from pathlib import Path

from .config import LongRNASeqConfig
from .utils import LongRNASeqLogger, check_tool
from .align import LongRNASeqAligner


def parse_arguments():
    """解析命令行参数|Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='三代转录组比对工具|Long RNA-seq Alignment Tool',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例|Example:
%(prog)s -i input.bam -r genome.fa -o output_dir
        '''
    )

    parser.add_argument('-i', '--input-file',
                        required=True,
                        help='输入文件或文件夹（BAM/FASTQ）|Input file or directory (BAM/FASTQ)')

    parser.add_argument('-r', '--ref-genome',
                        required=True,
                        help='参考基因组文件|Reference genome file')

    parser.add_argument('-o', '--output-dir',
                        required=True,
                        help='输出目录|Output directory')

    parser.add_argument('-s', '--sample-name',
                        required=False,
                        default=None,
                        help='样本名称 (可选，默认从输入文件名提取)|Sample name (optional, auto-extracted from input filename)')

    parser.add_argument('-t', '--threads',
                        type=int,
                        default=64,
                        help='线程数 (默认: 64)|Number of threads (default: 64)')

    parser.add_argument('--max-intron',
                        type=int,
                        default=100000,
                        help='最大intron长度，默认: 100000|Maximum intron length, default: 100000')

    parser.add_argument('--min-mapq',
                        type=int,
                        default=20,
                        help='最小mapping quality，默认: 20|Minimum mapping quality, default: 20')

    parser.add_argument('--no-secondary',
                        action='store_true',
                        help='不输出次优比对|Do not output secondary alignments')

    parser.add_argument('--minimap2-path',
                        default='minimap2',
                        help='minimap2可执行文件路径 (默认: minimap2)|minimap2 executable path (default: minimap2)')

    parser.add_argument('--samtools-path',
                        default='samtools',
                        help='samtools可执行文件路径 (默认: samtools)|samtools executable path (default: samtools)')

    parser.add_argument('--version',
                        action='version',
                        version='%(prog)s 1.3.0')

    return parser.parse_args()


def main():
    """主函数|Main function"""
    args = parse_arguments()

    try:
        # 创建基础配置|Create base configuration (用于检测是否为文件夹|for detecting if directory)
        base_config = LongRNASeqConfig(
            input_file=args.input_file,
            ref_genome=args.ref_genome,
            output_dir=args.output_dir,
            sample_name=args.sample_name,
            threads=args.threads,
            max_intron=args.max_intron,
            min_mapq=args.min_mapq,
            secondary=not args.no_secondary,
            minimap2_path=args.minimap2_path,
            samtools_path=args.samtools_path
        )

        base_config.validate()

        # 如果是文件夹，批量处理|If directory, batch process
        if base_config.is_directory:
            return _process_directory(base_config, args)

        # 如果是单个文件，处理单个文件|If single file, process single file
        return _process_single_file(base_config)

    except KeyboardInterrupt:
        print("用户中断|User interrupted")
        sys.exit(1)
    except ValueError as e:
        print(f"配置错误|Configuration error: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"程序执行出错|Program execution error: {str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


def _process_directory(base_config, args):
    """
    批量处理文件夹内的所有文件|Process all files in directory

    Args:
        base_config: 基础配置|Base configuration
        args: 命令行参数|Command line arguments
    """
    import datetime

    input_files = base_config.input_files
    total_files = len(input_files)

    def log_print(msg):
        """带时间戳的输出|Output with timestamp"""
        timestamp = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        print(f"{timestamp} - INFO - {msg}")

    log_print("=" * 60)
    log_print(f"文件夹模式|Directory mode")
    log_print(f"找到|Found {total_files} 个文件|files")
    log_print("=" * 60)

    success_count = 0
    failed_files = []

    for idx, input_file in enumerate(input_files, 1):
        log_print(f"处理文件 {idx}/{total_files}|Processing file {idx}/{total_files}: {Path(input_file).name}")
        log_print("-" * 60)

        try:
            # 为每个文件创建单独的配置|Create separate config for each file
            file_config = LongRNASeqConfig(
                input_file=input_file,
                ref_genome=base_config.ref_genome,
                output_dir=base_config.output_dir,
                sample_name=None,  # 自动从文件名提取|Auto-extract from filename
                threads=base_config.threads,
                max_intron=base_config.max_intron,
                min_mapq=base_config.min_mapq,
                secondary=base_config.secondary,
                minimap2_path=base_config.minimap2_path,
                samtools_path=base_config.samtools_path
            )

            # 处理单个文件|Process single file
            if _process_single_file(file_config):
                success_count += 1
            else:
                failed_files.append(input_file)

        except Exception as e:
            log_print(f"处理文件失败|Failed to process file: {e}")
            failed_files.append(input_file)

    # 打印汇总|Print summary
    log_print("=" * 60)
    log_print("批量处理完成|Batch processing completed")
    log_print("=" * 60)
    log_print(f"成功|Success: {success_count}/{total_files}")
    if failed_files:
        log_print(f"失败|Failed: {len(failed_files)}/{total_files}")
        log_print("失败文件列表|Failed files:")
        for f in failed_files:
            log_print(f"  - {f}")

    return 0 if not failed_files else 1


def _process_single_file(config):
    """
    处理单个文件|Process single file

    Args:
        config: 配置对象|Configuration object

    Returns:
        bool: 是否成功|Whether successful
    """
    try:
        # 创建日志|Create logger
        logger_manager = LongRNASeqLogger(config.log_file)
        logger = logger_manager.get_logger()

        # 记录输入类型|Log input type
        logger.info(f"检测到输入类型|Detected input type: {config.input_type}")
        if config.input_type == 'fastq_paired':
            r2_file = config._find_r2_file()
            logger.info(f"配对的R2文件|Paired R2 file: {r2_file}")

        # 检查依赖工具|Check dependency tools
        logger.info("检查依赖工具|Checking dependency tools")
        if not check_tool(config.minimap2_path, "minimap2", logger):
            return False
        if not check_tool(config.samtools_path, "samtools", logger):
            return False

        # 检查参考基因组索引|Check reference genome index
        ref_fai = Path(config.ref_genome).with_suffix('.fa.fai')
        if not ref_fai.exists():
            ref_fai = Path(config.ref_genome + '.fai')

        if not ref_fai.exists():
            logger.warning("参考基因组索引不存在，正在创建|Reference genome index not found, creating...")
            from .utils import run_command
            cmd = [config.samtools_path, "faidx", config.ref_genome]
            try:
                run_command(cmd, logger, check=True)
                logger.info("索引创建完成|Index created")
            except Exception as e:
                logger.error(f"创建索引失败|Failed to create index: {e}")
                return False

        # 记录开始时间|Record start time
        start_time = time.time()

        # 运行比对流程|Run alignment pipeline
        aligner = LongRNASeqAligner(config, logger)
        success = aligner.align()

        # 计算运行时间|Calculate runtime
        end_time = time.time()
        runtime = end_time - start_time
        runtime_min = int(runtime // 60)
        runtime_sec = int(runtime % 60)

        if success:
            logger.info("=" * 60)
            logger.info("分析完成|Analysis completed successfully")
            logger.info(f"总运行时间|Total runtime: {runtime_min}分{runtime_sec}秒|minutes{runtime_sec}seconds")
            logger.info("=" * 60)

            # 生成汇总报告|Generate summary report
            _generate_summary_report(config, runtime_min, runtime_sec)

            return True
        else:
            logger.error("分析失败|Analysis failed")
            return False

    except Exception as e:
        print(f"处理出错|Processing error: {str(e)}")
        import traceback
        traceback.print_exc()
        return False


def _generate_summary_report(config, runtime_min: int, runtime_sec: int):
    """
    生成汇总报告|Generate summary report

    Args:
        config: 配置对象|Configuration object
        runtime_min: 运行时间（分钟）|Runtime (minutes)
        runtime_sec: 运行时间（秒）|Runtime (seconds)
    """
    summary_file = config.result_dir / "alignment_summary.txt"

    try:
        with open(summary_file, 'w', encoding='utf-8') as f:
            f.write("=" * 60 + "\n")
            f.write("三代转录组比对汇总报告|Long RNA-seq Alignment Summary Report\n")
            f.write("=" * 60 + "\n\n")

            f.write(f"分析日期|Analysis date: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"运行时间|Runtime: {runtime_min}分|minutes{runtime_sec}秒|seconds\n\n")

            f.write(f"参考基因组|Reference genome: {config.ref_genome}\n")
            f.write(f"样本名称|Sample name: {config.sample_name}\n")
            f.write(f"输入文件|Input file: {config.input_file}\n")
            f.write(f"输入类型|Input type: {config.input_type}\n")
            f.write(f"线程数|Threads: {config.threads}\n")
            f.write(f"最大Intron长度|Max intron length: {config.max_intron} bp\n")
            f.write(f"最小MAPQ|Min MAPQ: {config.min_mapq}\n\n")

            f.write("-" * 60 + "\n")
            f.write("比对参数|Alignment Parameters\n")
            f.write("-" * 60 + "\n")
            params = config.get_minimap2_params()
            f.write(f"{' '.join(params)}\n\n")

            f.write("-" * 60 + "\n")
            f.write("输出文件|Output Files\n")
            f.write("-" * 60 + "\n")
            f.write(f"  BAM文件|BAM file: {config.result_dir / f'{config.sample_name}.sorted.bam'}\n")
            f.write(f"  BAM索引|BAM index: {config.result_dir / f'{config.sample_name}.sorted.bam.bai'}\n")
            f.write(f"  统计文件|Statistics: {config.stats_dir / f'{config.sample_name}_stats.txt'}\n")
            f.write(f"  详细统计|Detail stats: {config.stats_dir / f'{config.sample_name}_detail_stats.txt'}\n")
            f.write(f"  日志文件|Log file: {config.log_file}\n\n")

            f.write("-" * 60 + "\n")
            f.write("下游分析建议|Downstream Analysis Suggestions\n")
            f.write("-" * 60 + "\n")
            f.write("1. 转录本组装|Transcript assembly: stringtie/cufflinks\n")
            f.write("2. 基因表达定量|Gene expression quantification: featureCounts/HTSeq\n")
            f.write("3. 可变剪接分析|Alternative splicing: rMATS/SUPPA\n")
            f.write("4. 融合基因检测|Fusion gene detection: STAR-Fusion/FusionCatcher\n")
            f.write("5. 长reads特异分析|Long-read specific analysis: TALON/SQANTI\n")

            f.write("\n" + "=" * 60 + "\n")

        print(f"汇总报告已生成|Summary report generated: {summary_file}")
        print(f"查看汇总报告|View summary report: cat {summary_file}")

    except Exception as e:
        print(f"生成汇总报告失败|Failed to generate summary report: {e}")


if __name__ == "__main__":
    main()
