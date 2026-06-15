"""插入检测主函数|Insert detection main function"""

import argparse
import sys
from pathlib import Path

from .config import InsertDetectionConfig
from .utils import InsertDetectionLogger, check_dependencies, identify_samples
from .detector import InsertDetector
from .output import write_results, write_summary
from .pipeline_info import PipelineInfoGenerator


def parse_arguments():
    """解析命令行参数|Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="插入序列位点检测工具|Insert sequence insertion site detection tool",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    # 必需参数|Required parameters
    parser.add_argument(
        '-i', '--genome',
        required=True,
        help='参考基因组FASTA文件|Reference genome FASTA file'
    )

    parser.add_argument(
        '--insert',
        required=True,
        help='插入序列FASTA文件|Insert sequence FASTA file'
    )

    parser.add_argument(
        '--fastq-dir',
        required=True,
        help='FASTQ文件目录|FASTQ files directory'
    )

    parser.add_argument(
        '-o', '--output-dir',
        default='./insert_detection_output',
        help='输出目录|Output directory (default: ./insert_detection_output)'
    )

    # 可选参数|Optional parameters
    parser.add_argument(
        '-t', '--threads',
        type=int,
        default=12,
        help='线程数|Threads (default: 12)'
    )

    parser.add_argument(
        '--min-clip',
        type=int,
        default=20,
        help='最小soft-clip长度|Minimum soft-clip length (default: 20)'
    )

    parser.add_argument(
        '--min-support',
        type=int,
        default=5,
        help='最小支持reads数|Minimum supporting reads (default: 5)'
    )

    parser.add_argument(
        '--score-threshold',
        type=int,
        default=1000,
        help='得分阈值|Score threshold (default: 1000)'
    )

    # 工具路径|Tool paths
    parser.add_argument(
        '--bowtie2-path',
        default='bowtie2',
        help='Bowtie2可执行文件路径|Bowtie2 executable path'
    )

    parser.add_argument(
        '--samtools-path',
        default='samtools',
        help='samtools可执行文件路径|samtools executable path'
    )

    # FASTQ文件模式|FASTQ file patterns
    parser.add_argument(
        '--read1-suffix',
        default='_1.clean.fq.gz',
        help='R1文件后缀（包含扩展名，默认匹配fastp输出）|Read 1 file suffix with extension (default: _1.clean.fq.gz, matches fastp output)'
    )

    parser.add_argument(
        '--read2-suffix',
        default='_2.clean.fq.gz',
        help='R2文件后缀（包含扩展名，默认匹配fastp输出）|Read 2 file suffix with extension (default: _2.clean.fq.gz, matches fastp output)'
    )

    # 流程控制|Process control
    parser.add_argument(
        '--force',
        action='store_true',
        help='强制重新运行所有步骤|Force rerun all steps'
    )

    parser.add_argument(
        '--verbose',
        action='store_true',
        help='显示详细日志|Show verbose logs'
    )

    parser.add_argument(
        '--quiet',
        action='store_true',
        help='仅显示错误日志|Show error logs only'
    )

    return parser.parse_args()


def main():
    """主函数|Main function"""
    args = parse_arguments()

    try:
        # 创建配置|Create configuration
        config = InsertDetectionConfig(
            genome=args.genome,
            insert_sequence=args.insert,
            fastq_dir=args.fastq_dir,
            output_dir=args.output_dir,
            threads=args.threads,
            skip_existing=not args.force,  # force=True时skip_existing=False
            min_clip=args.min_clip,
            min_support=args.min_support,
            score_threshold=args.score_threshold,
            bowtie2_path=args.bowtie2_path,
            samtools_path=args.samtools_path,
            read1_suffix=args.read1_suffix,
            read2_suffix=args.read2_suffix
        )

        # 验证配置|Validate configuration
        config.validate()

        # 设置日志级别|Setup log level
        if args.verbose:
            log_level = "DEBUG"
        elif args.quiet:
            log_level = "ERROR"
        else:
            log_level = "INFO"

        # 创建日志目录|Create log directory
        log_dir = Path(config.output_dir) / "99_logs"
        log_dir.mkdir(parents=True, exist_ok=True)
        log_file = log_dir / "insert_detection.log"

        logger_manager = InsertDetectionLogger(str(log_file), log_level)
        logger = logger_manager.get_logger()

        # 检查依赖|Check dependencies
        if not check_dependencies(config, logger):
            logger.error("依赖检查失败，请安装所需工具|Dependency check failed, please install required tools")
            sys.exit(1)

        # 生成全局流程信息|Generate global pipeline info
        logger.info("=" * 60)
        logger.info("生成流程信息|Generating pipeline information")
        logger.info("=" * 60)
        pipeline_info_generator = PipelineInfoGenerator(logger, config)
        pipeline_info_generator.generate_pipeline_info(config.output_dir)

        # 识别样品|Identify samples
        logger.info("=" * 60)
        logger.info("识别样品信息|Identifying sample information")
        logger.info("=" * 60)
        samples = identify_samples(
            config.fastq_dir,
            config.read1_suffix,
            config.read2_suffix,
            logger
        )

        if not samples:
            logger.error("未找到有效样品|No valid samples found")
            sys.exit(1)

        # 运行检测流程|Run detection pipeline
        logger.info("=" * 60)
        logger.info("开始检测流程|Starting detection pipeline")
        logger.info("=" * 60)

        detector = InsertDetector(config, logger)
        results = detector.run(samples)

        # 为每个样本生成pipeline_info|Generate pipeline info for each sample
        logger.info("=" * 60)
        logger.info("生成样本流程信息|Generating sample pipeline information")
        logger.info("=" * 60)
        for sample_id in samples.keys():
            pipeline_info_generator.generate_pipeline_info(config.output_dir, sample_id)

        # 输出结果|Output results
        if results:
            logger.info("=" * 60)
            logger.info("输出结果|Outputting results")
            logger.info("=" * 60)
            results_dir = Path(config.output_dir) / "04_results"
            output_file = write_results(results, config.output_dir)
            write_summary(results, config.output_dir)
            logger.info(f"找到 {len(results)} 个插入位点|Found {len(results)} insertion sites")
            logger.info(f"结果保存到|Results saved to: {output_file}")
            logger.info(f"汇总文件|Summary file: {results_dir / 'summary.txt'}")
        else:
            logger.warning("未检测到插入位点|No insertion sites detected")

        logger.info("=" * 60)
        logger.info("检测完成|Detection completed")
        logger.info(f"结果保存在|Results saved to: {config.output_dir}")
        logger.info("=" * 60)
        sys.exit(0)

    except ValueError as e:
        print(f"配置错误|Configuration error: {e}")
        sys.exit(1)

    except KeyboardInterrupt:
        print("\n用户中断|User interrupted")
        sys.exit(130)

    except Exception as e:
        print(f"错误|Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
