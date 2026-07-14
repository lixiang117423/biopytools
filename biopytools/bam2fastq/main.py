"""
BAM to FASTQ转换主程序模块|BAM to FASTQ Conversion Main Module
"""

import argparse
import sys
from .config import BAM2FASTQConfig
from .utils import BAM2FASTQLogger, Bam2FQChecker
from .converter import BAMConverter


def parse_arguments():
    """解析命令行参数|Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='批量将BAM文件转换为FASTQ格式|Batch convert BAM files to FASTQ format',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例|Examples:
  %(prog)s -i ./sample.bam -o ./sample.fq.gz -t 32
  %(prog)s -i ./bam_dir -o ./fastq_dir -j 4 -t 16
        '''
    )

    parser.add_argument('-i', '--input',
                       required=True,
                       help='输入BAM文件或文件夹路径(包含BAM文件)|Input BAM file or directory path (containing BAM files)')

    parser.add_argument('-o', '--output-dir',
                       required=True,
                       help='输出路径(文件或目录,自动识别)|Output path (file or directory, auto-detected)')

    parser.add_argument('-t', '--threads',
                       type=int,
                       default=64,
                       help='每个BAM文件转换使用的线程数 (默认: 64)|Threads per BAM file conversion (default: 64)')

    parser.add_argument('-j', '--jobs',
                       type=int,
                       default=1,
                       help='并行处理的BAM文件数量 (默认: 1)|Number of parallel BAM file processing (default: 1)')

    parser.add_argument('--bam2fastq-path',
                       default='bam2fastq',
                       help='bam2fastq可执行文件路径 (默认: bam2fastq)|bam2fastq executable path (default: bam2fastq)')

    return parser.parse_args()


def main():
    """主函数|Main function"""
    args = parse_arguments()

    try:
        # 创建配置|Create configuration
        config = BAM2FASTQConfig(
            input_dir=args.input,
            output_dir=args.output_dir,
            threads=args.threads,
            jobs=args.jobs,
            bam2fastq_path=args.bam2fastq_path
        )
        config.validate()

        # 初始化日志|Initialize logging
        logger_manager = BAM2FASTQLogger(config.output_path)
        logger = logger_manager.get_logger()

        logger.info("=" * 50)
        logger.info("BAM to FASTQ转换工具|BAM to FASTQ Conversion Tool")
        logger.info("=" * 50)

        # 检查bam2fastq|Check bam2fastq
        checker = Bam2FQChecker(logger, config.bam2fastq_path)
        if not checker.check_bam2fastq():
            sys.exit(1)

        # 创建转换器并执行转换|Create converter and execute conversion
        converter = BAMConverter(config, logger)
        results = converter.convert_all_bams()

        # 根据结果退出程序|Exit program based on results
        if results['failed'] > 0:
            logger.warning(f"部分文件转换失败|Some files failed to convert")
            sys.exit(1)
        else:
            logger.info("所有文件转换成功|All files converted successfully")
            sys.exit(0)

    except ValueError as e:
        print(f"配置错误|Configuration error: {e}")
        sys.exit(1)

    except Exception as e:
        print(f"程序执行出错|Program execution error: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    main()
