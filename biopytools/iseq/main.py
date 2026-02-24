"""
iSeq下载主程序模块|iSeq Download Main Module
"""

import argparse
import sys
from pathlib import Path
from .config import ISeqConfig
from .utils import ISeqLogger, CondaCommandRunner
from .calculator import ISeqCalculator


class ISeqDownloader:
    """iSeq下载主类|Main iSeq Downloader Class"""

    def __init__(self, **kwargs):
        # 初始化配置|Initialize configuration
        self.config = ISeqConfig(**kwargs)
        self.config.validate()

        # 初始化日志|Initialize logging
        self.logger_manager = ISeqLogger(self.config.output_path)
        self.logger = self.logger_manager.get_logger()

        # 初始化命令执行器|Initialize command runner
        self.cmd_runner = CondaCommandRunner(
            self.logger,
            self.config.conda_env,
            self.config.output_path
        )

        # 初始化计算器|Initialize calculator
        self.calculator = ISeqCalculator(
            self.config,
            self.logger,
            self.cmd_runner
        )

    def run_analysis(self):
        """运行分析|Run analysis"""
        try:
            success = self.calculator.run_analysis()
            if not success:
                self.logger.error("下载失败|Download failed")
                sys.exit(1)
        except Exception as e:
            self.logger.error(f"程序执行出错|Program execution error: {str(e)}")
            import traceback
            self.logger.debug(traceback.format_exc())
            sys.exit(1)


def main():
    """主函数|Main function"""
    parser = argparse.ArgumentParser(
        description='基于iSeq的公共测序数据下载工具|Public sequencing data download tool based on iSeq',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # 必需参数|Required arguments
    parser.add_argument('-i', '--accession',
                       required=True,
                       help='项目/样本/实验ID|Project/Sample/Experiment ID (e.g., PRJNA1014406)')

    # 可选参数|Optional arguments
    parser.add_argument('-o', '--output-dir',
                       default='./iseq_output',
                       help='输出目录|Output directory')

    parser.add_argument('-a', '--iseq-path',
                       default='/share/org/YZWL/yzwl_lixg/miniforge3/envs/iseq_v.1.9.8/bin/iseq',
                       help='iSeq软件路径|iSeq software path')

    parser.add_argument('-c', '--conda-env',
                       default='iseq_v.1.9.8',
                       help='Conda环境名|Conda environment name')

    # 下载选项|Download options
    parser.add_argument('-m', '--metadata-only',
                       action='store_true',
                       help='仅下载元数据|Only download metadata')

    parser.add_argument('-g', '--gzip',
                       action='store_true',
                       help='下载gzip格式FASTQ|Download FASTQ in gzip format')

    parser.add_argument('-q', '--fastq',
                       action='store_true',
                       help='转换为FASTQ格式|Convert to FASTQ format')

    parser.add_argument('-e', '--merge',
                       choices=['ex', 'sa', 'st'],
                       help='合并选项|Merge option (ex/sa/st)')

    # 性能参数|Performance parameters
    parser.add_argument('-t', '--threads',
                       type=int, default=16,
                       help='线程数|Number of threads')

    parser.add_argument('-p', '--parallel',
                       type=int, default=10,
                       help='并行连接数|Number of parallel connections')

    parser.add_argument('-s', '--speed',
                       type=int,
                       help='下载速度限制(MB/s)|Download speed limit (MB/s)')

    # 数据库选项|Database options
    parser.add_argument('-d', '--database',
                       default='ena',
                       choices=['ena', 'sra'],
                       help='数据库选择|Database selection')

    parser.add_argument('--protocol',
                       default='ftp',
                       choices=['ftp', 'https'],
                       help='协议选择|Protocol selection')

    # 高级选项|Advanced options
    parser.add_argument('--use-aspera',
                       action='store_true',
                       help='使用Aspera下载|Use Aspera for download')

    parser.add_argument('--skip-md5',
                       action='store_true',
                       help='跳过MD5校验|Skip MD5 check')

    parser.add_argument('--quiet',
                       action='store_true',
                       help='静默模式|Quiet mode (suppress progress bars)')

    args = parser.parse_args()

    # 创建下载器并运行|Create downloader and run
    downloader = ISeqDownloader(
        accession=args.accession,
        output_dir=args.output_dir,
        iseq_path=args.iseq_path,
        conda_env=args.conda_env,
        metadata_only=args.metadata_only,
        gzip=args.gzip,
        fastq=args.fastq,
        merge=args.merge,
        threads=args.threads,
        parallel=args.parallel,
        speed=args.speed,
        database=args.database,
        protocol=args.protocol,
        use_aspera=args.use_aspera,
        skip_md5=args.skip_md5,
        quiet=args.quiet
    )

    downloader.run_analysis()


if __name__ == "__main__":
    main()
