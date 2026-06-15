"""
Cactus命令行接口|Cactus Command Line Interface
"""

import sys
import argparse
from pathlib import Path

from .config import CactusConfig
from .pangenome import CactusPangenomeRunner


def parse_arguments():
    """
    解析命令行参数|Parse command line arguments

    Returns:
        参数对象|Arguments object
    """
    parser = argparse.ArgumentParser(
        description='Cactus泛基因组批量分析工具|Cactus Pangenome Batch Analysis Tool',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例|Examples:
  # 基本用法|Basic usage
  %(prog)s -s seqfile.txt -o output/ -r ref_genome

  # 指定输出格式|Specify output formats
  %(prog)s -s seqfile.txt -o output/ -r ref_genome --formats gfa gbz odgi

  # 自定义Singularity路径|Custom Singularity path
  %(prog)s -s seqfile.txt -o output/ -r ref_genome --singularity ~/bin/singularity --cactus-sif ~/images/cactus.sif

  # 保留jobstore用于调试|Keep jobstore for debugging
  %(prog)s -s seqfile.txt -o output/ -r ref_genome --no-cleanup

序列文件格式|Sequence file format (seqfile.txt):
  sample_name1  /path/to/genome1.fa
  sample_name2  /path/to/genome2.fa
  sample_name3  /path/to/genome3.fa
  ...

  两列格式：样本名 + FASTA路径（制表符或空格分隔）|Two columns: sample_name + FASTA path (tab or space separated)
  参考基因组名称必须使用-r指定|Reference genome name must be specified with -r
        '''
    )

    # 必需参数|Required arguments
    parser.add_argument('-s', '--seqfile',
                        required=True,
                        help='序列文件|Sequence file (两列格式：样本名 + 路径|Two columns: sample_name + path)')
    parser.add_argument('-o', '--output',
                        required=True,
                        help='输出目录|Output directory')
    parser.add_argument('-r', '--reference',
                        required=True,
                        help='参考基因组名称|Reference genome name (必须与seqfile第一个基因组名称匹配|Must match first genome in seqfile)')

    # Singularity相关|Singularity related
    parser.add_argument('--singularity',
                        default='~/miniforge3/envs/singularity_v.3.8.7/bin/singularity',
                        help='Singularity可执行文件路径|Singularity executable path (default: ~/miniforge3/envs/singularity_v.3.8.7/bin/singularity)')
    parser.add_argument('--cactus-sif',
                        default='~/software/singularity/cactus_v3.1.4.sif',
                        help='Cactus SIF镜像路径|Cactus SIF image path (default: ~/software/singularity/cactus_v3.1.4.sif)')

    # 流程参数|Pipeline parameters
    parser.add_argument('--jobstore',
                        default='cactus-jobstore',
                        help='Toil jobstore目录名称|Toil jobstore directory name (default: cactus-jobstore)')
    parser.add_argument('--out-name',
                        default='cactus_output',
                        help='输出文件前缀|Output file prefix (default: cactus_output)')
    parser.add_argument('--no-cleanup',
                        action='store_true',
                        help='保留jobstore不删除|Keep jobstore without deleting (useful for debugging)')

    # 输出格式|Output formats
    # 注意：HAL是默认输出（{out_name}.full.hal），不需要在列表中指定
    # Note: HAL is default output ({out_name}.full.hal), no need to specify in list
    parser.add_argument('--formats',
                        nargs='+',
                        choices=['gfa', 'gbz', 'odgi', 'hal', 'vg', 'vcf', 'xg'],
                        default=['gfa', 'gbz'],
                        help='输出格式|Output formats (default: gfa gbz)')

    # 性能参数|Performance parameters
    parser.add_argument('-t', '--threads',
                        type=int,
                        default=12,
                        help='CPU核心数|Number of CPU cores (default: 12)')
    parser.add_argument('-m', '--max-memory',
                        default='100G',
                        help='最大内存|Maximum memory (default: 100G)')

    # 绑定目录|Bind directories
    parser.add_argument('--bind',
                        action='append',
                        help='绑定目录到容器|Bind directory to container (可多次使用|can be used multiple times)')

    # 工作目录|Work directory
    parser.add_argument('--workDir',
                        help='工作目录（用于临时文件）|Work directory for temporary files (default: output_dir/work)')

    # 日志选项|Logging options
    parser.add_argument('--log-level',
                        choices=['DEBUG', 'INFO', 'WARN', 'ERROR'],
                        default='INFO',
                        help='日志级别|Log level (default: INFO)')
    parser.add_argument('--quiet',
                        action='store_true',
                        help='静默模式|Quiet mode')

    return parser.parse_args()


def main():
    """主函数|Main function"""
    try:
        # 解析参数|Parse arguments
        args = parse_arguments()

        # 创建输出目录|Create output directory
        output_dir = Path(args.output).expanduser()
        output_dir.mkdir(parents=True, exist_ok=True)

        # 设置日志|Setup logging
        from .utils import CactusLogger
        log_file = output_dir / f"{args.out_name}.log"
        logger_manager = CactusLogger(log_file, args.log_level)
        logger = logger_manager.get_logger()

        # 验证seqfile存在|Validate seqfile existence
        seqfile = Path(args.seqfile).expanduser()
        if not seqfile.exists():
            logger.error(f"序列文件不存在|Sequence file does not exist: {seqfile}")
            sys.exit(1)

        # 创建配置|Create config
        config = CactusConfig(
            seqfile=str(seqfile),
            output_dir=str(output_dir),
            reference=args.reference,
            singularity_path=args.singularity,
            cactus_sif=args.cactus_sif,
            jobstore=args.jobstore,
            out_name=args.out_name,
            cleanup=not args.no_cleanup,
            output_formats=args.formats,
            threads=args.threads,
            max_memory=args.max_memory,
            bind_paths=args.bind,
            work_dir=args.workDir,
            log_level=args.log_level,
            verbose=not args.quiet
        )

        # 验证配置|Validate config
        try:
            config.validate()
        except ValueError as e:
            logger.error(f"配置错误|Configuration error: {e}")
            sys.exit(1)

        # 创建分析器|Create analyzer
        analyzer = CactusPangenomeRunner(config)

        # 运行分析|Run analysis
        success = analyzer.run()

        if success:
            sys.exit(0)
        else:
            logger.error("分析失败|Analysis failed")
            sys.exit(1)

    except KeyboardInterrupt:
        sys.stderr.write("用户中断操作|User interrupted operation\n")
        sys.exit(1)
    except Exception as e:
        sys.stderr.write(f"程序错误|Program error: {e}\n")
        sys.exit(1)


if __name__ == '__main__':
    main()
