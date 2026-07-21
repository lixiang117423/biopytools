"""
HiTE 单基因组转座子检测与注释|HiTE single-genome TE detection and annotation

命令行入口:通过 singularity 直接挂载调用 HiTE
CLI entry: calls HiTE via singularity direct-mount
"""

import argparse
import sys

from .config import HiteConfig
from .logger import HiteLoggerManager
from .hite_runner import HiteRunner
from .results import HiteResultsProcessor


def parse_arguments():
    """解析命令行参数|Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="HiTE转座子检测与注释|HiTE TE detection and annotation"
    )
    parser.add_argument(
        '-i', '--input', dest='genome', required=True,
        help='基因组FASTA文件|Genome FASTA file',
    )
    parser.add_argument(
        '-o', '--output-dir', dest='output_dir', default='./hite_output',
        help='输出目录|Output directory (default: ./hite_output)',
    )
    parser.add_argument(
        '-t', '--threads', type=int, default=12,
        help='线程数|Number of threads (default: 12)',
    )
    parser.add_argument(
        '--singularity-path', dest='singularity_path', default=None,
        help='Singularity可执行文件路径|Singularity executable path',
    )
    parser.add_argument(
        '--sif-file', dest='sif_file', default=None,
        help='HiTE SIF镜像路径|HiTE SIF image path',
    )
    parser.add_argument(
        '--plant', type=int, default=1, choices=[0, 1],
        help='是否植物基因组(1/0)|Plant genome (1/0, default: 1)',
    )
    parser.add_argument(
        '--annotate', type=int, default=0, choices=[0, 1],
        help='是否注释基因组(1/0)|Annotate genome (1/0, default: 0)',
    )
    parser.add_argument(
        '--recover', type=int, default=0, choices=[0, 1],
        help='HiTE断点续跑(1/0)|HiTE recovery mode (1/0, default: 0)',
    )
    parser.add_argument(
        '--domain', type=int, default=0, choices=[0, 1],
        help='预测TE蛋白结构域(1/0)|Predict TE domains (1/0, default: 0)',
    )
    parser.add_argument(
        '--te-type', dest='te_type', default='all',
        choices=['ltr', 'tir', 'helitron', 'non-ltr', 'all'],
        help='TE类型|TE type (default: all)',
    )
    parser.add_argument('--chunk-size', dest='chunk_size', type=int, default=400,
                        help='基因组分块MB|Genome chunk size MB (default: 400)')
    parser.add_argument('--miu', type=float, default=1.3e-8,
                        help='中性突变率|Neutral mutation rate (default: 1.3e-8)')
    parser.add_argument('--min-te-len', dest='min_te_len', type=int, default=80,
                        help='最小TE长度bp|Min TE length bp (default: 80)')
    parser.add_argument(
        '--remove-nested', dest='remove_nested', type=int, default=1, choices=[0, 1],
        help='移除嵌套TE(1/0)|Remove nested TE (1/0, default: 1)',
    )
    parser.add_argument('--curated-lib', dest='curated_lib', default=None,
                        help='可信curated TE库|Trusted curated TE library')
    parser.add_argument('--debug', type=int, default=0, choices=[0, 1],
                        help='HiTE debug模式(保留临时文件)|Debug mode (1/0)')
    return parser.parse_args()


def main() -> int:
    """
    主函数|Main function

    装配 config → logger → runner → results,执行完整流程
    Returns:
        0 成功,1 失败|0 success, 1 failure
    """
    args = parse_arguments()

    try:
        # 构建 config(只传非 None 的可选路径)|Build config
        config_kwargs = {
            'genome': args.genome,
            'output_dir': args.output_dir,
            'threads': args.threads,
            'plant': bool(args.plant),
            'annotate': bool(args.annotate),
            'recover': bool(args.recover),
            'domain': bool(args.domain),
            'te_type': args.te_type,
            'chunk_size': args.chunk_size,
            'miu': args.miu,
            'min_te_len': args.min_te_len,
            'remove_nested': bool(args.remove_nested),
            'curated_lib': args.curated_lib,
            'debug': bool(args.debug),
        }
        if args.singularity_path:
            config_kwargs['singularity_path'] = args.singularity_path
        if args.sif_file:
            config_kwargs['sif_file'] = args.sif_file

        config = HiteConfig(**config_kwargs)
        config.validate()

        # 日志写入 99_logs/|Logger writes to 99_logs/
        logger_manager = HiteLoggerManager(
            config.logs_dir, log_prefix="hite", log_level="INFO",
        )
        logger = logger_manager.get_logger()

        logger.info("=" * 80)
        logger.info("HiTE 模块启动|HiTE module started")
        logger.info(f"基因组|Genome: {config.genome}")
        logger.info(f"输出目录|Output: {config.output_dir}")
        logger.info(f"线程|Threads: {config.threads}")
        logger.info(f"SIF|SIF: {config.sif_file}")
        logger.info("=" * 80)

        # 运行|Run
        runner = HiteRunner(config, logger)
        runner.run()

        # 结果处理|Process results
        results = HiteResultsProcessor(config, logger)
        results.process_results()
        results.write_software_versions()
        results.generate_summary()
        results.print_summary()

        logger.info("HiTE 模块完成|HiTE module finished")
        return 0

    except Exception as e:
        print(f"错误|Error: {e}", file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
