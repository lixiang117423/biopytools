"""
EDTA转座子注释主程序模块|EDTA TE Annotation Main Module
"""

import argparse
import sys
import os
from pathlib import Path

from .config import EDTAConfig, PanEDTAConfig
from .utils import EDTALogger
from .edta_runner import EDTARunner
from .panedta_runner import PanEDTARunner


def main():
    """主函数|Main function"""
    parser = argparse.ArgumentParser(
        description='EDTA转座子注释工具|EDTA TE Annotation Tool',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    # 子命令|Subcommands
    subparsers = parser.add_subparsers(dest='mode', help='运行模式|Run mode')

    # EDTA模式|EDTA mode
    edta_parser = subparsers.add_parser('edta', help='单基因组转座子注释|Single-genome TE annotation')
    edta_parser.add_argument('-i', '--genome', required=True,
                           help='基因组FASTA文件|Genome FASTA file')
    edta_parser.add_argument('--species', default='others', choices=['Rice', 'Maize', 'others'],
                           help='物种类型|Species type')
    edta_parser.add_argument('--step', default='all', choices=['all', 'filter', 'final', 'anno'],
                           help='运行步骤|Step to run')
    edta_parser.add_argument('--overwrite', type=int, default=0, choices=[0, 1],
                           help='覆盖已有结果|Overwrite existing results')
    edta_parser.add_argument('--cds',
                           help='CDS序列文件|CDS sequences file')
    edta_parser.add_argument('--curatedlib',
                           help='筛选TE库|Curated TE library')
    edta_parser.add_argument('--rmlib',
                           help='RepeatModeler库|RepeatModeler library')
    edta_parser.add_argument('--sensitive', type=int, default=0, choices=[0, 1],
                           help='使用RepeatModeler|Use RepeatModeler')
    edta_parser.add_argument('--anno', type=int, default=0, choices=[0, 1],
                           help='执行全基因组注释|Perform whole-genome annotation')
    edta_parser.add_argument('--rmout',
                           help='RepeatMasker输出文件|RepeatMasker output file')
    edta_parser.add_argument('--maxdiv', type=int, default=40,
                           help='最大分歧度|Maximum divergence')
    edta_parser.add_argument('--evaluate', type=int, default=0, choices=[0, 1],
                           help='评估注释一致性|Evaluate annotation consistency')
    edta_parser.add_argument('--exclude',
                           help='排除区域BED文件|Exclude regions BED file')
    edta_parser.add_argument('--force', type=int, default=0, choices=[0, 1],
                           help='强制使用水稻TE|Force to use rice TEs')
    edta_parser.add_argument('--u', type=float, default=1.3e-8,
                           help='中性突变率|Neutral mutation rate')
    edta_parser.add_argument('-t', '--threads', type=int, default=12,
                           help='线程数|Number of threads')
    edta_parser.add_argument('--debug', type=int, default=0, choices=[0, 1],
                           help='保留中间文件|Retain intermediate files')
    edta_parser.add_argument('-o', '--output-dir', default='./edta_output',
                           help='输出目录|Output directory')
    edta_parser.add_argument('--edta-path',
                           help='EDTA安装路径|EDTA installation path')

    # panEDTA模式|panEDTA mode
    panedta_parser = subparsers.add_parser('panedta', help='泛基因组转座子注释|Pan-genome TE annotation')
    panedta_parser.add_argument('-i', '--genome-list', required=True,
                              help='基因组列表文件|Genome list file')
    panedta_parser.add_argument('-c', '--cds',
                              help='CDS序列文件|CDS sequences file')
    panedta_parser.add_argument('-l', '--curatedlib',
                              help='筛选TE库|Curated TE library')
    panedta_parser.add_argument('-f', '--fl-copy', type=int, default=3,
                              help='全长拷贝数阈值|Full-length copy number cutoff')
    panedta_parser.add_argument('-a', '--anno', type=int, default=1, choices=[0, 1],
                              help='执行全基因组注释|Perform whole-genome annotation')
    panedta_parser.add_argument('--overwrite', type=int, default=0, choices=[0, 1],
                              help='覆盖已有结果|Overwrite existing results')
    panedta_parser.add_argument('-t', '--threads', type=int, default=12,
                              help='线程数|Number of threads')
    panedta_parser.add_argument('-o', '--output-dir', default='./panedta_output',
                              help='输出目录|Output directory')
    panedta_parser.add_argument('--edta-path',
                              help='EDTA安装路径|EDTA installation path')

    args = parser.parse_args()

    # 检查模式|Check mode
    if not args.mode:
        parser.print_help()
        sys.exit(1)

    try:
        if args.mode == 'edta':
            # 创建EDTA配置|Create EDTA config
            config = EDTAConfig(
                genome=args.genome,
                species=args.species,
                step=args.step,
                overwrite=args.overwrite,
                cds=args.cds,
                curatedlib=args.curatedlib,
                rmlib=args.rmlib,
                sensitive=args.sensitive,
                anno=args.anno,
                rmout=args.rmout,
                maxdiv=args.maxdiv,
                evaluate=args.evaluate,
                exclude=args.exclude,
                force=args.force,
                u=args.u,
                threads=args.threads,
                debug=args.debug,
                output_dir=args.output_dir,
                edta_path=args.edta_path
            )

            # 验证配置|Validate config
            config.validate()

            # 创建运行器并运行|Create runner and run
            runner = EDTARunner(config)
            success = runner.run()

            sys.exit(0 if success else 1)

        elif args.mode == 'panedta':
            # 创建PanEDTA配置|Create PanEDTA config
            config = PanEDTAConfig(
                genome_list=args.genome_list,
                cds=args.cds,
                curatedlib=args.curatedlib,
                fl_copy=args.fl_copy,
                anno=args.anno,
                overwrite=args.overwrite,
                threads=args.threads,
                output_dir=args.output_dir,
                edta_path=args.edta_path
            )

            # 验证配置|Validate config
            config.validate()

            # 创建运行器并运行|Create runner and run
            runner = PanEDTARunner(config)
            success = runner.run()

            sys.exit(0 if success else 1)

    except KeyboardInterrupt:
        print("\n用户中断|User interrupted")
        sys.exit(1)
    except Exception as e:
        print(f"错误|Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
