"""
YaHS主程序模块|YaHS Main Module

Hi-C scaffolding流程的命令行入口
Command-line entry point for Hi-C scaffolding pipeline
"""

import argparse
import sys
from pathlib import Path

from .config import YaHSConfig
from .pipeline import YaHSPipeline


def parse_arguments():
    """
    解析命令行参数|Parse command line arguments

    Returns:
        解析后的参数|Parsed arguments
    """
    parser = argparse.ArgumentParser(
        description='YaHS Hi-C Scaffolding Pipeline|YaHS Hi-C scaffolding pipeline',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例|Examples:
  # 完整流程|Complete pipeline
  %(prog)s -r genome.fa -1 hic_R1.fq.gz -2 hic_R2.fq.gz

  # 自定义参数|Custom parameters
  %(prog)s -r genome.fa -1 hic_R1.fq.gz -2 hic_R2.fq.gz -e GATC -t 24

  # 运行单个步骤|Run single step
  %(prog)s -r genome.fa -1 hic_R1.fq.gz -2 hic_R2.fq.gz -s 3
        '''
    )

    # 必需参数|Required parameters
    required = parser.add_argument_group('必需参数|Required parameters')
    required.add_argument(
        '-r', '--ref',
        required=True,
        dest='ref_fa',
        help='参考基因组FASTA文件|Reference genome FASTA file'
    )
    required.add_argument(
        '-1', '--hic-r1',
        required=True,
        dest='hic_r1',
        help='Hi-C R1测序文件|Hi-C R1 sequencing file'
    )
    required.add_argument(
        '-2', '--hic-r2',
        required=True,
        dest='hic_r2',
        help='Hi-C R2测序文件|Hi-C R2 sequencing file'
    )

    # 输出配置|Output configuration
    output = parser.add_argument_group('输出配置|Output configuration')
    output.add_argument(
        '-o', '--output-dir',
        default='./yahs_output',
        dest='output_dir',
        help='输出目录|Output directory (default: ./yahs_output)'
    )

    # 资源配置|Resource configuration
    resources = parser.add_argument_group('资源配置|Resource configuration')
    resources.add_argument(
        '-t', '--threads',
        type=int,
        default=12,
        help='线程数|Number of threads (default: 12)'
    )
    resources.add_argument(
        '--java-ram',
        default='32G',
        help='Java内存|Java memory (default: 32G)'
    )
    resources.add_argument(
        '--sam-ram',
        default='4G',
        help='Samtools排序内存|Samtools sort memory (default: 4G)'
    )

    # YaHS 核心参数|YaHS core parameters
    yahs_params = parser.add_argument_group('YaHS参数|YaHS parameters')
    yahs_params.add_argument(
        '-e', '--enzyme',
        default='GATC',
        dest='enzyme_seq',
        help='限制性酶切位点|Restriction enzyme sequence (default: GATC)'
    )
    yahs_params.add_argument(
        '--min-len',
        type=int,
        default=10000,
        dest='min_len',
        help='最小contig长度|Minimum contig length (default: 10000)'
    )
    yahs_params.add_argument(
        '--min-mapq',
        type=int,
        default=30,
        dest='min_mapq',
        help='最小MAPQ值|Minimum MAPQ value (default: 30)'
    )
    yahs_params.add_argument(
        '--no-contig-ec',
        action='store_true',
        help='跳过contig错误校正|Skip contig error correction'
    )
    yahs_params.add_argument(
        '--no-scaffold-ec',
        action='store_true',
        help='跳过scaffold错误校正|Skip scaffold error correction'
    )
    yahs_params.add_argument(
        '--resolutions',
        help='分辨率列表(逗号分隔)|Resolution list (comma-separated)'
    )
    yahs_params.add_argument(
        '--rounds',
        type=int,
        default=1,
        dest='rounds_per_resolution',
        help='每分辨率运行轮数|Rounds per resolution (default: 1)'
    )
    yahs_params.add_argument(
        '--telo-motif',
        dest='telo_motif',
        help='端粒序列模体|Telomeric sequence motif'
    )

    # 工具路径|Tool paths
    tools = parser.add_argument_group('工具路径|Tool paths')
    tools.add_argument(
        '--yahs-bin',
        default='~/miniforge3/envs/yahs_v.1.2.2/bin/yahs',
        help='YaHS可执行文件路径|YaHS executable path (default: ~/miniforge3/envs/yahs_v.1.2.2/bin/yahs)'
    )
    tools.add_argument(
        '--juicer-bin',
        default='~/miniforge3/envs/yahs_v.1.2.2/bin/juicer',
        help='juicer可执行文件路径|juicer executable path (default: ~/miniforge3/envs/yahs_v.1.2.2/bin/juicer)'
    )
    tools.add_argument(
        '--juicer-jar',
        default='~/software/juicer/scripts/juicer_tools.jar',
        help='juicer_tools.jar文件路径|juicer_tools.jar file path (default: ~/software/juicer/scripts/juicer_tools.jar)'
    )
    tools.add_argument(
        '--bwa-bin',
        default='~/miniforge3/envs/Population_genetics/bin/bwa',
        help='BWA可执行文件路径|BWA executable path (default: ~/miniforge3/envs/Population_genetics/bin/bwa)'
    )
    tools.add_argument(
        '--samtools-bin',
        default='~/miniforge3/envs/GATK_v.4.6.2.0/bin/samtools',
        help='samtools可执行文件路径|samtools executable path (default: ~/miniforge3/envs/GATK_v.4.6.2.0/bin/samtools)'
    )
    tools.add_argument(
        '--java-cmd',
        default='java',
        help='Java可执行文件路径|Java executable path (default: java)'
    )

    # 执行控制|Execution control
    execution = parser.add_argument_group('执行控制|Execution control')
    execution.add_argument(
        '-s', '--step',
        choices=['1', '2', '3', '4', '5', '6'],
        help='运行指定步骤(1-6)|Run specified step (1-6)'
    )
    execution.add_argument(
        '--force-rerun',
        action='store_true',
        help='强制重新运行所有步骤|Force rerun all steps'
    )
    execution.add_argument(
        '--keep-temp',
        action='store_true',
        help='保留临时文件|Keep temporary files'
    )

    return parser.parse_args()


def main():
    """
    主函数|Main function

    命令行入口点|Command-line entry point
    """
    # 解析参数|Parse arguments
    args = parse_arguments()

    try:
        # 创建配置|Create configuration
        config = YaHSConfig(
            ref_fa=args.ref_fa,
            hic_r1=args.hic_r1,
            hic_r2=args.hic_r2,
            output_dir=args.output_dir,
            threads=args.threads,
            java_ram=args.java_ram,
            sam_ram=args.sam_ram,
            enzyme_seq=args.enzyme_seq,
            min_len=args.min_len,
            min_mapq=args.min_mapq,
            no_contig_ec=args.no_contig_ec,
            no_scaffold_ec=args.no_scaffold_ec,
            resolutions=args.resolutions,
            rounds_per_resolution=args.rounds_per_resolution,
            telo_motif=args.telo_motif,
            yahs_bin=args.yahs_bin,
            juicer_bin=args.juicer_bin,
            juicer_jar=args.juicer_jar,
            bwa_bin=args.bwa_bin,
            samtools_bin=args.samtools_bin,
            java_cmd=args.java_cmd,
            step=args.step,
            force_rerun=args.force_rerun,
            keep_temp=args.keep_temp
        )

        # 运行流程|Run pipeline
        pipeline = YaHSPipeline(config)
        exit_code = pipeline.run()

        sys.exit(exit_code)

    except ValueError as e:
        print(f"配置错误|Configuration error: {e}", file=sys.stderr)
        sys.exit(1)

    except KeyboardInterrupt:
        print("\n流程被用户中断|Pipeline interrupted by user", file=sys.stderr)
        sys.exit(130)

    except Exception as e:
        print(f"未知错误|Unknown error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
