"""
CPhasing主程序模块|CPhasing Main Module

命令行入口，解析参数并执行CPhasing命令
Command line entry point, parse arguments and execute CPhasing commands
"""

import argparse
import sys
from pathlib import Path

from .config import CPhasingConfig
from .runner import CPhasingRunner


# CPhasing所有可用子命令|All available CPhasing subcommands
CPHASING_SUBCOMMANDS = [
    'pipeline', 'mapper', 'alleles', 'hypergraph', 'hg',
    'hyperpartition', 'hp', 'partition',
    'hyperoptimize', 'ho', 'optimize',
    'scaffolding', 'sf', 'scaf', 'scaffold',
    'chimeric', 'hcr', 'prepare', 'summary',
    'evaluator', 'eval',
    'pairs2contacts', 'pairs2contact',
    'pairs2cool', 'pair2cool', 'p2c',
    'plot', 'rename', 'sort-chromosomes', 'sort-chrom', 'sort-chroms',
    'chromsizes', 'contigsizes',
    'collapse', 'agp2fasta', 'prepartition', 'gfa2depth',
]


def parse_arguments():
    """
    解析命令行参数|Parse command line arguments

    Returns:
        解析后的参数|Parsed arguments
    """
    parser = argparse.ArgumentParser(
        description="CPhasing基因组分相和挂载工具|CPhasing genome phasing and scaffolding tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例|Examples:
  %(prog)s -i genome.fa --hic1 R1.fq.gz --hic2 R2.fq.gz -t 12 -n 16:2
        """
    )

    # 子命令（位置参数，可选，默认pipeline）
    # Subcommand (positional, optional, default pipeline)
    parser.add_argument(
        'subcommand',
        nargs='?',
        default='pipeline',
        help=f'CPhasing子命令|CPhasing subcommand (默认|default: pipeline)'
    )

    # --- 透传参数分隔符 ---
    # Pass-through argument separator
    parser.add_argument(
        '--',
        dest='extra_args_start',
        action='store_true',
        help=argparse.SUPPRESS,
    )

    # --- pipeline 专用参数 ---
    # pipeline-specific parameters
    parser.add_argument(
        '-i', '--input', dest='fasta',
        help='基因组FASTA文件|Genome FASTA file (pipeline模式必需|required for pipeline)'
    )
    parser.add_argument(
        '--hic1',
        help='Hi-C R1 reads文件|Hi-C R1 reads file (pipeline模式必需|required for pipeline)'
    )
    parser.add_argument(
        '--hic2',
        help='Hi-C R2 reads文件|Hi-C R2 reads file (pipeline模式必需|required for pipeline)'
    )
    parser.add_argument(
        '-n', '--groups',
        default='0',
        help='分组数|Number of groups (例如: "8:4" | e.g., "8:4", "0"=自动|auto)'
    )
    parser.add_argument(
        '-t', '--threads',
        type=int,
        default=12,
        help='线程数|Number of threads (默认|default: 12)'
    )
    parser.add_argument(
        '--mode',
        choices=['phasing', 'haploid', 'basal', 'basal_withprune'],
        default='phasing',
        help='分相模式|Phasing mode (默认|default: phasing)'
    )
    parser.add_argument(
        '--preset',
        choices=['precision', 'sensitive', 'very-sensitive', 'nofilter'],
        default='precision',
        help='分析预设|Analysis preset (默认|default: precision)'
    )

    # --- 通用参数 ---
    # Generic parameters
    parser.add_argument(
        '-o', '--output-dir',
        default='./cphasing_output',
        help='输出目录|Output directory (默认|default: ./cphasing_output)'
    )
    parser.add_argument(
        '--steps',
        help='运行指定步骤|Run specified steps (pipeline模式)'
    )
    parser.add_argument(
        '--skip-steps',
        help='跳过步骤|Skip steps (pipeline模式)'
    )
    parser.add_argument(
        '--hic-aligner',
        choices=['_chromap', 'chromap', 'minimap2', 'bwa-mem2'],
        default='_chromap',
        help='Hi-C比对器|Hi-C aligner (默认|default: _chromap)'
    )
    parser.add_argument(
        '--hic-mapper-k',
        type=int,
        help='Hi-C mapper kmer大小|Hi-C mapper kmer size'
    )
    parser.add_argument(
        '--hic-mapper-w',
        type=int,
        help='Hi-C mapper窗口大小|Hi-C mapper window size'
    )
    parser.add_argument(
        '--mapping-quality',
        type=int,
        default=0,
        help='最小比对质量|Minimum mapping quality (默认|default: 0)'
    )
    parser.add_argument(
        '--hcr',
        action='store_true',
        help='启用高置信区域|Enable high confidence regions'
    )
    parser.add_argument(
        '--pattern',
        help='酶切位点模式|Restriction enzyme pattern (例如: AAGCTT)'
    )
    parser.add_argument(
        '--low-memory',
        action='store_true',
        help='低内存模式|Low memory mode'
    )

    return parser.parse_args()


def main():
    """
    主函数|Main function
    """
    args = parse_arguments()

    try:
        # 收集透传参数|Collect pass-through arguments
        # sys.argv中在'--'之后的所有参数都透传给CPhasing
        # All arguments after '--' in sys.argv are passed through to CPhasing
        extra_args = []
        passthrough = False
        for arg in sys.argv[1:]:
            if arg == '--':
                passthrough = True
                continue
            if passthrough:
                extra_args.append(arg)

        config = CPhasingConfig(
            subcommand=args.subcommand,
            fasta=args.fasta,
            hic1=args.hic1,
            hic2=args.hic2,
            groups=args.groups,
            threads=args.threads,
            mode=args.mode,
            preset=args.preset,
            output_dir=args.output_dir,
            steps=args.steps,
            skip_steps=args.skip_steps,
            hic_aligner=args.hic_aligner,
            hic_mapper_k=args.hic_mapper_k,
            hic_mapper_w=args.hic_mapper_w,
            mapping_quality=args.mapping_quality,
            hcr=args.hcr,
            pattern=args.pattern,
            low_memory=args.low_memory,
            extra_args=extra_args,
        )

        config.validate()

        runner = CPhasingRunner(config)
        success = runner.run()

        sys.exit(0 if success else 1)

    except ValueError as e:
        print(f"配置错误|Configuration error: {e}", file=sys.stderr)
        sys.exit(1)
    except KeyboardInterrupt:
        print("\n用户中断|User interrupted", file=sys.stderr)
        sys.exit(130)
    except Exception as e:
        print(f"未知错误|Unknown error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
