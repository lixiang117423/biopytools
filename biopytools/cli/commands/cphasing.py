"""
CPhasing基因组分相和挂载工具|CPhasing Genome Phasing and Scaffolding Tool

支持所有CPhasing子命令：
pipeline, mapper, alleles, hypergraph, hyperpartition, plot, chimeric, collapse 等
"""

import click
import sys
import os


def _lazy_import_cphasing_main():
    """延迟加载cphasing主函数|Lazy load cphasing main function"""
    try:
        from ...cphasing.main import main as cphasing_main
        return cphasing_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


@click.command(
    short_help='CPhasing基因组分相和挂载工具|CPhasing phasing and scaffolding tool',
    context_settings=dict(
        help_option_names=['-h', '--help'],
        max_content_width=120,
        ignore_unknown_options=True,
    )
)
@click.argument('subcommand', nargs=-1, type=click.UNPROCESSED)
@click.option('-i', '--input', 'fasta',
              help='基因组FASTA文件|Genome FASTA file')
@click.option('--hic1',
              help='Hi-C R1 reads文件|Hi-C R1 reads file')
@click.option('--hic2',
              help='Hi-C R2 reads文件|Hi-C R2 reads file')
@click.option('-n', '--groups',
              default='0',
              help='分组数|Number of groups (例如: "8:4" | e.g., "8:4", "0"=自动|auto)')
@click.option('-t', '--threads',
              type=int,
              default=12,
              show_default=True,
              help='线程数|Number of threads')
@click.option('--mode',
              type=click.Choice(['phasing', 'haploid', 'basal', 'basal_withprune']),
              default='phasing',
              show_default=True,
              help='分相模式|Phasing mode')
@click.option('--preset',
              type=click.Choice(['precision', 'sensitive', 'very-sensitive', 'nofilter']),
              default='precision',
              show_default=True,
              help='分析预设|Analysis preset')
@click.option('-o', '--output-dir',
              default='./cphasing_output',
              show_default=True,
              help='输出目录|Output directory')
@click.option('--steps',
              help='运行指定步骤|Run specified steps (例如: "1,2,3" | e.g., "1,2,3")')
@click.option('--skip-steps',
              help='跳过步骤|Skip steps (例如: "1,2" | e.g., "1,2")')
@click.option('--hic-aligner',
              type=click.Choice(['_chromap', 'chromap', 'minimap2', 'bwa-mem2']),
              default='_chromap',
              show_default=True,
              help='Hi-C比对器|Hi-C aligner')
@click.option('--hic-mapper-k',
              type=int,
              help='Hi-C mapper kmer大小|Hi-C mapper kmer size')
@click.option('--hic-mapper-w',
              type=int,
              help='Hi-C mapper窗口大小|Hi-C mapper window size')
@click.option('--mapping-quality',
              type=int,
              default=0,
              show_default=True,
              help='最小比对质量|Minimum mapping quality')
@click.option('--hcr',
              is_flag=True,
              help='启用高置信区域|Enable high confidence regions')
@click.option('--pattern',
              help='酶切位点模式|Restriction enzyme pattern (例如: AAGCTT)')
@click.option('--low-memory',
              is_flag=True,
              help='低内存模式|Low memory mode')
def cphasing(subcommand, fasta, hic1, hic2, groups, threads, mode, preset,
             output_dir, steps, skip_steps, hic_aligner, hic_mapper_k,
             hic_mapper_w, mapping_quality, hcr, pattern, low_memory):
    """
    CPhasing基因组分相和挂载工具|CPhasing Genome Phasing and Scaffolding Tool

    支持所有CPhasing子命令|Supports all CPhasing subcommands

    示例|Examples: biopytools cphasing -i genome.fa --hic1 R1.fq.gz --hic2 R2.fq.gz -t 12 -n 16:2
    """
    cphasing_main = _lazy_import_cphasing_main()

    # 识别子命令：subcommand元组中的非选项参数作为CPhasing子命令
    # Identify subcommand: non-option args in subcommand tuple as CPhasing subcommand
    cphasing_subcommand = 'pipeline'
    subcommand_parts = []
    for arg in subcommand:
        if arg.startswith('-'):
            subcommand_parts.append(arg)
        elif not cphasing_subcommand or cphasing_subcommand == 'pipeline':
            cphasing_subcommand = arg
        else:
            subcommand_parts.append(arg)

    # 构建参数列表|Build argument list
    args = ['cphasing.py', cphasing_subcommand]

    # pipeline 专用参数|pipeline-specific parameters
    if cphasing_subcommand == 'pipeline':
        if fasta:
            args.extend(['-i', fasta])
        if hic1:
            args.extend(['--hic1', hic1])
        if hic2:
            args.extend(['--hic2', hic2])

    # 通用参数|Generic parameters
    if groups != '0':
        args.extend(['-n', groups])
    if threads != 12:
        args.extend(['-t', str(threads)])
    if mode != 'phasing':
        args.extend(['--mode', mode])
    if preset != 'precision':
        args.extend(['--preset', preset])
    if output_dir != './cphasing_output':
        args.extend(['-o', output_dir])

    if steps:
        args.extend(['--steps', steps])
    if skip_steps:
        args.extend(['--skip-steps', skip_steps])
    if hic_aligner != '_chromap':
        args.extend(['--hic-aligner', hic_aligner])
    if hic_mapper_k is not None:
        args.extend(['--hic-mapper-k', str(hic_mapper_k)])
    if hic_mapper_w is not None:
        args.extend(['--hic-mapper-w', str(hic_mapper_w)])
    if mapping_quality != 0:
        args.extend(['--mapping-quality', str(mapping_quality)])
    if hcr:
        args.append('--hcr')
    if pattern:
        args.extend(['--pattern', pattern])
    if low_memory:
        args.append('--low-memory')

    # 透传子命令参数|Pass-through subcommand arguments
    args.extend(subcommand_parts)

    original_argv = sys.argv
    sys.argv = args

    try:
        cphasing_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
