"""泛基因组Block构建工具CLI命令|Pan-Blocks Construction Tool CLI Commands"""

import click
import sys
import os


def _lazy_import_main():
    """延迟加载主函数|Lazy load main function"""
    try:
        from ...pan_blocks.main import main as pan_blocks_main
        return pan_blocks_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_path_exists(path):
    """验证路径存在|Validate path exists"""
    if not _is_help_request() and path and not os.path.exists(path):
        raise click.BadParameter(f"路径不存在|Path does not exist: {path}")
    return path


def _build_common_options():
    """构建公共选项装饰器列表|Build common option decorators"""
    return []


@click.group(
    short_help='泛基因组Block构建工具|Pan-genome Block Construction Tool',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
def pan_blocks():
    """泛基因组Block构建工具 - 通过迭代两两比对构建泛基因组共线性块|Pan-genome Block Construction Tool - Build pan-genome syntenic blocks via iterative pairwise alignments.

    示例|Examples: biopytools pan-blocks all -i genome_list.txt -o output_dir/
    """
    pass


@pan_blocks.command(short_help='两两基因组比对|Pairwise genome alignment')
@click.option('-i', '--genome-list', required=True,
              callback=lambda ctx, param, value: _validate_path_exists(value),
              help='基因组列表文件|Genome list file (name<TAB>path)')
@click.option('-o', '--output-dir', default='./pan_blocks_output',
              help='输出目录|Output directory')
@click.option('-t', '--threads', default=12, type=int,
              help='线程数|Number of threads')
@click.option('--parallel-alignments', default=4, type=int,
              help='并行比对数|Parallel alignments')
@click.option('--min-alignment-length', default=10000, type=int,
              help='最小比对长度|Min alignment length')
def align(**kwargs):
    """运行两两基因组比对|Run pairwise genome alignment step.

    示例|Examples: biopytools pan-blocks align -i genome_list.txt -o output_dir/
    """
    pan_blocks_main = _lazy_import_main()

    args = ['pan_blocks.py', '--step', 'align']
    args.extend(['-i', kwargs['genome_list']])
    if kwargs['output_dir'] != './pan_blocks_output':
        args.extend(['-o', kwargs['output_dir']])
    if kwargs['threads'] != 12:
        args.extend(['-t', str(kwargs['threads'])])
    if kwargs['parallel_alignments'] != 4:
        args.extend(['--parallel-alignments', str(kwargs['parallel_alignments'])])
    if kwargs['min_alignment_length'] != 10000:
        args.extend(['--min-alignment-length', str(kwargs['min_alignment_length'])])

    _execute_main(args, pan_blocks_main)


@pan_blocks.command(short_help='构建泛基因组Blocks|Build pan-genome blocks')
@click.option('-i', '--genome-list', required=True,
              callback=lambda ctx, param, value: _validate_path_exists(value),
              help='基因组列表文件|Genome list file (name<TAB>path)')
@click.option('-o', '--output-dir', default='./pan_blocks_output',
              help='输出目录|Output directory')
@click.option('--genome-order',
              help='基因组优先级顺序文件|Genome priority order file')
@click.option('--chromosome',
              help='指定染色体|Specific chromosome')
def build(**kwargs):
    """构建泛基因组Blocks|Build pan-genome blocks from pairwise alignments.

    示例|Examples: biopytools pan-blocks build -i genome_list.txt -o output_dir/
    """
    pan_blocks_main = _lazy_import_main()

    args = ['pan_blocks.py', '--step', 'build']
    args.extend(['-i', kwargs['genome_list']])
    if kwargs['output_dir'] != './pan_blocks_output':
        args.extend(['-o', kwargs['output_dir']])
    if kwargs['genome_order']:
        args.extend(['--genome-order', kwargs['genome_order']])
    if kwargs['chromosome']:
        args.extend(['--chromosome', kwargs['chromosome']])

    _execute_main(args, pan_blocks_main)


@pan_blocks.command(short_help='可视化泛基因组Blocks|Visualize pan-genome blocks')
@click.option('-i', '--genome-list', required=True,
              callback=lambda ctx, param, value: _validate_path_exists(value),
              help='基因组列表文件|Genome list file (name<TAB>path)')
@click.option('-o', '--output-dir', default='./pan_blocks_output',
              help='输出目录|Output directory')
@click.option('--genome-order',
              help='基因组优先级顺序文件|Genome priority order file')
@click.option('--plot-format', default='svg', type=click.Choice(['svg', 'png']),
              help='绘图格式|Plot format')
def plot(**kwargs):
    """可视化泛基因组Blocks|Visualize pan-genome blocks with synteny links.

    示例|Examples: biopytools pan-blocks plot -i genome_list.txt -o output_dir/
    """
    pan_blocks_main = _lazy_import_main()

    args = ['pan_blocks.py', '--step', 'plot']
    args.extend(['-i', kwargs['genome_list']])
    if kwargs['output_dir'] != './pan_blocks_output':
        args.extend(['-o', kwargs['output_dir']])
    if kwargs['genome_order']:
        args.extend(['--genome-order', kwargs['genome_order']])
    if kwargs['plot_format'] != 'svg':
        args.extend(['--plot-format', kwargs['plot_format']])

    _execute_main(args, pan_blocks_main)


@pan_blocks.command(short_help='运行完整流程|Run complete pipeline')
@click.option('-i', '--genome-list', required=True,
              callback=lambda ctx, param, value: _validate_path_exists(value),
              help='基因组列表文件|Genome list file (name<TAB>path)')
@click.option('-o', '--output-dir', default='./pan_blocks_output',
              help='输出目录|Output directory')
@click.option('-t', '--threads', default=12, type=int,
              help='线程数|Number of threads')
@click.option('--parallel-alignments', default=4, type=int,
              help='并行比对数|Parallel alignments')
@click.option('--min-alignment-length', default=10000, type=int,
              help='最小比对长度|Min alignment length')
@click.option('--genome-order',
              help='基因组优先级顺序文件|Genome priority order file')
@click.option('--chromosome',
              help='指定染色体|Specific chromosome')
@click.option('--plot-format', default='svg', type=click.Choice(['svg', 'png']),
              help='绘图格式|Plot format')
def all(**kwargs):
    """运行完整的泛基因组Block构建流程|Run complete pan-block construction pipeline (align + build + plot).

    示例|Examples: biopytools pan-blocks all -i genome_list.txt -o output_dir/
    """
    pan_blocks_main = _lazy_import_main()

    args = ['pan_blocks.py']
    args.extend(['-i', kwargs['genome_list']])
    if kwargs['output_dir'] != './pan_blocks_output':
        args.extend(['-o', kwargs['output_dir']])
    if kwargs['threads'] != 12:
        args.extend(['-t', str(kwargs['threads'])])
    if kwargs['parallel_alignments'] != 4:
        args.extend(['--parallel-alignments', str(kwargs['parallel_alignments'])])
    if kwargs['min_alignment_length'] != 10000:
        args.extend(['--min-alignment-length', str(kwargs['min_alignment_length'])])
    if kwargs['genome_order']:
        args.extend(['--genome-order', kwargs['genome_order']])
    if kwargs['chromosome']:
        args.extend(['--chromosome', kwargs['chromosome']])
    if kwargs['plot_format'] != 'svg':
        args.extend(['--plot-format', kwargs['plot_format']])

    _execute_main(args, pan_blocks_main)


def _execute_main(args, pan_blocks_main):
    """执行主程序|Execute main program"""
    original_argv = sys.argv
    sys.argv = args
    try:
        pan_blocks_main()
    except SystemExit as e:
        sys.exit(e.code)
    finally:
        sys.argv = original_argv


@pan_blocks.command(short_help='生成基因组列表文件|Generate genome list file')
@click.option('-i', '--input-dir', required=True,
              callback=lambda ctx, param, value: _validate_path_exists(value),
              help='FASTA文件目录|FASTA files directory')
@click.option('-o', '--output', default='./genome_list.txt',
              help='输出文件路径|Output file path')
def prepare(**kwargs):
    """从FASTA文件目录自动生成genome_list.txt|Auto-generate genome_list.txt from FASTA directory.

    示例|Examples: biopytools pan-blocks prepare -i /path/to/genomes/ -o genome_list.txt
    """
    pan_blocks_main = _lazy_import_main()

    args = ['pan_blocks.py', '--prepare', kwargs['input_dir']]
    if kwargs['output'] != './genome_list.txt':
        args.extend(['--prepare-output', kwargs['output']])

    _execute_main(args, pan_blocks_main)
