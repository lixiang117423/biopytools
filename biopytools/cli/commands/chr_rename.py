"""
染色体重命名命令|Chromosome Rename Command

"""

import click
import sys
import os


def _lazy_import_chr_rename_main():
    """延迟加载chr_rename主函数|Lazy load chr_rename main function"""
    try:
        from ...chr_rename.main import main as chr_rename_main
        return chr_rename_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_file_exists(file_path):
    """验证文件存在(仅在非帮助模式)|Validate file existence (only in non-help mode)"""
    if not _is_help_request() and not os.path.exists(file_path):
        raise click.BadParameter(f"文件不存在|File does not exist: {file_path}")
    return file_path


@click.command(
    short_help='基于minimap2的染色体重命名工具|Chromosome rename tool based on minimap2',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--ref', '-r',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='参考基因组FASTA文件|Reference genome FASTA file')
@click.option('--query', '-q',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='待重命名的基因组FASTA文件|Query genome FASTA file to rename')
@click.option('--output-dir', '-o',
              default='./chr_rename_output',
              show_default=True,
              type=click.Path(),
              help='输出目录|Output directory')
@click.option('--minimap2-path', '-a',
              default='minimap2',
              show_default=True,
              help='minimap2软件路径|minimap2 software path')
@click.option('--preset', '-x',
              default='asm5',
              show_default=True,
              type=click.Choice(['asm5', 'asm10', 'asm20']),
              help='minimap2预设模式|minimap2 preset mode')
@click.option('--threads', '-t',
              default=12,
              show_default=True,
              type=int,
              help='线程数|Number of threads')
@click.option('--min-identity', '-i',
              default=0.9,
              show_default=True,
              type=float,
              help='最小序列一致性阈值(0-1)|Minimum identity threshold (0-1)')
@click.option('--min-alignment-length', '-l',
              default=100000,
              show_default=True,
              type=int,
              help='最小比对长度(bp)|Minimum alignment length (bp)')
def chr_rename(ref, query, output_dir, minimap2_path, preset, threads, min_identity, min_alignment_length):
    """
    基于minimap2的染色体重命名工具|Chromosome Rename Tool Based on Minimap2 Alignment

    通过minimap2全基因组比对，将group命名方式的染色体转换为标准Chr命名（自动迭代映射策略）|Convert chromosome names from group format to standard Chr format using minimap2 whole genome alignment (auto iterative mapping strategy)

    示例|Examples: biopytools chr_rename -r ref.fa -q query.fa -o output
    """

    # 延迟加载|Lazy loading
    chr_rename_main = _lazy_import_chr_rename_main()

    # 构建参数列表|Build argument list
    args = ['chr_rename.py']

    # 必需参数|Required parameters
    args.extend(['-r', ref])
    args.extend(['-q', query])

    # 可选参数|Optional parameters
    if output_dir != './chr_rename_output':
        args.extend(['-o', output_dir])

    if minimap2_path != 'minimap2':
        args.extend(['-a', minimap2_path])

    if preset != 'asm5':
        args.extend(['-x', preset])

    if threads != 12:
        args.extend(['-t', str(threads)])

    if min_identity != 0.9:
        args.extend(['-i', str(min_identity)])

    if min_alignment_length != 100000:
        args.extend(['-l', str(min_alignment_length)])

    # 执行主程序|Execute main program
    original_argv = sys.argv
    sys.argv = args

    try:
        chr_rename_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
