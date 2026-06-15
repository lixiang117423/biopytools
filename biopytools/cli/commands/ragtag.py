"""
RagTag基因组scaffolding命令|RagTag Genome Scaffolding Command
"""

import click
import sys
import os


def _lazy_import_ragtag_main():
    """延迟加载ragtag主函数|Lazy load ragtag main function"""
    try:
        from ...ragtag.main import main as ragtag_main
        return ragtag_main
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
    short_help='RagTag基因组scaffolding工具|RagTag genome scaffolding tool',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('-r', '--reference',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='参考基因组FASTA文件|Reference genome FASTA file')
@click.option('-q', '--query',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='查询基因组FASTA文件|Query genome FASTA file')
@click.option('-s', '--sample-name',
              required=True,
              help='样品名称，用于输出文件命名|Sample name for output file naming')
@click.option('-p', '--prefix',
              default=None,
              show_default=True,
              help='序列ID前缀，会添加到所有输出序列的ID前面|Sequence ID prefix to add to all output sequence IDs')
@click.option('-t', '--threads',
              type=int,
              default=12,
              show_default=True,
              help='线程数|Number of threads')
@click.option('-o', '--output-dir',
              default='./ragtag_output',
              show_default=True,
              type=click.Path(),
              help='输出目录|Output directory')
@click.option('--aligner',
              type=click.Choice(['minimap2', 'unimap', 'nucmer']),
              default='minimap2',
              show_default=True,
              help='比对器|Aligner to use')
@click.option('-C', '--concatenate-unplaced',
              is_flag=True,
              help='将未定位的contigs合并为chr0|Concatenate unplaced contigs into chr0')
@click.option('-R', '--infer-gaps',
              is_flag=True,
              help='推断gap大小|Infer gap sizes')
def ragtag(reference, query, sample_name, prefix, threads, output_dir, aligner,
           concatenate_unplaced, infer_gaps):
    """
    RagTag基因组scaffolding工具|RagTag Genome Scaffolding Tool

    基于参考基因组进行序列scaffolding，自动重命名序列ID并分类输出|Perform reference-based genome scaffolding with automatic sequence ID renaming and categorized output

    示例|Example: biopytools ragtag -r ref.fa -q query.fa -s Sample1 -o output_dir/
    """

    # 延迟加载|Lazy loading
    ragtag_main = _lazy_import_ragtag_main()

    # 构建参数列表|Build argument list
    args = ['ragtag.py']

    # 必需参数|Required parameters
    args.extend(['-r', reference])
    args.extend(['-q', query])
    args.extend(['-s', sample_name])

    # 可选参数|Optional parameters
    if prefix is not None:
        args.extend(['-p', prefix])

    if threads != 12:
        args.extend(['-t', str(threads)])

    if output_dir != './ragtag_output':
        args.extend(['-o', output_dir])

    if aligner != 'minimap2':
        args.extend(['--aligner', aligner])

    # 布尔选项|Boolean options
    if concatenate_unplaced:
        args.append('-C')

    if infer_gaps:
        args.append('-R')

    # 执行主程序|Execute main program
    original_argv = sys.argv
    sys.argv = args

    try:
        ragtag_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
