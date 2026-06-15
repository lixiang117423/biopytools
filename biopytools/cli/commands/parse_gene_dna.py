"""
基因序列提取命令|Gene Sequence Extraction Command
"""

import click
import sys
import os


def _lazy_import_gene_main():
    """延迟加载基因提取主函数|Lazy load gene extraction main function"""
    try:
        from ...parse_gene_dna.main import main as gene_main
        return gene_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_file_exists(file_path):
    """验证文件存在(仅在非帮助模式)|Validate file exists (only in non-help mode)"""
    if not _is_help_request() and file_path and not os.path.exists(file_path):
        raise click.BadParameter(f"文件不存在|File does not exist: {file_path}")
    return file_path


@click.command(
    short_help='基因序列提取工具|Gene Sequence Extraction Tool',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--genome', '-g',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='基因组FASTA文件路径|Genome FASTA file path')
@click.option('--gff', '-f',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='GFF注释文件路径|GFF annotation file path')
@click.option('--output', '-o',
              required=True,
              type=click.Path(),
              help='输出FASTA文件路径|Output FASTA file path')
@click.option('--feature-type',
              default='gene',
              show_default=True,
              help='要提取的特征类型|Feature type to extract')
@click.option('--min-length',
              type=int,
              default=0,
              show_default=True,
              help='最小基因长度过滤|Minimum gene length filter')
@click.option('--threads', '-t',
              type=int,
              default=12,
              show_default=True,
              help='线程数|Number of threads')
@click.option('--line-width',
              type=int,
              default=60,
              show_default=True,
              help='FASTA序列行宽度|FASTA sequence line width')
@click.option('--verbose', '-v',
              is_flag=True,
              help='显示详细信息|Show verbose output')
def parse_gene_dna(genome, gff, output, feature_type, min_length, threads, line_width, verbose):
    """
    基因序列提取工具|Gene Sequence Extraction Tool

    从基因组FASTA文件中根据GFF注释提取基因序列|Extract gene sequences from genome FASTA based on GFF annotation

    示例|Examples: biopytools parse-gene-dna -g genome.fasta -f annotation.gff -o genes.fasta
    """

    # 延迟加载|Lazy load: import only when actually called
    gene_main = _lazy_import_gene_main()

    # 构建主函数的参数列表|Build argument list for main function
    args = ['parse_gene_dna.py']

    # 必需参数|Required parameters
    args.extend(['-g', genome])
    args.extend(['-f', gff])
    args.extend(['-o', output])

    # 可选参数|Optional parameters (add only when non-default)
    if feature_type != 'gene':
        args.extend(['--feature-type', feature_type])
    if min_length != 0:
        args.extend(['--min-length', str(min_length)])
    if threads != 88:
        args.extend(['-t', str(threads)])
    if line_width != 60:
        args.extend(['--line-width', str(line_width)])
    if verbose:
        args.append('-v')

    # 保存并恢复sys.argv|Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        # 调用原始主函数|Call original main function
        gene_main()
    except SystemExit as e:
        # 处理正常程序退出|Handle normal program exit
        if e.code != 0:
            sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n用户中断|User interrupted", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
