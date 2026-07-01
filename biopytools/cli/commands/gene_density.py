"""基因密度计算命令|Gene density calculation command"""

import click
import sys
import os


def _lazy_import_main():
    """延迟加载主函数|Lazy load main function"""
    try:
        from ...gene_density.main import main as gene_density_main
        return gene_density_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    return any(arg in {'-h', '--help'} for arg in sys.argv)


def _validate_file_exists(file_path):
    """验证文件存在(仅在非帮助模式)|Validate file exists (only in non-help mode)"""
    if not _is_help_request() and file_path and not os.path.exists(file_path):
        raise click.BadParameter(f"文件不存在|File does not exist: {file_path}")
    return file_path


@click.command(short_help='基因密度计算(每窗口基因数+基因/Mb)|Gene density (genes/window + genes/Mb)',
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('-i', '--gff',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value),
              type=click.Path(),
              help='GFF3注释文件|GFF3 annotation file')
@click.option('-o', '--output-dir',
              default='./gene_density_output',
              show_default=True,
              help='输出目录|Output directory')
@click.option('-w', '--window-size',
              type=int,
              default=100000,
              show_default=True,
              help='窗口大小(bp)|Window size (bp)')
@click.option('--feature-type',
              default='gene',
              show_default=True,
              help='统计的GFF feature类型|GFF feature type to count')
@click.option('-g', '--genome',
              default=None,
              help='染色体长度来源(.fai或FASTA)|Chromosome length source (.fai or FASTA)')
@click.option('--prefix',
              default=None,
              help='输出文件前缀(默认GFF文件名)|Output file prefix (default: GFF stem)')
@click.option('--no-plot',
              is_flag=True,
              help='不绘制密度图|Skip density plot')
def gene_density(gff, output_dir, window_size, feature_type, genome, prefix, no_plot):
    """基因密度计算工具|Gene density calculation tool

    按固定大小窗口统计每条染色体各区间的基因数量与基因密度(基因/Mb)
    |Count genes and density (genes/Mb) per fixed-size window along each chromosome

    示例|Examples: biopytools gene-density -i genome.gff3 -g genome.fa.fai -w 100000
    """

    gene_density_main = _lazy_import_main()

    # 构建参数列表|Build argument list
    args = ['gene_density.py']
    args.extend(['-i', gff])
    args.extend(['-o', output_dir])
    args.extend(['-w', str(window_size)])
    args.extend(['--feature-type', feature_type])
    if genome:
        args.extend(['-g', genome])
    if prefix:
        args.extend(['--prefix', prefix])
    if no_plot:
        args.append('--no-plot')

    # 临时修改sys.argv|Temporarily modify sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        gene_density_main()
    except SystemExit as e:
        if e.code != 0:
            sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n用户中断操作|Operation interrupted by user", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
