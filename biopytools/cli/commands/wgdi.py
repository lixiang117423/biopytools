"""
WGDI命令行接口|WGDI Command Line Interface
"""

import click
import sys


def _lazy_import_wgdi_main():
    """延迟加载wgdi主函数|Lazy load wgdi main function"""
    try:
        from ...wgdi.main import main as wgdi_main
        return wgdi_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


@click.command(
    short_help='WGDI比较基因组学分析工具|WGDI comparative genomics tool',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.argument('subcommand',
                type=click.Choice(['dotplot', 'collinearity', 'calks']),
                required=True)
@click.option('-b', '--blast',
              help='BLAST结果文件|BLAST result file')
@click.option('--gff1',
              help='物种1 GFF文件|Species 1 GFF file')
@click.option('--gff2',
              help='物种2 GFF文件|Species 2 GFF file')
@click.option('--lens1',
              help='物种1 LENS文件|Species 1 LENS file')
@click.option('--lens2',
              help='物种2 LENS文件|Species 2 LENS file')
@click.option('-c', '--collinearity',
              help='共线性分析结果文件|Collinearity analysis result file')
@click.option('--fasta1',
              help='物种1 FASTA文件|Species 1 FASTA file')
@click.option('--fasta2',
              help='物种2 FASTA文件|Species 2 FASTA file')
@click.option('-o', '--output-dir',
              default='./wgdi_output',
              show_default=True,
              help='输出目录|Output directory')
@click.option('-t', '--threads',
              type=int,
              default=8,
              show_default=True,
              help='线程数|Number of threads')
@click.option('--wgdi-path',
              help='WGDI软件路径|WGDI software path')
def wgdi(subcommand, blast, gff1, gff2, lens1, lens2, collinearity,
         fasta1, fasta2, output_dir, threads, wgdi_path):
    """
    WGDI比较基因组学分析工具|WGDI Comparative Genomics Analysis Tool

    示例|Examples: biopytools wgdi dotplot -b blast.txt --gff1 s1.gff --gff2 s2.gff --lens1 s1.lens --lens2 s2.lens
    """
    # 延迟加载|Lazy loading
    wgdi_main = _lazy_import_wgdi_main()

    # 构建参数列表|Build argument list
    args = ['wgdi.py', subcommand]

    # 添加通用参数|Add common parameters
    if output_dir != './wgdi_output':
        args.extend(['-o', output_dir])

    if threads != 8:
        args.extend(['-t', str(threads)])

    if wgdi_path:
        args.extend(['--wgdi-path', wgdi_path])

    # 添加子命令特定参数|Add subcommand-specific parameters
    if subcommand in ['dotplot', 'collinearity']:
        if blast:
            args.extend(['-b', blast])
        if gff1:
            args.extend(['--gff1', gff1])
        if gff2:
            args.extend(['--gff2', gff2])
        if lens1:
            args.extend(['--lens1', lens1])
        if lens2:
            args.extend(['--lens2', lens2])

    elif subcommand == 'calks':
        if collinearity:
            args.extend(['-c', collinearity])
        if fasta1:
            args.extend(['--fasta1', fasta1])
        if fasta2:
            args.extend(['--fasta2', fasta2])

    # 执行主程序|Execute main program
    original_argv = sys.argv
    sys.argv = args

    try:
        wgdi_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
