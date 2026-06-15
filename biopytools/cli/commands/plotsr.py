"""
PlotSR命令|PlotSR Command
"""

import click
import sys
import os


def _lazy_import_plotsr_main():
    """延迟加载plotsr主函数|Lazy load plotsr main function"""
    try:
        from ...plotsr.main import main as plotsr_main
        return plotsr_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_path(path):
    """验证路径|Validate path"""
    if not _is_help_request() and path:
        # path可能是元组(multiple=True)或字符串
        # path can be tuple (multiple=True) or string
        if isinstance(path, tuple):
            for p in path:
                if not os.path.exists(p):
                    raise click.BadParameter(f"路径不存在|Path does not exist: {p}")
        else:
            if not os.path.exists(path):
                raise click.BadParameter(f"路径不存在|Path does not exist: {path}")
    return path


@click.command(
    short_help="多基因组共线性可视化|Multi-genome synteny visualization",
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('-i', '--input',
              multiple=True,
              required=True,
              callback=lambda ctx, param, value: _validate_path(value) if value else value,
              help='输入基因组FASTA文件（可多次使用）、包含基因组的文件夹或map文件|'
                   'Input genome FASTA files (multiple), folder containing genomes, or map file')
@click.option('-o', '--output-dir',
              required=True,
              help='输出目录|Output directory')
@click.option('-n', '--names',
              help='基因组名称（逗号分隔）|Genome names (comma-separated)')
@click.option('-t', '--threads',
              type=int,
              default=12,
              show_default=True,
              help='线程数|Number of threads')
@click.option('--minimap2-preset',
              type=click.Choice(['asm5', 'asm10', 'asm20']),
              default='asm5',
              show_default=True,
              help='minimap2预设参数|minimap2 preset')
@click.option('-s', '--min-sr-size',
              type=int,
              default=10000,
              show_default=True,
              help='最小结构变异大小|Minimum structural variant size')
@click.option('--output-format',
              type=click.Choice(['pdf', 'png', 'svg']),
              default='pdf',
              show_default=True,
              help='输出格式|Output format')
@click.option('-f', '--font-size',
              type=int,
              default=6,
              show_default=True,
              help='字体大小|Font size')
@click.option('-d', '--dpi',
              type=int,
              default=300,
              show_default=True,
              help='图片DPI|Image DPI')
@click.option('--space-ratio',
              type=float,
              default=0.7,
              show_default=True,
              help='同源染色体间距(0.1-0.75)|Space for homologous chromosomes')
@click.option('-v', '--vertical',
              is_flag=True,
              default=False,
              help='垂直排列染色体|Plot vertical chromosomes')
@click.option('--itx',
              is_flag=True,
              default=False,
              help='染色体间交互模式|Inter-chromosomal plotting mode')
@click.option('--nosyn',
              is_flag=True,
              default=False,
              help='不绘制同源区域|Do not plot syntenic regions')
@click.option('--noinv',
              is_flag=True,
              default=False,
              help='不绘制倒位|Do not plot inversions')
@click.option('--notr',
              is_flag=True,
              default=False,
              help='不绘制易位|Do not plot translocations')
@click.option('--nodup',
              is_flag=True,
              default=False,
              help='不绘制重复|Do not plot duplications')
@click.option('-c', '--chr',
              multiple=True,
              help='指定要显示的染色体（可多次使用，支持数字如1或名称如Chr1）|'
                   'Specify chromosomes to display (can be used multiple times, supports number like 1 or name like Chr1)')
@click.option('--skip-existing',
              is_flag=True,
              default=True,
              help='跳过已完成的步骤（默认启用）|Skip completed steps (default: enabled)')
@click.option('--force-run',
              is_flag=True,
              default=False,
              help='强制重新运行所有步骤|Force re-run all steps')
def plotsr(input, output_dir, names, threads, minimap2_preset, min_sr_size, output_format,
           font_size, dpi, space_ratio, vertical, itx, nosyn, noinv, notr, nodup, chr,
           skip_existing, force_run):
    """多基因组共线性可视化工具|Multi-Genome Synteny Visualization Tool

    自动运行minimap2 + SyRI + PlotSR完整流程
    Automatically run minimap2 + SyRI + PlotSR pipeline

    示例|Examples: biopytools plotsr -i genome1.fa -i genome2.fa -o output/
    """

    # 延迟加载|Lazy load: import only when actually called
    plotsr_main_func = _lazy_import_plotsr_main()

    # 构建主函数的参数列表|Build argument list for main function
    args = ['plotsr.py']

    # 输入基因组|Input genomes
    for genome in input:
        args.extend(['-i', genome])

    # 输出目录|Output directory
    args.extend(['-o', output_dir])

    # 基因组名称|Genome names
    if names:
        args.extend(['-n', names])

    # 线程和性能参数|Threads and performance parameters
    args.extend(['-t', str(threads)])
    args.extend(['--minimap2-preset', minimap2_preset])

    # 结构变异参数|Structural variant parameters
    args.extend(['-s', str(min_sr_size)])

    # 输出格式|Output format
    args.extend(['--output-format', output_format])

    # 可视化参数|Visualization parameters
    args.extend(['-f', str(font_size)])
    args.extend(['-d', str(dpi)])
    args.extend(['--space-ratio', str(space_ratio)])

    # 标志参数|Flag parameters
    if vertical:
        args.append('-v')
    if itx:
        args.append('--itx')
    if nosyn:
        args.append('--nosyn')
    if noinv:
        args.append('--noinv')
    if notr:
        args.append('--notr')
    if nodup:
        args.append('--nodup')

    # 染色体过滤参数|Chromosome filtering parameters
    if chr:
        for c in chr:
            args.extend(['-c', c])

    # 流程控制参数|Pipeline control parameters
    if not skip_existing or force_run:
        args.append('--force-run')

    # 保存并恢复sys.argv|Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        # 调用原始主函数|Call original main function
        plotsr_main_func()
    except SystemExit as e:
        # 处理正常程序退出|Handle normal program exit
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
