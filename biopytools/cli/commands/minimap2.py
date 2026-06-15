"""
Minimap2序列比对|Minimap2 Sequence Alignment Command
"""

import click
import sys


def _lazy_import_minimap2_main():
    """延迟加载Minimap2主函数|Lazy load Minimap2 main function"""
    try:
        from ...minimap2.main import main as minimap2_main
        return minimap2_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


@click.command(
    short_help='Minimap2序列比对和未比对区间提取|Minimap2 alignment and unmapped region extraction',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--target', '-t',
              required=True,
              type=click.Path(exists=True),
              help='目标基因组文件路径|Target genome file path')
@click.option('--query', '-q',
              required=True,
              type=click.Path(exists=True),
              help='查询基因组文件路径|Query genome file path')
@click.option('--output-dir', '-o',
              default='./minimap2_output',
              show_default=True,
              type=click.Path(),
              help='输出目录|Output directory')
@click.option('--preset', '-x',
              default='asm5',
              show_default=True,
              type=click.Choice(['asm5', 'asm10', 'asm20', 'map-ont', 'map-pb']),
              help='Minimap2预设参数|Minimap2 preset parameters')
@click.option('--threads', '-p',
              default=12,
              show_default=True,
              type=int,
              help='线程数|Number of threads')
@click.option('--min-match', '-m',
              default=1000,
              show_default=True,
              type=int,
              help='最小匹配长度阈值|Minimum match length threshold')
@click.option('--min-unmapped', '-u',
              default=1000,
              show_default=True,
              type=int,
              help='最小未比对区间长度阈值|Minimum unmapped region length threshold')
@click.option('--tp-type',
              default='P',
              show_default=True,
              type=click.Choice(['S', 'P', 'SP']),
              help='保留的比对类型(tp类型)|Alignment type to keep: S(secondary), P(primary), SP(both)')
@click.option('--minimap2-path', '-M',
              default='minimap2',
              show_default=True,
              type=str,
              help='minimap2可执行文件路径|minimap2 executable path')
@click.option('--seqkit-path', '-S',
              default='seqkit',
              show_default=True,
              type=str,
              help='seqkit可执行文件路径|seqkit executable path')
def minimap2(target, query, output_dir, preset, threads, min_match, min_unmapped,
             tp_type, minimap2_path, seqkit_path):
    """
    Minimap2序列比对和未比对区间提取工具|Minimap2 Alignment and Unmapped Region Extraction Tool

    使用Minimap2进行序列比对,提取未比对区间和序列
    Perform sequence alignment with Minimap2, extract unmapped regions and sequences

    示例|Examples: biopytools minimap2 -t target_genome.fasta -q query_genome.fasta -o results/
    """

    # 延迟加载|Lazy loading
    minimap2_main = _lazy_import_minimap2_main()

    # 构建主函数参数列表|Build argument list for main function
    args = ['minimap2.py']

    # 必需参数|Required parameters
    args.extend(['-t', target])
    args.extend(['-q', query])

    # 可选参数(仅在非默认时添加)|Optional parameters (add only when non-default)
    if output_dir != './minimap2_output':
        args.extend(['-o', output_dir])

    if preset != 'asm5':
        args.extend(['-x', preset])

    if threads != 8:
        args.extend(['-p', str(threads)])

    if min_match != 1000:
        args.extend(['-m', str(min_match)])

    if min_unmapped != 1000:
        args.extend(['-u', str(min_unmapped)])

    if tp_type != 'P':
        args.extend(['--tp-type', tp_type])

    if minimap2_path != 'minimap2':
        args.extend(['-M', minimap2_path])

    if seqkit_path != 'seqkit':
        args.extend(['-S', seqkit_path])

    # 保存并恢复sys.argv|Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        # 调用主函数|Call main function
        minimap2_main()
    except SystemExit as e:
        # 处理系统退出|Handle system exit
        if e.code != 0:
            sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n用户中断|Analysis pipeline interrupted by user", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"运行时错误|Runtime error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
