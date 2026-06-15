"""
TeloComp端粒鉴定命令|TeloComp Telomere Identification Command
"""

import click
import sys
import os


def _lazy_import_telocomp_main():
    """延迟加载telocomp主函数|Lazy load telocomp main function"""
    try:
        from ...telocomp.main import main as telocomp_main
        return telocomp_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_file_exists(file_path):
    """验证文件存在(仅在非帮助模式)|Validate file existence (only in non-help mode)"""
    if not _is_help_request() and file_path and not os.path.exists(file_path):
        raise click.BadParameter(f"文件不存在|File does not exist: {file_path}")
    return file_path


@click.command(
    short_help='TeloComp端粒鉴定和可视化工具|TeloComp telomere identification and visualization tool',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('-g', '--genome',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='基因组FASTA文件|Genome FASTA file')
@click.option('-o', '--output-dir',
              required=True,
              type=click.Path(),
              help='输出目录|Output directory')
@click.option('--ont',
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='ONT数据文件|ONT data file')
@click.option('--hifi',
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='HiFi数据文件|HiFi data file')
@click.option('-t', '--threads',
              type=int,
              default=12,
              show_default=True,
              help='线程数|Number of threads')
@click.option('-c', '--coverage',
              type=int,
              default=100,
              show_default=True,
              help='覆盖度参数(0-100)|Coverage parameter (0-100)')
@click.option('-m', '--motif',
              default='CCCTAAA',
              show_default=True,
              help='端粒重复序列(植物默认CCCTAAA, 动物TTAGGG)|Telomeric repeat sequence (plant: CCCTAAA, animal: TTAGGG)')
@click.option('-M', '--motif-num',
              type=int,
              default=7,
              show_default=True,
              help='端粒重复序列碱基数|Number of bases in telomere motif')
@click.option('--skip-filter',
              is_flag=True,
              help='跳过Filter步骤|Skip filter steps')
@click.option('--no-visualization',
              is_flag=True,
              help='跳过可视化|Skip visualization')
def telocomp(genome, output_dir, ont, hifi, threads, coverage, motif,
             motif_num, skip_filter, no_visualization):
    """
    TeloComp端粒鉴定和可视化工具|TeloComp Telomere Identification and Visualization Tool

    基于ONT/HiFi数据进行端粒序列鉴定、提取和可视化分析|Identify, extract and visualize telomere sequences using ONT/HiFi data

    示例|Example: biopytools telocomp -g genome.fa -o output --ont ont.fastq.gz --hifi hifi.fastq.gz
    """

    # 延迟加载|Lazy loading
    telocomp_main = _lazy_import_telocomp_main()

    # 构建参数列表|Build argument list
    args = ['telocomp.py']

    # 必需参数|Required parameters
    args.extend(['-g', genome])
    args.extend(['-o', output_dir])

    # 可选参数|Optional parameters
    if ont:
        args.extend(['--ont', ont])

    if hifi:
        args.extend(['--hifi', hifi])

    if threads != 12:
        args.extend(['-t', str(threads)])

    if coverage != 100:
        args.extend(['-c', str(coverage)])

    if motif != 'CCCTAAA':
        args.extend(['-m', motif])

    if motif_num != 7:
        args.extend(['-M', str(motif_num)])

    # 布尔选项|Boolean options
    if skip_filter:
        args.append('--skip-filter')

    if no_visualization:
        args.append('--no-visualization')

    # 执行主程序|Execute main program
    original_argv = sys.argv
    sys.argv = args

    try:
        telocomp_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
