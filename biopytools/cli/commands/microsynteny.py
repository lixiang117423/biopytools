"""
微观共线性分析CLI包装器|Microsynteny Analysis CLI Wrapper
"""

import click
import sys
import os


def _lazy_import_microsynteny_main():
    """延迟加载microsynteny主函数|Lazy load microsynteny main function"""
    try:
        from ...microsynteny.main import main as microsynteny_main
        return microsynteny_main
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


def _validate_dir_exists(dir_path):
    """验证目录存在(仅在非帮助模式)|Validate directory exists (only in non-help mode)"""
    if not _is_help_request() and dir_path and not os.path.exists(dir_path):
        raise click.BadParameter(f"目录不存在|Directory does not exist: {dir_path}")
    return dir_path


@click.command(
    short_help='微观共线性分析工具|Microsynteny analysis tool',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('-i', '--genome-folder',
              required=True,
              callback=lambda ctx, param, value: _validate_dir_exists(value) if value else None,
              help='基因组文件夹路径|Genome folder path (包含A.fa, A.gff等|containing A.fa, A.gff, etc.)')
@click.option('-g', '--gene-list',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='目标基因列表文件|Target gene list file (两列|two cols: species_id gene_id)')
@click.option('-o', '--output-dir',
              default='./microsynteny_output',
              show_default=True,
              help='输出目录路径|Output directory path')
@click.option('-j', '--jcvi-path',
              default='~/miniforge3/envs/jcvi_v.1.5.7',
              show_default=True,
              help='JCVI环境路径|JCVI environment path')
@click.option('-t', '--threads',
              type=int,
              default=12,
              show_default=True,
              help='线程数|Number of threads')
@click.option('--extend-genes',
              type=int,
              default=30,
              show_default=True,
              help='延伸基因数|Number of genes to extend on each side')
@click.option('--cscore',
              type=float,
              default=0.99,
              show_default=True,
              help='共线性分数阈值|Synteny score threshold (0-1)')
@click.option('--step',
              type=click.Choice(['1', '2', '3', '4']),
              help='运行指定步骤|Run specific step only:\n'
                   '1: 数据预处理|Data preprocessing\n'
                   '2: 共线性分析|Synteny analysis\n'
                   '3: 区块提取|Block extraction\n'
                   '4: 绘图|Plotting')
@click.option('--log-level',
              type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL']),
              default='INFO',
              show_default=True,
              help='日志级别|Log level')
def microsynteny(genome_folder, gene_list, output_dir, jcvi_path,
                 threads, extend_genes, cscore, step, log_level):
    """
    微观共线性分析工具|Microsynteny Analysis Tool

    基于JCVI的自动化微观共线性分析和可视化
    Automated microsynteny analysis and visualization based on JCVI

    示例|Examples: biopytools microsynteny -i ./genome_data -g genes.txt
    """
    # 延迟加载|Lazy loading
    microsynteny_main = _lazy_import_microsynteny_main()

    # 构建参数列表|Build argument list
    args = ['microsynteny.py']

    # 必需参数|Required parameters
    args.extend(['-i', genome_folder])
    args.extend(['-g', gene_list])

    # 可选参数|Optional parameters
    if output_dir != './microsynteny_output':
        args.extend(['-o', output_dir])

    if jcvi_path != '~/miniforge3/envs/jcvi_v.1.5.7':
        args.extend(['-j', jcvi_path])

    if threads != 12:
        args.extend(['-t', str(threads)])

    if extend_genes != 30:
        args.extend(['--extend-genes', str(extend_genes)])

    if cscore != 0.99:
        args.extend(['--cscore', str(cscore)])

    if step:
        args.extend(['--step', step])

    if log_level != 'INFO':
        args.extend(['--log-level', log_level])

    # 执行主程序|Execute main program
    original_argv = sys.argv
    sys.argv = args

    try:
        microsynteny_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
