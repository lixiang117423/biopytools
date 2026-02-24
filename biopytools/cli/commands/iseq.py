"""
iSeq下载命令|iSeq Download Command

"""

import click
import sys
import os


def _lazy_import_iseq_main():
    """延迟加载iseq主函数|Lazy load iseq main function"""
    try:
        from ...iseq.main import main as iseq_main
        return iseq_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


@click.command(
    short_help='公共测序数据下载工具|Public sequencing data download tool',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('-i', '--accession',
              required=True,
              help='项目/样本/实验ID|Project/Sample/Experiment ID (e.g., PRJNA1014406)')
@click.option('-o', '--output-dir',
              default='./iseq_output',
              show_default=True,
              type=click.Path(),
              help='输出目录|Output directory')
@click.option('--iseq-path',
              default='/share/org/YZWL/yzwl_lixg/miniforge3/envs/iseq_v.1.9.8/bin/iseq',
              show_default=True,
              help='iSeq软件路径|iSeq software path')
@click.option('-c', '--conda-env',
              default='iseq_v.1.9.8',
              show_default=True,
              help='Conda环境名|Conda environment name')
@click.option('-m', '--metadata-only',
              is_flag=True,
              help='仅下载元数据|Only download metadata')
@click.option('-g', '--gzip',
              is_flag=True,
              help='下载gzip格式FASTQ|Download FASTQ in gzip format')
@click.option('-q', '--fastq',
              is_flag=True,
              help='转换为FASTQ格式|Convert to FASTQ format')
@click.option('-e', '--merge',
              type=click.Choice(['ex', 'sa', 'st']),
              help='合并选项|Merge option (ex/sa/st)')
@click.option('-t', '--threads',
              default=16,
              show_default=True,
              type=int,
              help='线程数|Number of threads')
@click.option('-p', '--parallel',
              default=10,
              show_default=True,
              type=int,
              help='并行连接数|Number of parallel connections')
@click.option('-s', '--speed',
              type=int,
              help='下载速度限制(MB/s)|Download speed limit (MB/s)')
@click.option('-d', '--database',
              default='ena',
              show_default=True,
              type=click.Choice(['ena', 'sra']),
              help='数据库选择|Database selection')
@click.option('--protocol',
              default='ftp',
              show_default=True,
              type=click.Choice(['ftp', 'https']),
              help='协议选择|Protocol selection')
@click.option('-a', '--use-aspera',
              is_flag=True,
              help='使用Aspera高速下载|Use Aspera for high-speed download')
@click.option('--skip-md5',
              is_flag=True,
              help='跳过MD5校验|Skip MD5 check')
@click.option('--quiet',
              is_flag=True,
              help='静默模式|Quiet mode')
def iseq(accession, output_dir, iseq_path, conda_env, metadata_only, gzip, fastq, merge,
         threads, parallel, speed, database, protocol, use_aspera, skip_md5, quiet):
    """
    基于iSeq的公共测序数据下载工具|Public Sequencing Data Download Tool Based on iSeq

    支持从GSA/SRA/ENA/DDBJ数据库下载测序数据和元数据，支持Aspera高速下载
    Download sequencing data and metadata from GSA/SRA/ENA/DDBJ databases with Aspera high-speed download support

    示例|Examples: biopytools iseq -i PRJNA1014406 -g -a -p 10 -t 16 -o ./data
    """

    # 延迟加载|Lazy loading
    iseq_main = _lazy_import_iseq_main()

    # 构建参数列表|Build argument list
    args = ['iseq.py']

    # 必需参数|Required parameters
    args.extend(['-i', accession])

    # 可选参数|Optional parameters
    if output_dir != './iseq_output':
        args.extend(['-o', output_dir])

    if iseq_path != '/share/org/YZWL/yzwl_lixg/miniforge3/envs/iseq_v.1.9.8/bin/iseq':
        args.extend(['--iseq-path', iseq_path])

    if conda_env != 'iseq_v.1.9.8':
        args.extend(['-c', conda_env])

    # 下载选项|Download options
    if metadata_only:
        args.append('-m')

    if gzip:
        args.append('-g')

    if fastq:
        args.append('-q')

    if merge:
        args.extend(['-e', merge])

    # 性能参数|Performance parameters
    if threads != 16:
        args.extend(['-t', str(threads)])

    if parallel != 10:
        args.extend(['-p', str(parallel)])

    if speed:
        args.extend(['-s', str(speed)])

    # 数据库选项|Database options
    if database != 'ena':
        args.extend(['-d', database])

    if protocol != 'ftp':
        args.extend(['--protocol', protocol])

    # 高级选项|Advanced options
    if use_aspera:
        args.append('--use-aspera')

    if skip_md5:
        args.append('--skip-md5')

    if quiet:
        args.append('--quiet')

    # 执行主程序|Execute main program
    original_argv = sys.argv
    sys.argv = args

    try:
        iseq_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
