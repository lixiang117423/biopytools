"""
ENA数据下载命令|ENA Data Download Command
"""

import click
import sys
import os


def _lazy_import_ena_main():
    """延迟加载ena主函数|Lazy load ena main function"""
    try:
        from ...ena_downloader.main import main as ena_main
        return ena_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_aspera_key(file_path):
    """验证Aspera密钥文件|Validate aspera key file existence"""
    if not _is_help_request() and file_path and not os.path.exists(file_path):
        raise click.BadParameter(f"Aspera密钥文件不存在|Aspera key file does not exist: {file_path}")
    return file_path


@click.command(short_help='ENA数据下载|ENA Data Download',
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--accession', '-a',
              required=True,
              help='ENA项目编号|ENA accession number')
@click.option('--output-dir', '-o',
              type=click.Path(),
              help='输出目录|Output directory')
@click.option('--create-dir', '-d',
              is_flag=True,
              help='创建专门输出目录|Create dedicated output directory')
@click.option('--metadata-format', '-f',
              type=click.Choice(['tsv', 'csv', 'xlsx']),
              default='tsv',
              show_default=True,
              help='元数据文件格式|Metadata file format')
@click.option('--protocol', '-p',
              type=click.Choice(['ftp', 'aspera']),
              default='ftp',
              show_default=True,
              help='下载协议类型|Download protocol type')
@click.option('--aspera-key', '-k',
              callback=lambda ctx, param, value: _validate_aspera_key(value) if value else None,
              help='Aspera私钥路径|Path to aspera private key')
@click.option('--method', '-m',
              type=click.Choice(['save', 'run']),
              default='save',
              show_default=True,
              help='执行模式|Execution mode')
@click.option('--metadata-only', '-M',
              is_flag=True,
              help='仅下载元数据|Only download metadata')
@click.option('--fields', '-F',
              multiple=True,
              help='自定义元数据字段|Custom metadata fields')
@click.option('--max-retries', '-r',
              type=int,
              default=3,
              show_default=True,
              help='API请求最大重试次数|Maximum API request retries')
def ena_downloader(accession, output_dir, create_dir, metadata_format, protocol,
                  aspera_key, method, metadata_only, fields, max_retries):
    """
    ENA数据下载工具|ENA Data Download Tool

    从ENA数据库下载FASTQ数据及元信息|Download FASTQ data and metadata from ENA database

    示例|Examples: biopytools ena-download -a PRJNA661210
    """

    # 延迟加载|Lazy loading
    ena_main = _lazy_import_ena_main()

    # 构建参数列表|Build argument list
    args = ['ena_downloader.py']

    # 必需参数|Required parameters
    args.extend(['--accession', accession])

    # 可选参数|Optional parameters
    if output_dir:
        args.extend(['--output-dir', output_dir])

    if create_dir:
        args.append('--create-dir')

    if metadata_format != 'tsv':
        args.extend(['--metadata-format', metadata_format])

    if protocol != 'ftp':
        args.extend(['--protocol', protocol])

    if aspera_key:
        args.extend(['--aspera-key', aspera_key])

    if method != 'save':
        args.extend(['--method', method])

    if metadata_only:
        args.append('--metadata-only')

    if fields:
        args.extend(['--fields'] + list(fields))

    if max_retries != 3:
        args.extend(['--max-retries', str(max_retries)])

    # 保存和恢复sys.argv|Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        # 调用原始主函数|Call original main function
        ena_main()
    except SystemExit as e:
        # 处理正常程序退出|Handle normal program exit
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
