"""
ENA数据下载命令 | ENA Data Download Command
优化版本：使用懒加载解决响应速度问题
"""

import click
import sys
import os


def _lazy_import_ena_main():
    """懒加载ena main函数 | Lazy load ena main function"""
    try:
        from ...ena_downloader.main import main as ena_main
        return ena_main
    except ImportError as e:
        click.echo(f"❌ 导入错误 | Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否是帮助请求 | Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_aspera_key(file_path):
    """验证Aspera密钥文件是否存在（仅在非帮助模式下）| Validate aspera key file existence (only in non-help mode)"""
    if not _is_help_request() and file_path and not os.path.exists(file_path):
        raise click.BadParameter(f"Aspera密钥文件不存在 | Aspera key file does not exist: {file_path}")
    return file_path


@click.command(short_help='从ENA下载测序数据样品信息和下载链接',
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--accession', '-a',
              required=True,
              help='ENA项目编号 | ENA accession number (e.g., PRJNA661210, SRP000123)')
@click.option('--output-dir', '-o',
              type=click.Path(),
              help='输出目录 | Output directory')
@click.option('--create-dir', '-d',
              is_flag=True,
              help='创建专门的输出目录 | Create dedicated output directory')
@click.option('--metadata-format', '-f',
              type=click.Choice(['tsv', 'csv', 'xlsx']),
              default='tsv',
              help='元数据文件格式 (默认: tsv) | Metadata file format (default: tsv)')
@click.option('--protocol', '-p',
              type=click.Choice(['ftp', 'aspera']),
              default='ftp',
              help='下载协议类型 (默认: ftp) | Download protocol type (default: ftp)')
@click.option('--aspera-key', '-k',
              callback=lambda ctx, param, value: _validate_aspera_key(value) if value else None,
              help='Aspera私钥路径 | Path to aspera private key')
@click.option('--method', '-m',
              type=click.Choice(['save', 'run']),
              default='save',
              help='执行模式 (默认: save) | Execution mode (default: save)')
@click.option('--metadata-only', '-M',
              is_flag=True,
              help='仅下载元数据，不处理FASTQ文件 | Only download metadata, do not process FASTQ files')
@click.option('--fields', '-F',
              multiple=True,
              help='自定义元数据字段 | Custom metadata fields')
@click.option('--max-retries', '-r',
              type=int,
              default=3,
              help='API请求最大重试次数 (默认: 3) | Maximum API request retries (default: 3)')
def ena_downloader(accession, output_dir, create_dir, metadata_format, protocol, 
                  aspera_key, method, metadata_only, fields, max_retries):
    """
    ENA数据下载工具
    
    从ENA数据库下载测序数据的元数据和FASTQ文件，支持FTP和Aspera协议，
    可生成下载脚本或直接执行下载。
    
    示例 | Examples:
    
    \b
    # 仅下载元数据
    biopytools ena-download -a PRJNA661210 -M
    
    \b
    # 下载元数据并生成FTP下载脚本
    biopytools ena-download -a PRJNA661210 -p ftp -m save
    
    \b
    # 使用Aspera协议直接下载
    biopytools ena-download -a PRJNA661210 -p aspera \\
        -k ~/.aspera/connect/etc/asperaweb_id_dsa.openssh -m run
    
    \b
    # 自定义输出和格式
    biopytools ena-download -a PRJNA661210 -o my_results -f xlsx -d
    
    \b
    # 自定义字段
    biopytools ena-download -a PRJNA661210 -F fastq_ftp -F study_title -f csv
    """
    
    # 🚀 懒加载：只有在实际调用时才导入模块 | Lazy loading: import only when actually called
    ena_main = _lazy_import_ena_main()
    
    # 构建参数列表传递给原始main函数
    args = ['ena_downloader.py']
    
    # 必需参数
    args.extend(['--accession', accession])
    
    # 可选参数（只在非默认值时添加）
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
    
    # 保存并恢复sys.argv
    original_argv = sys.argv
    sys.argv = args
    
    try:
        # 调用原始的main函数
        ena_main()
    except SystemExit as e:
        # 处理程序正常退出
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv