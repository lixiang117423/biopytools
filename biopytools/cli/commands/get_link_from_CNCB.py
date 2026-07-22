"""
CNCB测序数据下载链接提取命令|CNCB Sequencing Data Download Links Command
"""

import click
import sys
import os


def _lazy_import_cncb_main():
    """延迟加载CNCB主函数|Lazy load CNCB main function"""
    try:
        from ...get_link_from_CNCB.main import CNCLinkExtractor
        return CNCLinkExtractor
    except ImportError as e:
        click.echo(f"Import Error: {e}", err=True)
        sys.exit(1)


def _validate_file_exists(file_path):
    """验证文件存在|Validate file existence"""
    if not file_path or not os.path.exists(file_path):
        raise click.BadParameter(f"File does not exist: {file_path}")
    return file_path


@click.command(
    short_help='CNCB下载链接提取|CNCB Download Links Extraction',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--input', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              type=click.Path(exists=True, file_okay=True, dir_okay=False, resolve_path=True),
              help='输入文件|Input file path (ProjectID and Run ID)')
@click.option('--output', '-o',
              type=click.Path(file_okay=True, dir_okay=False, resolve_path=True),
              help='输出链接文件|Output links file path')
@click.option('--failed',
              type=click.Path(file_okay=True, dir_okay=False, resolve_path=True),
              help='失败记录文件|Failed records file path')
@click.option('--download-script',
              default='download.sh',
              show_default=True,
              type=click.Path(file_okay=True, dir_okay=False, resolve_path=True),
              help='下载脚本名|Download script filename')
@click.option('--ftp-host',
              default='download2.cncb.ac.cn',
              show_default=True,
              help='FTP服务器|FTP server host')
@click.option('--ftp-timeout',
              default=60,
              show_default=True,
              type=int,
              help='FTP连接超时|FTP connection timeout in seconds')
@click.option('--retry-attempts',
              default=3,
              show_default=True,
              type=int,
              help='FTP重试次数|FTP connection retry attempts')
@click.option('--verbose', '-v',
              is_flag=True,
              help='详细输出模式|Verbose output mode')
@click.option('--log-file',
              type=click.Path(file_okay=True, dir_okay=False, resolve_path=True),
              help='日志文件路径|Log file path')
@click.option('--no-download-script',
              is_flag=True,
              help='不生成下载脚本|Don\'t generate download script')
@click.option('--no-executable',
              is_flag=True,
              help='不设置可执行权限|Don\'t make script executable')
def get_link_from_CNCB(input, output, failed, download_script, ftp_host,
                       ftp_timeout, retry_attempts, verbose,
                       log_file, no_download_script, no_executable):
    """
    CNCB测序数据下载链接提取工具|CNCB Sequencing Data Download Links Extraction Tool

    从CNCB FTP服务器批量搜索和获取测序数据下载链接|Batch search and retrieve download links from CNCB FTP server

    示例|Examples: biopytools get-link-from-CNCB -i projects.txt
    """

    try:
        # 延迟加载|Lazy loading
        CNCLinkExtractor = _lazy_import_cncb_main()

        # 创建提取器|Create extractor
        extractor = CNCLinkExtractor(
            input_file=str(input),
            output_file=str(output) if output else None,
            failed_file=str(failed) if failed else None,
            download_script=str(download_script),
            ftp_host=ftp_host,
            ftp_timeout=ftp_timeout,
            retry_attempts=retry_attempts,
            verbose=verbose,
            log_file=str(log_file) if log_file else None,
            generate_download_script=not no_download_script,
            script_executable=not no_executable
        )

        # 运行提取|Run extraction
        success = extractor.run_extraction()

        if success:
            click.echo("CNCB链接提取完成!|CNCB link extraction completed successfully!")
        else:
            click.echo("CNCB链接提取失败!|CNCB link extraction failed!", err=True)
            sys.exit(1)

    except KeyboardInterrupt:
        click.echo("\n用户中断|Interrupted by user", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"发生错误|Error occurred: {e}", err=True)
        sys.exit(1)
