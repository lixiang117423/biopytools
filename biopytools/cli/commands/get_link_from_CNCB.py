"""
从CNCB获取测序数据下载链接命令 | CNCB Sequencing Data Download Links Command
"""

import click
import sys
import os

def _lazy_import_cncb_main():
    """懒加载CNCB main函数 | Lazy load CNCB main function"""
    try:
        from ...get_link_from_CNCB.main import CNCLinkExtractor
        return CNCLinkExtractor
    except ImportError as e:
        click.echo(f"❌ 导入错误 | Import Error: {e}", err=True)
        sys.exit(1)

def _validate_file_exists(file_path):
    """验证文件是否存在 | Validate file existence"""
    if not file_path or not os.path.exists(file_path):
        raise click.BadParameter(f"文件不存在 | File does not exist: {file_path}")
    return file_path

@click.command(
    short_help='从CNCB批量获取测序数据下载链接',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--input', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              type=click.Path(exists=True, file_okay=True, dir_okay=False, resolve_path=True),
              help='📁 输入文件路径 (ProjectID和RunID，Tab分隔) | Input file path (ProjectID and Run ID, tab separated)')
@click.option('--output', '-o',
              type=click.Path(file_okay=True, dir_okay=False, resolve_path=True),
              help='📄 输出链接文件路径 | Output links file path (default: [input]_links.txt)')
@click.option('--failed',
              type=click.Path(file_okay=True, dir_okay=False, resolve_path=True),
              help='❌ 失败记录文件路径 | Failed records file path (default: [input]_failed.txt)')
@click.option('--download-script',
              default='download.sh',
              type=click.Path(file_okay=True, dir_okay=False, resolve_path=True),
              help='📜 下载脚本文件名 | Download script filename (default: download.sh)')
@click.option('--ftp-host',
              default='download2.cncb.ac.cn',
              help='🌐 FTP服务器地址 | FTP server host (default: download2.cncb.ac.cn)')
@click.option('--ftp-timeout',
              default=60,
              type=int,
              help='⏱️ FTP连接超时时间(秒) | FTP connection timeout in seconds (default: 60)')
@click.option('--retry-attempts',
              default=3,
              type=int,
              help='🔄 FTP连接重试次数 | FTP connection retry attempts (default: 3)')
@click.option('--threads',
              default=4,
              type=int,
              help='🧵 并发线程数 | Concurrent threads (default: 4)')
@click.option('--verbose', '-v',
              is_flag=True,
              help='📝 详细输出模式 | Verbose output mode')
@click.option('--log-file',
              type=click.Path(file_okay=True, dir_okay=False, resolve_path=True),
              help='📊 日志文件路径 | Log file path')
@click.option('--no-download-script',
              is_flag=True,
              help='⏭️ 不生成下载脚本 | Don\'t generate download script')
@click.option('--no-executable',
              is_flag=True,
              help='🔒 不设置脚本执行权限 | Don\'t make script executable')
def get_link_from_CNCB(input, output, failed, download_script, ftp_host,
                       ftp_timeout, retry_attempts, threads, verbose,
                       log_file, no_download_script, no_executable):
    """
    从CNCB数据库批量获取测序数据下载链接

    这个工具用于从CNCB（国家基因组科学数据中心）的FTP服务器批量查找和获取测序数据的下载链接。
    支持SRA、ERA、DRA数据库的Run ID，并自动生成下载脚本。

    This tool is used to batch search and retrieve sequencing data download links from CNCB (National Genomics Data Center) FTP servers.
    Supports Run IDs from SRA, ERA, DRA databases and automatically generates download scripts.

    📋 输入文件格式 | Input File Format:
    输入文件应为Tab分隔的两列格式，包含ProjectID和RunID：
    The input file should be tab-separated with two columns containing ProjectID and RunID:

    ```
    ProjectID    RunID
    PRJNA123456  SRR12345678
    PRJNA123456  SRR12345679
    PRJNA789012  ERR12345678
    PRJNA789012  DRR12345678
    ```

    🎯 支持的数据库前缀 | Supported Database Prefixes:
    - SRR*: SRA (Sequence Read Archive)
    - ERR*: ERA (European Nucleotide Archive)
    - DRR*: DRA (DNA Data Bank of Japan)

    📁 输出文件 | Output Files:
    - 链接文件 | Links file: 包含所有找到的下载链接
    - 失败记录 | Failed records: 包含未能找到文件的Run ID列表
    - 下载脚本 | Download script: 自动生成的wget下载脚本（可选）

    ⚡ 性能优化 | Performance Optimization:
    - 使用路径缓存减少FTP查询
    - 支持重试机制和超时控制
    - 并发处理提高效率
    - 智能错误处理和恢复

    🔍 搜索策略 | Search Strategy:
    - 智能路径搜索算法
    - 多级目录缓存机制
    - 多文件名模板匹配
    - 自动适应不同的目录结构

    💡 使用技巧 | Usage Tips:
    - 输入文件中可以包含注释行（以#开头）
    - 建议预先筛选Run ID以提高成功率
    - 大规模查询时可以适当增加超时时间
    - 使用详细模式查看处理进度

    📖 应用场景 | Use Cases:
    - 批量下载公共数据库的测序数据
    - 构建本地数据集和数据库
    - 元数据分析和数据管理
    - 生物信息学流水线的数据获取

    示例 | Examples:

    \b
    # 基本用法 | Basic usage
    biopytools get-link-from-CNCB -i projects.txt

    \b
    # 指定输出文件 | Specify output files
    biopytools get-link-from-CNCB -i projects.txt -o links.txt -f failed.txt

    \b
    # 详细输出和日志记录 | Verbose output and logging
    biopytools get-link-from-CNCB -i projects.txt -v --log-file process.log

    \b
    # 自定义下载脚本名称 | Custom download script name
    biopytools get-link-from-CNCB -i projects.txt --download-script my_download.sh

    \b
    # 高级配置 | Advanced configuration
    biopytools get-link-from-CNCB -i large_project_list.txt \\
        -o comprehensive_links.txt \\
        -f comprehensive_failed.txt \\
        --download-script comprehensive_download.sh \\
        --threads 8 \\
        --ftp-timeout 120 \\
        --retry-attempts 5 \\
        -v \\
        --log-file comprehensive.log

    \b
    # 只生成链接文件，不生成下载脚本 | Generate links file only
    biopytools get-link-from-CNCB -i quick_links.txt --no-download-script

    🛠️ 故障排除 | Troubleshooting:
    - 连接超时：增加 --ftp-timeout 参数
    - 权限问题：检查网络和防火墙设置
    - 找不到文件：验证Run ID格式和拼写
    - 内存不足：减少并发线程数

    ⚠️ 注意事项 | Notes:
    - 需要稳定的网络连接
    - 大规模查询可能需要较长时间
    - 建议在网络空闲时运行
    - 尊重数据库的使用政策
    """

    try:
        # 懒加载 | Lazy loading
        CNCLinkExtractor = _lazy_import_cncb_main()

        # 创建提取器 | Create extractor
        extractor = CNCLinkExtractor(
            input_file=str(input),
            output_file=str(output) if output else None,
            failed_file=str(failed) if failed else None,
            download_script=str(download_script),
            ftp_host=ftp_host,
            ftp_timeout=ftp_timeout,
            retry_attempts=retry_attempts,
            max_threads=threads,
            verbose=verbose,
            log_file=str(log_file) if log_file else None,
            generate_download_script=not no_download_script,
            script_executable=not no_executable
        )

        # 运行提取 | Run extraction
        success = extractor.run_extraction()

        if success:
            click.echo("✅ CNCB链接提取完成 | CNCB link extraction completed successfully!")
        else:
            click.echo("❌ CNCB链接提取失败 | CNCB link extraction failed!", err=True)
            sys.exit(1)

    except KeyboardInterrupt:
        click.echo("\n⚠️ 用户中断操作 | Operation interrupted by user", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"💥 发生错误 | Error occurred: {e}", err=True)
        sys.exit(1)