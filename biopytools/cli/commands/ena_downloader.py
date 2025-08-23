"""
📥 ENA下载工具命令 | ENA Downloader Command
"""

import click
import sys
# In your actual project, you would use a relative import like this:
# from ...ena_downloader.main import main as ena_downloader_main

# --- Placeholder for the original main function to make this snippet runnable ---
# --- In your project, delete this section and use the import statement above ---
# START: Placeholder for demonstration
def get_original_main_for_demo():
    def main_placeholder():
        print("--- 🚀 Original main function called (simulated) ---")
        print(f"Received sys.argv: {sys.argv}")
        # The original main() would create a parser and run the analysis.
        import argparse
        parser = argparse.ArgumentParser()
        parser.add_argument('--accession', '-a', required=True)
        parser.add_argument('--output-dir', '-o')
        parser.add_argument('--create-dir', '-d', action='store_true')
        parser.add_argument('--metadata-format', '-f', choices=['tsv', 'csv', 'xlsx'], default='tsv')
        parser.add_argument('--protocol', '-p', choices=['ftp', 'aspera'], default='ftp')
        parser.add_argument('--aspera-key', '-k')
        parser.add_argument('--method', '-m', choices=['save', 'run'], default='save')
        parser.add_argument('--metadata-only', '-M', action='store_true')
        parser.add_argument('--fields', '-F', nargs='+')
        parser.add_argument('--max-retries', '-r', type=int, default=3)
        try:
            args = parser.parse_args()
            print(f"Argparse would have parsed arguments as: {args}")
        except Exception as e:
            print(f"Argparse simulation failed: {e}")
        
        print("--- ✅ Process finished (simulated) ---")
    return main_placeholder
ena_downloader_main = get_original_main_for_demo()
# END: Placeholder

@click.command(context_settings=dict(help_option_names=['-h', '--help']), short_help = " ENA样品信息和下载链接获取工具")
# --- Required arguments ---
@click.option('--accession', '-a',
              required=True,
              help='🎯 ENA项目编号 (e.g., PRJNA661210) | ENA accession number.')
# --- Output settings ---
@click.option('--output-dir', '-o',
              type=click.Path(file_okay=False, resolve_path=True),
              help='📁 输出目录 (默认: 当前目录) | Output directory (default: current).')
@click.option('--create-dir', '-d',
              is_flag=True,
              help='📂 创建专门的输出目录 ([accession].ena.download) | Create a dedicated output directory.')
@click.option('--metadata-format', '-f',
              type=click.Choice(['tsv', 'csv', 'xlsx'], case_sensitive=False),
              default='tsv', show_default=True,
              help='📋 元数据文件格式 | Metadata file format.')
# --- Download protocol settings ---
@click.option('--protocol', '-p',
              type=click.Choice(['ftp', 'aspera'], case_sensitive=False),
              default='ftp', show_default=True,
              help='🌐 下载协议 (ftp 或 aspera) | Download protocol.')
@click.option('--aspera-key', '-k',
              type=click.Path(exists=True, dir_okay=False, resolve_path=True),
              help='🔐 Aspera私钥路径 (使用aspera时必需) | Path to Aspera private key.')
@click.option('--method', '-m',
              type=click.Choice(['save', 'run'], case_sensitive=False),
              default='save', show_default=True,
              help='⚙️ 执行模式 (save: 生成脚本, run: 直接下载) | Execution mode.')
# --- Special modes ---
@click.option('--metadata-only', '-M',
              is_flag=True,
              help='📊 仅下载元数据 | Only download metadata.')
# --- Advanced options ---
@click.option('--fields', '-F',
              multiple=True,
              help='🔧 自定义元数据字段 (可多次使用) | Custom metadata fields (use multiple times).')
@click.option('--max-retries', '-r',
              type=int, default=3, show_default=True,
              help='🔄 API请求最大重试次数 | Maximum API request retries.')
def ena_downloader(**kwargs):
    """
    ENA样品信息和下载链接获取工具.

    根据ENA项目编号(Accession)下载相关的元数据和FASTQ文件。
    支持FTP和高速Aspera协议，可以选择直接下载或生成下载脚本。
    
    🌟 示例 | Examples:
    
    \b
    # 📊 仅下载元数据
    biopytools ena-downloader -a PRJNA123456 -M
    
    \b
    # 💾 生成FTP下载脚本
    biopytools ena-downloader -a PRJNA123456 -p ftp -m save
        
    \b
    # 🚀 使用Aspera直接下载
    biopytools ena-downloader -a PRJNA123456 -p aspera -k ~/.aspera/key.openssh -m run
    """
    
    # 构建参数列表以传递给原始的main函数 🔄 | Build argument list for original main function
    args = ['biopytools', 'ena-downloader']
    
    # 遍历所有参数
    for key, value in kwargs.items():
        if value is None:
            continue
            
        param_name = '--' + key.replace('_', '-')
        
        if isinstance(value, bool) and value:
            args.append(param_name)
        elif isinstance(value, tuple) and value: # For 'multiple=True' options like --fields
            args.append(param_name)
            args.extend(value)
        elif not isinstance(value, (bool, tuple)):
            default_val = ena_downloader.params_by_name[key].default
            if value != default_val:
                args.append(param_name)
                args.append(str(value))
    
    # 确保必需参数和有默认值的参数总是存在（如果它们是核心功能的一部分）
    # Required
    args.extend(['-a', kwargs['accession']])
    # Defaults that are always passed to the underlying script
    if 'protocol' not in args: args.extend(['-p', kwargs['protocol']])
    if 'method' not in args: args.extend(['-m', kwargs['method']])

    # 保存并恢复sys.argv 💾 | Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args
    
    try:
        # 调用原始的main函数 🚀 | Call original main function
        ena_downloader_main()
    except SystemExit as e:
        # 处理程序正常退出 ✅ | Handle normal program exit
        if e.code != 0:
            click.secho(f"❌ 脚本执行被终止，退出码: {e.code}", fg='red', err=True)
        sys.exit(e.code)
    except Exception as e:
        click.secho(f"💥 发生未知错误 | An unexpected error occurred: {e}", fg='red', err=True)
        sys.exit(1)
    finally:
        # 无论如何都要恢复原始的 sys.argv | Restore original sys.argv regardless of outcome
        sys.argv = original_argv

# 如果直接运行此文件用于测试 | If running this file directly for testing
if __name__ == '__main__':
    ena_downloader()