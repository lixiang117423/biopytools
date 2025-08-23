"""
🧬 FASTA ID分割命令 | FASTA ID Splitting Command
"""

import click
import sys
# In your actual project, you would use a relative import like this:
# from ...split_fasta_id.main import main as split_fasta_id_main

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
        parser.add_argument('-i', '--input', required=True)
        parser.add_argument('-o', '--output', default='output.fasta')
        parser.add_argument('-p', '--position', type=int, default=0)
        parser.add_argument('-d', '--delimiter', default='auto')
        parser.add_argument('--keep-original', action='store_true')
        parser.add_argument('--no-skip-empty', action='store_true')
        parser.add_argument('--preserve-comments', action='store_true')
        try:
            args = parser.parse_args()
            print(f"Argparse would have parsed arguments as: {args}")
        except Exception as e:
            print(f"Argparse simulation failed: {e}")
        
        print("--- ✅ Splitting finished (simulated) ---")
    return main_placeholder
split_fasta_id_main = get_original_main_for_demo()
# END: Placeholder

@click.command(context_settings=dict(help_option_names=['-h', '--help']))
# --- Required arguments ---
@click.option('--input', '-i', 'input_file', # Rename to avoid conflict with Python's 'input'
              required=True,
              type=click.Path(exists=True, dir_okay=False, resolve_path=True),
              help='📂 输入FASTA文件路径 | Input FASTA file path.')
# --- Optional arguments ---
@click.option('--output', '-o', 'output_file',
              default='output.fasta', show_default=True,
              type=click.Path(dir_okay=False, resolve_path=True),
              help='📁 输出FASTA文件路径 | Output FASTA file path.')
@click.option('--position', '-p',
              type=int, default=0, show_default=True,
              help='🎯 提取位置 (0表示第一个元素) | Extract position (0 means first element).')
@click.option('--delimiter', '-d',
              default='auto', show_default=True,
              help='🔧 分隔符 (auto, space, tab, or a character) | Delimiter.')
# --- Processing options ---
@click.option('--keep-original',
              is_flag=True,
              help='💾 保留原始文件作为备份 | Keep original file as backup.')
@click.option('--no-skip-empty',
              is_flag=True,
              help='⚠️ 不跳过空的序列名称行 | Do not skip empty sequence name lines.')
@click.option('--preserve-comments',
              is_flag=True,
              help='📝 保留序列ID中的注释部分 | Preserve comments in sequence IDs.')
def split_fasta_id(input_file, output_file, position, delimiter, keep_original, no_skip_empty, preserve_comments):
    """
    FASTA序列ID分割工具.

    根据指定的分隔符分割FASTA文件中的序列ID，并提取特定位置的子字符串作为新的ID。
    
    🌟 示例 | Examples:
    
    \b
    # 🎯 提取ID的第一部分 (默认以空格/tab分割)
    biopytools split-fasta-id -i seqs.fa -o clean.fa -p 0
    
    \b
    # 🔧 使用'|'作为分隔符，提取第二部分
    biopytools split-fasta-id -i db.fa -o simple.fa -p 1 -d "|"
        
    \b
    # 💾 保留原始文件作为备份
    biopytools split-fasta-id -i data.fa -o out.fa --keep-original
    """
    
    # 构建参数列表以传递给原始的main函数 🔄 | Build argument list for original main function
    args = ['biopytools', 'split-fasta-id']
    
    # 必需参数 📋 | Required parameters
    args.extend(['-i', input_file])
    
    # 可选参数 ⚙️ | Optional parameters
    if output_file != 'output.fasta':
        args.extend(['-o', output_file])
    else: # Always pass output for simplicity, argparse will handle the default
        args.extend(['-o', output_file])
        
    if position != 0:
        args.extend(['-p', str(position)])
        
    if delimiter != 'auto':
        args.extend(['-d', delimiter])
        
    if keep_original:
        args.append('--keep-original')
        
    if no_skip_empty:
        args.append('--no-skip-empty')
        
    if preserve_comments:
        args.append('--preserve-comments')
        
    # 保存并恢复sys.argv 💾 | Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args
    
    try:
        # 调用原始的main函数 🚀 | Call original main function
        split_fasta_id_main()
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
    split_fasta_id()