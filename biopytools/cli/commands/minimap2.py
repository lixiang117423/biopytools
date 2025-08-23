"""
🧬 Minimap2比对与分析命令 | Minimap2 Alignment and Analysis Command
"""

import click
import sys
# In your actual project, you would use a relative import like this:
# from ...minimap2.main import main as minimap2_main

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
        parser.add_argument('-t', '--target', required=True)
        parser.add_argument('-q', '--query', required=True)
        parser.add_argument('-o', '--output-dir', default='./minimap2_output')
        parser.add_argument('-x', '--preset', default='asm5', choices=['asm5', 'asm10', 'asm20', 'map-ont', 'map-pb'])
        parser.add_argument('-p', '--threads', type=int, default=8)
        parser.add_argument('-m', '--min-match', type=int, default=1000)
        parser.add_argument('-u', '--min-unmapped', type=int, default=1000)
        parser.add_argument('--tp-type', default='P', choices=['S', 'P', 'SP'])
        parser.add_argument('-M', '--minimap2-path', default='minimap2')
        parser.add_argument('-S', '--seqkit-path', default='seqkit')
        try:
            args = parser.parse_args()
            print(f"Argparse would have parsed arguments as: {args}")
        except Exception as e:
            print(f"Argparse simulation failed: {e}")
        
        print("--- ✅ Analysis finished (simulated) ---")
    return main_placeholder
minimap2_main = get_original_main_for_demo()
# END: Placeholder

@click.command(context_settings=dict(help_option_names=['-h', '--help']), short_help = "Minimap2全基因组比对和未比对区间提取工具")
# --- Required arguments ---
@click.option('--target', '-t',
              required=True,
              type=click.Path(exists=True, dir_okay=False, resolve_path=True),
              help='🎯 目标基因组文件路径 | Target genome file path.')
@click.option('--query', '-q',
              required=True,
              type=click.Path(exists=True, dir_okay=False, resolve_path=True),
              help='🧬 查询基因组文件路径 | Query genome file path.')
# --- Optional arguments ---
@click.option('--output-dir', '-o',
              default='./minimap2_output', show_default=True,
              type=click.Path(file_okay=False, resolve_path=True),
              help='📂 输出目录 | Output directory.')
@click.option('--preset', '-x',
              type=click.Choice(['asm5', 'asm10', 'asm20', 'map-ont', 'map-pb'], case_sensitive=False),
              default='asm5', show_default=True,
              help='🔧 Minimap2预设参数 | Minimap2 preset parameters.')
@click.option('--threads', '-p',
              type=int, default=8, show_default=True,
              help='🚀 线程数 | Number of threads.')
@click.option('--min-match', '-m',
              type=int, default=1000, show_default=True,
              help='📏 最小匹配长度阈值 | Minimum match length threshold.')
@click.option('--min-unmapped', '-u',
              type=int, default=1000, show_default=True,
              help='📏 最小未比对区间长度阈值 | Minimum unmapped region length threshold.')
@click.option('--tp-type',
              type=click.Choice(['S', 'P', 'SP'], case_sensitive=False),
              default='P', show_default=True,
              help='🔗 保留的比对类型 (P: primary) | Alignment type to keep (P: primary).')
@click.option('--minimap2-path', '-M',
              default='minimap2', show_default=True,
              help='⚙️ minimap2可执行文件路径 | minimap2 executable path.')
@click.option('--seqkit-path', '-S',
              default='seqkit', show_default=True,
              help='⚙️ seqkit可执行文件路径 | seqkit executable path.')
def minimap2(target, query, output_dir, preset, threads, min_match, min_unmapped, tp_type, minimap2_path, seqkit_path):
    """
    Minimap2全基因组比对和未比对区间提取工具.

    该工具使用Minimap2进行全基因组比对，然后解析PAF结果，
    识别并提取在查询(query)基因组中相对于目标(target)基因组的未比对区域。
    
    🌟 示例 | Examples:
    
    \b
    # 🎯 基本用法
    biopytools minimap2 -t target.fa -q query.fa -o ./results
    
    \b
    # 🔧 使用不同的预设和过滤参数
    biopytools minimap2 -t ref.fa -q sample.fa -o ./analysis \\
        -x asm10 -m 5000 -u 2000
        
    \b
    # 🚀 使用更多线程
    biopytools minimap2 -t target.fa -q query.fa -o ./out -p 32
    """
    
    # 构建参数列表以传递给原始的main函数 🔄 | Build argument list for original main function
    args = ['biopytools', 'minimap2']
    
    # 必需参数 📋 | Required parameters
    args.extend(['-t', target])
    args.extend(['-q', query])
    
    # 可选参数 ⚙️ | Optional parameters
    if output_dir != './minimap2_output':
        args.extend(['-o', output_dir])
    else: # 总是传递输出目录
        args.extend(['-o', output_dir])
        
    if preset != 'asm5':
        args.extend(['-x', preset])
        
    if threads != 8:
        args.extend(['-p', str(threads)])
        
    if min_match != 1000:
        args.extend(['-m', str(min_match)])
        
    if min_unmapped != 1000:
        args.extend(['-u', str(min_unmapped)])
        
    if tp_type.upper() != 'P':
        args.extend(['--tp-type', tp_type])
        
    if minimap2_path != 'minimap2':
        args.extend(['-M', minimap2_path])
        
    if seqkit_path != 'seqkit':
        args.extend(['-S', seqkit_path])
        
    # 保存并恢复sys.argv 💾 | Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args
    
    try:
        # 调用原始的main函数 🚀 | Call original main function
        minimap2_main()
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
    minimap2()