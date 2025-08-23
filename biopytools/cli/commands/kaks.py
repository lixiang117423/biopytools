"""
🏭 Ka/Ks Calculator Command
"""

import click
import sys

# 1. 使用绝对路径导入原始模块的 main 函数
#    假设您的项目根目录是 biopytools
from biopytools.kakscalc.main import main as kaks_main

# 2. 使用绝对路径导入原始模块的 config 文件
#    这样它就能准确找到 KaKsConfig 类
from biopytools.kakscalc.config import KaKsConfig

# --- Placeholder for the original main function to make this snippet runnable ---
# --- In your project, delete this section and use the import statement above ---
# START: Placeholder for demonstration
def get_original_main_for_demo():
    def main_placeholder():
        print("--- 🚀 Original main function called (simulated) ---")
        print(f"Received sys.argv: {sys.argv}")
        # The original argparse would parse sys.argv and run the KaKsAnalyzer.
        import argparse
        parser = argparse.ArgumentParser()
        parser.add_argument('-1', '--fasta1', required=True)
        parser.add_argument('-2', '--fasta2', required=True)
        parser.add_argument('-p', '--pairs', required=True)
        parser.add_argument('-o', '--output', required=True)
        parser.add_argument('-m', '--method', default=KaKsConfig.DEFAULT_METHOD, choices=KaKsConfig.SUPPORTED_METHODS)
        parser.add_argument('--kaks-path', default='KaKs_Calculator')
        parser.add_argument('-v', '--verbose', action='store_true')
        parser.add_argument('--temp-dir')
        parser.add_argument('--keep-temp', action='store_true')
        try:
            args = parser.parse_args()
            print(f"Argparse would have parsed arguments as: {args}")
        except Exception as e:
            print(f"Argparse simulation failed: {e}")
        
        print("--- ✅ Analysis finished (simulated) ---")
    return main_placeholder
kaks_main = get_original_main_for_demo()
# END: Placeholder

@click.command(context_settings=dict(help_option_names=['-h', '--help']), short_help = "计算Ka/Ks")
@click.version_option(version='1.0.0', prog_name='Ka/Ks Calculator')
# --- Required arguments ---
@click.option('--fasta1', '-1',
              required=True,
              type=click.Path(exists=True, dir_okay=False, resolve_path=True),
              help='🧬 第一个FASTA文件 (物种1 CDS序列) | First FASTA file (species 1 CDS).')
@click.option('--fasta2', '-2',
              required=True,
              type=click.Path(exists=True, dir_okay=False, resolve_path=True),
              help='🧬 第二个FASTA文件 (物种2 CDS序列) | Second FASTA file (species 2 CDS).')
@click.option('--pairs', '-p',
              required=True,
              type=click.Path(exists=True, dir_okay=False, resolve_path=True),
              help='🔗 序列配对文件 (TSV/CSV格式) | Sequence pair file (TSV/CSV format).')
@click.option('--output', '-o',
              required=True,
              type=click.Path(file_okay=False, resolve_path=True),
              help='📂 输出目录 | Output directory for results.')
# --- Analysis parameters ---
@click.option('--method', '-m',
              type=click.Choice(KaKsConfig.SUPPORTED_METHODS, case_sensitive=False),
              default=KaKsConfig.DEFAULT_METHOD, show_default=True,
              help='🧮 计算方法 | Calculation method.')
@click.option('--kaks-path',
              default='KaKs_Calculator', show_default=True,
              help='🛠️ KaKs_Calculator可执行文件路径 | Path to KaKs_Calculator executable.')
# --- Runtime options ---
@click.option('--verbose', '-v',
              is_flag=True,
              help='🔍 启用详细日志记录 | Enable verbose logging.')
@click.option('--temp-dir',
              type=click.Path(file_okay=False, resolve_path=True),
              help='🗂️ 自定义临时目录 | Custom temporary directory.')
@click.option('--keep-temp',
              is_flag=True,
              help='🗑️ 保留临时文件 (用于调试) | Keep temporary files (for debugging).')
def kaks(fasta1, fasta2, pairs, output, method, kaks_path, verbose, temp_dir, keep_temp):
    """
    计算Ka/Ks.

    一个用于批量计算编码序列对之间非同义替换率(Ka)和同义替换率(Ks)的工具，
    以推断分子演化中的选择压力。
    
    📋 示例 | Examples:
    
    \b
    # 🎯 基本用法
    biopytools kaks-calculator -1 sp1.fa -2 sp2.fa -p pairs.txt -o ./results
    
    \b
    # 🔧 指定计算方法并启用详细日志
    biopytools kaks-calculator -1 human.fa -2 mouse.fa -p orthologs.tsv \\
        -o ./analysis -m YN -v
    
    \b
    # 🛠️ 使用自定义的KaKs_Calculator路径
    biopytools kaks-calculator -1 cds1.fa -2 cds2.fa -p pairs.txt \\
        -o ./out --kaks-path /opt/kaks/KaKs_Calculator
    """
    
    # 构建参数列表以传递给原始的main函数 🔄 | Build argument list for original main function
    args = ['biopytools', 'kaks-calculator']
    
    # 必需参数 📌 | Required parameters
    args.extend(['-1', fasta1])
    args.extend(['-2', fasta2])
    args.extend(['-p', pairs])
    args.extend(['-o', output])
    
    # 可选参数 ⚙️ | Optional parameters
    if method.upper() != KaKsConfig.DEFAULT_METHOD:
        args.extend(['-m', method])
    
    if kaks_path != 'KaKs_Calculator':
        args.extend(['--kaks-path', kaks_path])
        
    if verbose:
        args.append('-v')
        
    if temp_dir:
        args.extend(['--temp-dir', temp_dir])
        
    if keep_temp:
        args.append('--keep-temp')
        
    # 保存并恢复sys.argv 💾 | Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args
    
    try:
        # 调用原始的main函数 🚀 | Call original main function
        kaks_main()
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
    kaks()