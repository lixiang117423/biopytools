"""
🧬 Bismark甲基化分析命令 | Bismark Methylation Analysis Command
"""

import click
import sys
# 在您的实际项目中，您应该使用这样的相对导入
# from ...bismark.main import main as bismark_main
# from ...bismark import __version__

# --- 为了让此代码块能独立运行，我们在此处包含占位符 ---
# --- 在您的项目中，您应该使用上面的 import 语句，并删除这部分 ---
# START: Placeholder for demonstration
def get_original_main_for_demo():
    def main_placeholder():
        print("--- 🚀 Original main function called (simulated) ---")
        print(f"Received sys.argv: {sys.argv}")
        # 原始的 main() 会创建解析器并运行分析
        import argparse
        parser = argparse.ArgumentParser()
        parser.add_argument('-g', '--genome-fa', required=True)
        parser.add_argument('-r', '--raw-dir', required=True)
        parser.add_argument('-o', '--output-dir', required=True)
        parser.add_argument('-p', '--pattern', default='_1_clean.fq.gz')
        parser.add_argument('-j', '--threads', type=int, default=88)
        parser.add_argument('--sort-buffer', type=str, default='400G')
        parser.add_argument('--no-no-overlap', dest='no_overlap', action='store_false')
        try:
            args = parser.parse_args()
            print(f"Argparse would have parsed arguments as: {args}")
        except Exception as e:
            print(f"Argparse simulation failed: {e}")
        
        print("--- ✅ Analysis finished (simulated) ---")
    return main_placeholder
bismark_main = get_original_main_for_demo()
__version__ = "1.0.0" # Placeholder version
# END: Placeholder

@click.command(context_settings=dict(help_option_names=['-h', '--help']), short_help = "Bismark甲基化分析流程")
@click.version_option(version=__version__, prog_name='Bismark Pipeline')
# --- Required arguments ---
@click.option('--genome-fa', '-g',
              required=True,
              type=click.Path(exists=True, dir_okay=False, resolve_path=True),
              help='🧬 基因组FASTA文件路径 | Path to genome FASTA file.')
@click.option('--raw-dir', '-r',
              required=True,
              type=click.Path(exists=True, file_okay=False, dir_okay=True, resolve_path=True),
              help='📁 原始FASTQ数据目录 | Raw FASTQ data directory.')
@click.option('--output-dir', '-o',
              required=True,
              type=click.Path(file_okay=False, resolve_path=True),
              help='📂 主输出目录 | Main output directory.')
# --- Optional arguments ---
@click.option('--pattern', '-p',
              default='_1_clean.fq.gz', show_default=True,
              help='🔗 R1文件的后缀模式 | Suffix pattern for R1 files.')
@click.option('--threads', '-j',
              type=int, default=88, show_default=True,
              help='🚀 使用的线程数 | Number of threads to use.')
@click.option('--sort-buffer',
              type=str, default='400G', show_default=True,
              help='💾 提取步骤的排序缓存大小 | Sort buffer size for extraction step.')
@click.option('--include-overlap',
              is_flag=True,
              help='🧪 在提取时【包含】重叠的reads (默认忽略) | Include overlapping reads during extraction.')
def bismark(genome_fa, raw_dir, output_dir, pattern, threads, sort_buffer, include_overlap):
    """🧬 Bismark甲基化分析流程。
    
    一个自动化运行Bismark进行全基因组重亚硫酸盐测序(WGBS)数据分析的流程，
    包括索引构建、比对和甲基化提取。
    
    示例 | Examples:
    
    \b
    # 🎯 基础用法
    biopytools bismark -g genome.fa -r ./raw_data -o ./bismark_results
    
    \b
    # 🚀 使用更多线程和自定义模式
    biopytools bismark -g hg38.fa -r ./cleandata -o ./results -j 96 -p "_1.fq.gz"
        
    \b
    # 🧪 包含重叠的reads
    biopytools bismark -g genome.fa -r ./data -o ./out --include-overlap
    """
    
    # 构建参数列表以传递给原始的main函数 🔄 | Build argument list for original main function
    args = ['biopytools', 'bismark']
    
    # 必需参数 📋 | Required parameters
    args.extend(['-g', genome_fa])
    args.extend(['-r', raw_dir])
    args.extend(['-o', output_dir])
    
    # 可选参数 ⚙️ | Optional parameters
    if pattern != '_1_clean.fq.gz':
        args.extend(['-p', pattern])
    if threads != 88:
        args.extend(['-j', str(threads)])
    if sort_buffer != '400G':
        args.extend(['--sort-buffer', sort_buffer])
        
    # 处理反向布尔标志 | Handle inverted boolean flag
    # argparse的 '--no-no-overlap' (action='store_false') 意味着默认 no_overlap=True (忽略重叠)
    # 当用户在click中指定 --include-overlap 时，我们必须传递 --no-no-overlap 来让 no_overlap 变为 False
    if include_overlap:
        args.append('--no-no-overlap')

    # 保存并恢复sys.argv 💾 | Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args
    
    try:
        # 调用原始的main函数 🚀 | Call original main function
        bismark_main()
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
    bismark()