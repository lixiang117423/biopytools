"""
🚀 K-mer丰度分析命令 | K-mer Abundance Analysis Command
"""

import click
import sys
# In your actual project, you would use a relative import like this:
# from ...kmer_count.main import main as kmer_count_main

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
        parser.add_argument('-p', '--pattern', required=True)
        parser.add_argument('-k', '--kmer-lib', required=True)
        parser.add_argument('-o', '--output', required=True)
        parser.add_argument('-b', '--bed-file')
        parser.add_argument('-m', '--kmer-size', type=int, default=51)
        parser.add_argument('-s', '--hash-size', default='1000M')
        parser.add_argument('-t', '--threads', type=int, default=8)
        parser.add_argument('-w', '--window-size', type=int, default=500000)
        parser.add_argument('--step-size', type=int)
        parser.add_argument('-C', '--canonical', action='store_true')
        parser.add_argument('--keep-temp', action='store_true')
        parser.add_argument('--keep-binary', action='store_true')
        parser.add_argument('--jellyfish-path', default='jellyfish')
        parser.add_argument('-v', '--verbose', action='store_true')
        try:
            args = parser.parse_args()
            print(f"Argparse would have parsed arguments as: {args}")
        except Exception as e:
            print(f"Argparse simulation failed: {e}")
        
        print("--- ✅ Analysis finished (simulated) ---")
    return main_placeholder
kmer_count_main = get_original_main_for_demo()
# END: Placeholder

@click.command(context_settings=dict(help_option_names=['-h', '--help']))
# --- Required arguments ---
@click.option('--input', '-i',
              required=True,
              type=click.Path(exists=True, file_okay=False, dir_okay=True, resolve_path=True),
              help='📂 输入文件目录 (FASTQ/FASTA) | Input files directory (FASTQ/FASTA).')
@click.option('--pattern', '-p',
              required=True,
              help='📁 文件模式 (e.g., "*_1.fq.gz", "*.fa") | File pattern.')
@click.option('--kmer-lib', '-k',
              required=True,
              type=click.Path(exists=True, dir_okay=False, resolve_path=True),
              help='🧬 K-mer库文件 (FASTA格式) | K-mer library file (FASTA format).')
@click.option('--output', '-o',
              required=True,
              type=click.Path(file_okay=False, resolve_path=True),
              help='📂 输出目录 | Output directory.')
# --- Optional arguments ---
@click.option('--bed-file', '-b',
              type=click.Path(exists=True, dir_okay=False, resolve_path=True),
              help='📋 BED文件路径 (用于窗口分析) | BED file path (for window analysis).')
@click.option('--kmer-size', '-m',
              type=int, default=51, show_default=True,
              help='📏 K-mer长度 | K-mer size.')
@click.option('--hash-size', '-s',
              default='1000M', show_default=True,
              help='🗂️ Jellyfish哈希表大小 | Jellyfish hash table size.')
@click.option('--threads', '-t',
              type=int, default=8, show_default=True,
              help='🧵 线程数 | Number of threads.')
@click.option('--window-size', '-w',
              type=int, default=500000, show_default=True,
              help='🪟 滑动窗口大小(bp) | Sliding window size (bp).')
@click.option('--step-size',
              type=int,
              help='👣 滑动窗口步长(bp) (默认: window-size/5) | Sliding window step size (bp).')
@click.option('--canonical', '-C',
              is_flag=True,
              help='🔄 统计正向和反向互补链 | Count both forward and reverse complement.')
@click.option('--keep-temp',
              is_flag=True,
              help='💾 保留临时文件 | Keep temporary files.')
@click.option('--keep-binary',
              is_flag=True,
              help='🔢 保留0/1存在/缺失矩阵 | Keep 0/1 presence/absence matrix.')
@click.option('--jellyfish-path',
              default='jellyfish', show_default=True,
              help='🐙 Jellyfish程序路径 | Jellyfish program path.')
@click.option('--verbose', '-v',
              is_flag=True,
              help='📝 详细输出 | Verbose output.')
def kmer_count(input, pattern, kmer_lib, output, bed_file, kmer_size, hash_size,
               threads, window_size, step_size, canonical, keep_temp, keep_binary,
               jellyfish_path, verbose):
    """
    K-mer丰度分析工具.

    使用Jellyfish计算指定k-mer文库在多个样本中的丰度，并生成丰度矩阵。
    
    🌟 示例 | Examples:
    
    \b
    # 🎯 基本用法 (FASTQ输入)
    biopytools kmer-count -i ./fastq -p "*_1.fq.gz" -k lib.fa -o results
    
    \b
    # 🪟 带滑动窗口分析 (FASTA输入)
    biopytools kmer-count -i ./fasta -p "*.fa" -k lib.fa -b lib.bed \\
        -w 250000 -o results_window
        
    \b
    # ⚡️ 使用更多线程和自定义Jellyfish路径
    biopytools kmer-count -i ./data -p "*_R1.fastq" -k lib.fa -o out \\
        -t 32 --jellyfish-path /opt/jellyfish/bin/jellyfish
    """
    
    # 构建参数列表以传递给原始的main函数 🔄 | Build argument list for original main function
    args = ['biopytools', 'kmer-count']
    
    # 必需参数 📋 | Required parameters
    args.extend(['-i', input])
    args.extend(['-p', pattern])
    args.extend(['-k', kmer_lib])
    args.extend(['-o', output])
    
    # 可选参数 ⚙️ | Optional parameters
    if bed_file:
        args.extend(['-b', bed_file])
    if kmer_size != 51:
        args.extend(['-m', str(kmer_size)])
    if hash_size != '1000M':
        args.extend(['-s', hash_size])
    if threads != 8:
        args.extend(['-t', str(threads)])
    if window_size != 500000:
        args.extend(['-w', str(window_size)])
    if step_size:
        args.extend(['--step-size', str(step_size)])
    if canonical:
        args.append('-C')
    if keep_temp:
        args.append('--keep-temp')
    if keep_binary:
        args.append('--keep-binary')
    if jellyfish_path != 'jellyfish':
        args.extend(['--jellyfish-path', jellyfish_path])
    if verbose:
        args.append('-v')

    # 保存并恢复sys.argv 💾 | Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args
    
    try:
        # 调用原始的main函数 🚀 | Call original main function
        kmer_count_main()
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
    kmer_count()