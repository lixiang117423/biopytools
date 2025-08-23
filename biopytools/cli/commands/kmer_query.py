"""
🚀 K-mer提取命令 | K-mer Extraction Command
"""

import click
import sys
# In your actual project, you would use a relative import like this:
# from ...kmer_extractor.main import main as kmer_extractor_main

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
        parser.add_argument('-i', '--input-files', nargs='+', required=True)
        parser.add_argument('-o', '--output-dir', default='./kmer_output')
        parser.add_argument('-k', '--kmer-length', type=int, default=51)
        parser.add_argument('-t', '--threads', type=int, default=88)
        parser.add_argument('-m', '--memory', type=int, default=880)
        parser.add_argument('--file-type', choices=['fasta', 'fastq'])
        parser.add_argument('--fastq-pattern')
        parser.add_argument('--no-canonical', action='store_true')
        parser.add_argument('--no-compress', action='store_true')
        parser.add_argument('--output-bed', action='store_true')
        parser.add_argument('--no-keep-binary', action='store_true')
        parser.add_argument('--unikmer-path', default='unikmer')
        parser.add_argument('--jellyfish-path', default='jellyfish')
        parser.add_argument('--jellyfish-hash-size', default='10000M')
        try:
            args = parser.parse_args()
            print(f"Argparse would have parsed arguments as: {args}")
        except Exception as e:
            print(f"Argparse simulation failed: {e}")
        
        print("--- ✅ Extraction finished (simulated) ---")
    return main_placeholder
kmer_extractor_main = get_original_main_for_demo()
# END: Placeholder

@click.command(context_settings=dict(help_option_names=['-h', '--help']), short_help = "K-mer提取工具")
# --- Required arguments ---
@click.option('--input-files', '-i',
              required=True,
              multiple=True,
              type=click.Path(exists=True, resolve_path=True),
              help='🎯 输入文件或目录 (可多个) | Input files or directories (multiple allowed).')
# --- Output arguments ---
@click.option('--output-dir', '-o',
              default='./kmer_output', show_default=True,
              type=click.Path(file_okay=False, resolve_path=True),
              help='📂 输出目录 | Output directory.')
# --- K-mer parameters ---
@click.option('--kmer-length', '-k',
              type=int, default=51, show_default=True,
              help='🧬 K-mer长度 (1-64) | K-mer length (1-64).')
# --- Performance parameters ---
@click.option('--threads', '-t',
              type=int, default=88, show_default=True,
              help='🚀 线程数 | Number of threads.')
@click.option('--memory', '-m',
              type=int, default=880, show_default=True,
              help='💾 内存限制(GB) | Memory limit (GB).')
# --- File type and pattern ---
@click.option('--file-type',
              type=click.Choice(['fasta', 'fastq'], case_sensitive=False),
              help='📁 文件类型 (默认自动检测) | File type (auto-detects if not specified).')
@click.option('--fastq-pattern',
              help='🔗 FASTQ配对模式 (e.g., "*_1.clean.fq.gz") | FASTQ pairing pattern.')
# --- Processing options ---
@click.option('--no-canonical',
              is_flag=True,
              help='🔄 不使用canonical k-mer | Do not use canonical k-mers.')
@click.option('--no-compress',
              is_flag=True,
              help='🗜️ 不压缩输出文件 | Do not compress output files.')
@click.option('--output-bed',
              is_flag=True,
              help='📋 输出BED文件 (仅FASTA输入) | Output BED file (FASTA input only).')
@click.option('--no-keep-binary',
              is_flag=True,
              help='🗑️ 不保留二进制文件 | Do not keep binary files.')
# --- Tool paths ---
@click.option('--unikmer-path',
              default='unikmer', show_default=True,
              help='⚙️ Unikmer软件路径 | Unikmer software path.')
@click.option('--jellyfish-path',
              default='jellyfish', show_default=True,
              help='🐟 Jellyfish软件路径 | Jellyfish software path.')
# --- Jellyfish-specific parameters ---
@click.option('--jellyfish-hash-size',
              default='10000M', show_default=True,
              help='🗂️ Jellyfish哈希表大小 | Jellyfish hash table size.')
def kmer_query(**kwargs):
    """
    K-mer提取工具.

    从FASTA或FASTQ文件中高效提取唯一的k-mer。
    FASTA输入使用unikmer，支持生成BED坐标；
    FASTQ输入使用jellyfish，支持处理配对文件。
    
    🌟 示例 | Examples:
    
    \b
    # 🎯 从单个FASTA文件提取并生成BED
    biopytools kmer-extractor -i genome.fa -o results --output-bed
    
    \b
    # 🧬 从多个配对FASTQ文件提取
    biopytools kmer-extractor -i ./fastq_dir -o results -k 31 \\
        --fastq-pattern "*_R1.fq.gz"
        
    \b
    # 🚀 高性能运行
    biopytools kmer-extractor -i sample.fq.gz -o out -t 64 -m 500
    """
    
    # 构建参数列表以传递给原始的main函数 🔄 | Build argument list for original main function
    args = ['biopytools', 'kmer-extractor']
    
    # 处理click传递的参数字典
    params = kwargs.copy()
    
    # 特殊处理 nargs='+' 的 input_files
    if 'input_files' in params and params['input_files']:
        args.append('-i')
        args.extend(params.pop('input_files'))
        
    # 遍历其余参数
    for key, value in params.items():
        if value is None:
            continue
            
        # 将Python风格的变量名转为命令行风格
        param_name = '--' + key.replace('_', '-')
        
        # 处理布尔标志
        if isinstance(value, bool) and value:
            args.append(param_name)
        # 处理非布尔、非默认值的常规参数
        elif not isinstance(value, bool):
            # 获取参数的默认值
            default_val = kmer_extractor.params_by_name[key].default
            # 只有当值不等于默认值时才添加
            if value != default_val:
                args.append(param_name)
                args.append(str(value))

    # 保存并恢复sys.argv 💾 | Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args
    
    try:
        # 调用原始的main函数 🚀 | Call original main function
        kmer_extractor_main()
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
    kmer_query()