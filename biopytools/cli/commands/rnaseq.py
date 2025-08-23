"""
🧬 RNA-seq分析命令 | RNA-seq Analysis Command
"""

import click
import sys
# In your actual project, you would use a relative import like this:
# from ...rnaseq.main import main as rnaseq_main

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
        parser.add_argument("-g", "--genome", required=True)
        parser.add_argument("-f", "--gtf", required=True)
        parser.add_argument("-i", "--input", required=True)
        parser.add_argument("-o", "--output", required=True)
        parser.add_argument("-p", "--pattern", default=None)
        parser.add_argument("-r", "--remove", default="no", choices=["yes", "y", "no", "n"])
        parser.add_argument("-t", "--threads", type=int, default=8)
        try:
            args = parser.parse_args()
            print(f"Argparse would have parsed arguments as: {args}")
        except Exception as e:
            print(f"Argparse simulation failed: {e}")
        
        print("--- ✅ Analysis finished (simulated) ---")
    return main_placeholder
rnaseq_main = get_original_main_for_demo()
# END: Placeholder

@click.command(context_settings=dict(help_option_names=['-h', '--help']), short_help = "RNA-seq分析流程: HISAT2 + StringTie")
# --- Required arguments ---
@click.option('--genome', '-g',
              required=True,
              type=click.Path(exists=True, dir_okay=False, resolve_path=True),
              help='🧬 基因组FASTA文件路径 | Genome FASTA file path.')
@click.option('--gtf', '-f',
              required=True,
              type=click.Path(exists=True, dir_okay=False, resolve_path=True),
              help='📜 基因注释GTF文件路径 | Gene annotation GTF file path.')
@click.option('--input', '-i', 'input_path', # 'input' is a keyword in Python, so we rename it
              required=True,
              type=click.Path(exists=True, resolve_path=True),
              help='📂 输入FASTQ目录或样本信息文件 | Input FASTQ directory or sample info file.')
@click.option('--output', '-o',
              required=True,
              type=click.Path(file_okay=False, resolve_path=True),
              help='📤 输出目录 | Output directory.')
# --- Optional arguments ---
@click.option('--pattern', '-p',
              help='🔗 FASTQ文件命名模式 (e.g., "*_1.fq.gz") | FASTQ file naming pattern.')
@click.option('--remove', '-r',
              type=click.Choice(['yes', 'y', 'no', 'n'], case_sensitive=False),
              default='no', show_default=True,
              help='🗑️ 处理后删除BAM文件 | Remove BAM files after processing.')
@click.option('--threads', '-t',
              type=int, default=8, show_default=True,
              help='🚀 线程数 | Number of threads.')
def rnaseq(genome, gtf, input_path, output, pattern, remove, threads):
    """
    RNA-seq分析流程: HISAT2 + StringTie.

    一个标准的RNA-seq分析流程，用于从FASTQ文件生成基因表达矩阵。
    它会自动构建索引、比对、定量，并合并所有样本的结果。
    
    🌟 示例 | Examples:
    
    \b
    # 🎯 基本用法 (自动发现配对文件)
    biopytools rnaseq -g genome.fa -f genes.gtf -i ./fastq_dir \\
        -o ./results -p "*_R1.fastq.gz"
    
    \b
    # 📁 使用样本信息文件作为输入
    biopytools rnaseq -g genome.fa -f genes.gtf -i samples.tsv -o ./results
        
    \b
    # 🚀 使用更多线程并保留BAM文件
    biopytools rnaseq -g genome.fa -f genes.gtf -i ./fastq_dir \\
        -o ./out -t 32 -r no
    """
    
    # 构建参数列表以传递给原始的main函数 🔄 | Build argument list for original main function
    args = ['biopytools', 'rnaseq']
    
    # 必需参数 📋 | Required parameters
    args.extend(['-g', genome])
    args.extend(['-f', gtf])
    args.extend(['-i', input_path])
    args.extend(['-o', output])
    
    # 可选参数 ⚙️ | Optional parameters
    if pattern:
        args.extend(['-p', pattern])
        
    if remove.lower() not in ['no', 'n']:
        args.extend(['-r', remove])
        
    if threads != 8:
        args.extend(['-t', str(threads)])

    # 保存并恢复sys.argv 💾 | Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args
    
    try:
        # 调用原始的main函数 🚀 | Call original main function
        rnaseq_main()
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
    rnaseq()