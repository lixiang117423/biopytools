"""
🧬 最长转录本提取命令 | Longest mRNA Extraction Command
"""

import click
import sys
# In your actual project, you would use a relative import like this:
# from ...longest_mrna.main import main as longest_mrna_main

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
        parser.add_argument('-g', '--genome', required=True)
        parser.add_argument('-f', '--gff3', required=True)
        parser.add_argument('-o', '--output', required=True)
        parser.add_argument('--gene-info')
        try:
            args = parser.parse_args()
            print(f"Argparse would have parsed arguments as: {args}")
        except Exception as e:
            print(f"Argparse simulation failed: {e}")
        
        print("--- ✅ Extraction finished (simulated) ---")
    return main_placeholder
longest_mrna_main = get_original_main_for_demo()
# END: Placeholder

@click.command(context_settings=dict(help_option_names=['-h', '--help']), short_help = "从GFF3注释中为每个基因提取最长的mRNA转录本")
# --- Required arguments ---
@click.option('--genome', '-g',
              required=True,
              type=click.Path(exists=True, dir_okay=False, resolve_path=True),
              help='🧬 输入基因组FASTA文件 | Input genome FASTA file.')
@click.option('--gff3', '-f',
              required=True,
              type=click.Path(exists=True, dir_okay=False, resolve_path=True),
              help='📜 输入GFF3注释文件 | Input GFF3 annotation file.')
@click.option('--output', '-o',
              required=True,
              type=click.Path(dir_okay=False, resolve_path=True),
              help='💾 输出蛋白质FASTA文件 | Output protein FASTA file.')
# --- Optional arguments ---
@click.option('--gene-info',
              type=click.Path(dir_okay=False, resolve_path=True),
              help='ℹ️ 基因信息输出文件 (默认自动生成) | Gene info output file (auto-generated).')
def longestmrna(genome, gff3, output, gene_info):
    """
    从GFF3注释中为每个基因提取最长的mRNA转录本.

    该工具解析GFF3文件，确定每个基因的最长编码序列(CDS)，
    然后从基因组FASTA文件中提取并翻译成对应的蛋白质序列。
    
    🌟 示例 | Examples:
    
    \b
    # 🎯 基本用法
    biopytools longest-mrna -g genome.fa -f annotation.gff3 -o longest_proteins.fa
    
    \b
    # ℹ️ 自定义基因信息输出文件
    biopytools longest-mrna -g genome.fa -f annotation.gff3 \\
        -o longest_proteins.fa --gene-info my_gene_info.tsv
    """
    
    # 构建参数列表以传递给原始的main函数 🔄 | Build argument list for original main function
    args = ['biopytools', 'longest-mrna']
    
    # 必需参数 📋 | Required parameters
    args.extend(['-g', genome])
    args.extend(['-f', gff3])
    args.extend(['-o', output])
    
    # 可选参数 ⚙️ | Optional parameters
    if gene_info:
        args.extend(['--gene-info', gene_info])
        
    # 保存并恢复sys.argv 💾 | Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args
    
    try:
        # 调用原始的main函数 🚀 | Call original main function
        longest_mrna_main()
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
    longestmrna()