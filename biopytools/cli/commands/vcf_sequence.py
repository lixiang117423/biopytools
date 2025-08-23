"""
🧬 序列提取命令 | Sequence Extraction Command
"""

import click
import sys
# In your actual project, you would use a relative import like this:
# from ...sequence_extractor.main import main as sequence_extractor_main

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
        parser.add_argument('-v', '--vcf', required=True)
        parser.add_argument('-g', '--genome', required=True)
        parser.add_argument('-c', '--chrom', required=True)
        parser.add_argument('-s', '--start', type=int, required=True)
        parser.add_argument('-e', '--end', type=int, required=True)
        parser.add_argument('-o', '--output-dir', default='./sequence_output')
        parser.add_argument('--format', choices=['tab', 'fasta', 'csv'], default='tab')
        parser.add_argument('--second-allele', action='store_true')
        parser.add_argument('--no-reference', action='store_true')
        parser.add_argument('--min-qual', type=int)
        parser.add_argument('--samples')
        parser.add_argument('--exclude-samples')
        try:
            args = parser.parse_args()
            # The original script does more processing after this, but we simulate the parsing.
            print(f"Argparse would have parsed arguments as: {args}")
        except Exception as e:
            print(f"Argparse simulation failed: {e}")
        
        print("--- ✅ Extraction finished (simulated) ---")
        return 0 # The original main returns an exit code
    return main_placeholder
sequence_extractor_main = get_original_main_for_demo()
# END: Placeholder

@click.command(context_settings=dict(help_option_names=['-h', '--help']))
# --- Required arguments ---
@click.option('--vcf', '-v',
              required=True,
              type=click.Path(exists=True, dir_okay=False, resolve_path=True),
              help='🧬 VCF文件路径 | VCF file path.')
@click.option('--genome', '-g',
              required=True,
              type=click.Path(exists=True, dir_okay=False, resolve_path=True),
              help='🧬 基因组FASTA文件路径 | Genome FASTA file path.')
@click.option('--chrom', '-c',
              required=True,
              help='🗺️ 染色体名称 | Chromosome name.')
@click.option('--start', '-s',
              type=int, required=True,
              help='➡️ 起始位置 (1-based) | Start position (1-based).')
@click.option('--end', '-e',
              type=int, required=True,
              help='⬅️ 结束位置 (1-based, inclusive) | End position (1-based, inclusive).')
# --- Optional arguments ---
@click.option('--output-dir', '-o',
              default='./sequence_output', show_default=True,
              type=click.Path(file_okay=False, resolve_path=True),
              help='📂 输出目录 | Output directory.')
@click.option('--format',
              type=click.Choice(['tab', 'fasta', 'csv'], case_sensitive=False),
              default='tab', show_default=True,
              help='📄 输出格式 | Output format.')
@click.option('--second-allele',
              is_flag=True,
              help='🥈 使用第二个等位基因 (默认使用第一个) | Use second allele instead of first.')
@click.option('--no-reference',
              is_flag=True,
              help='🚫 不在输出中包含参考序列 | Do not include reference sequence in output.')
@click.option('--min-qual',
              type=int,
              help='✅ 最小质量值过滤 | Minimum quality filter.')
@click.option('--samples',
              help='👥 指定样本 (文件路径或逗号分隔) | Specify samples (file path or comma-separated).')
@click.option('--exclude-samples',
              help='🗑️ 排除样本 (文件路径或逗号分隔) | Exclude samples (file path or comma-separated).')
def vcf_sequence(**kwargs):
    """
    VCF和基因组文件中提取特定区间的序列.

    根据VCF文件中的变异信息，为每个样本重建指定基因组区域的序列。
    
    🌟 示例 | Examples:
    
    \b
    # 🎯 基本用法
    biopytools sequence-extractor -v d.vcf -g g.fa -c chr1 -s 1000 -e 2000
    
    \b
    # 📄 输出为FASTA格式，并指定样本
    biopytools sequence-extractor -v d.vcf -g g.fa -c chr2 -s 500 -e 1500 \\
        --format fasta --samples "sample1,sample2"
        
    \b
    # 🥈 使用第二个等位基因，并进行质量过滤
    biopytools sequence-extractor -v d.vcf -g g.fa -c chr3 -s 1 -e 1000 \\
        --second-allele --min-qual 30
    """
    
    # 构建参数列表以传递给原始的main函数 🔄 | Build argument list for original main function
    args = ['biopytools', 'sequence-extractor']
    
    # 遍历所有参数
    for key, value in kwargs.items():
        if value is None or value is False: # Skip None and False flags
            continue
            
        param_name = '--' + key.replace('_', '-')
        
        # 处理布尔标志
        if isinstance(value, bool) and value:
            args.append(param_name)
        # 处理非布尔参数
        elif not isinstance(value, bool):
            args.append(param_name)
            args.append(str(value))

    # 保存并恢复sys.argv 💾 | Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args
    
    try:
        # 调用原始的main函数 🚀 | Call original main function
        exit_code = sequence_extractor_main()
        if exit_code != 0:
            sys.exit(exit_code)
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
    vcf_sequence()