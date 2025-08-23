"""
🧬 VCF基因型提取命令 | VCF Genotype Extraction Command
"""

import click
import sys
# In your actual project, you would use a relative import like this:
# from ...vcf_genotype.main import main as vcf_genotype_main

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
        parser.add_argument('-i', '--input', dest='vcf_file', required=True)
        parser.add_argument('-o', '--output', default='vcf_genotype')
        parser.add_argument('-s', '--samples', default='all')
        parser.add_argument('--biallelic-only', action='store_true')
        parser.add_argument('-e', '--each', choices=['yes', 'y', 'no', 'n'], default='n')
        parser.add_argument('-t', '--output-type', choices=['txt', 'csv', 'excel'], default='txt')
        parser.add_argument('--output-dir', default='./')
        try:
            args = parser.parse_args()
            print(f"Argparse would have parsed arguments as: {args}")
        except Exception as e:
            print(f"Argparse simulation failed: {e}")
        
        print("--- ✅ Extraction finished (simulated) ---")
    return main_placeholder
vcf_genotype_main = get_original_main_for_demo()
# END: Placeholder

@click.command(context_settings=dict(help_option_names=['-h', '--help']), short_help = "VCF基因型提取工具")
# --- Required arguments ---
@click.option('--input', '-i', 'vcf_file',
              required=True,
              type=click.Path(exists=True, dir_okay=False, resolve_path=True),
              help='📁 VCF文件路径 (支持.gz) | VCF file path (supports .gz).')
# --- Optional arguments ---
@click.option('--output', '-o', 'output_prefix',
              default='vcf_genotype', show_default=True,
              help='📝 输出文件前缀 | Output file prefix.')
@click.option('--samples', '-s',
              default='all', show_default=True,
              help='👥 样本选择 (all 或 逗号分隔) | Sample selection (all or comma-separated).')
@click.option('--biallelic-only',
              is_flag=True,
              help='🔬 只保留双等位位点 | Keep only biallelic sites.')
@click.option('--each', '-e', 'split_by_chromosome',
              type=click.Choice(['yes', 'y', 'no', 'n'], case_sensitive=False),
              default='n', show_default=True,
              help='🧬 按染色体拆分输出文件 | Split output files by chromosome.')
@click.option('--output-type', '-t',
              type=click.Choice(['txt', 'csv', 'excel'], case_sensitive=False),
              default='txt', show_default=True,
              help='📊 输出文件格式 | Output file format.')
@click.option('--output-dir',
              default='./', show_default=True,
              type=click.Path(file_okay=False, resolve_path=True),
              help='📂 输出目录 | Output directory.')
def vcf_genotype(vcf_file, output_prefix, samples, biallelic_only, split_by_chromosome, output_type, output_dir):
    """
    VCF基因型提取工具.

    从VCF文件中提取指定样本的基因型信息，并支持多种格式输出。
    
    🌟 示例 | Examples:
    
    \b
    # 🎯 提取所有样本的基因型
    biopytools vcf-genotype -i variants.vcf.gz -o all_samples
    
    \b
    # 👥 提取指定样本，并按染色体拆分输出
    biopytools vcf-genotype -i data.vcf -s "sample1,sample2" -e yes
        
    \b
    # 🔬 只提取双等位位点，并输出为Excel格式
    biopytools vcf-genotype -i data.vcf -o biallelic_results -t excel --biallelic-only
    """
    
    # 构建参数列表以传递给原始的main函数 🔄 | Build argument list for original main function
    args = ['biopytools', 'vcf-genotype']
    
    # 必需参数 📋 | Required parameters
    args.extend(['-i', vcf_file])
    
    # 可选参数 ⚙️ | Optional parameters
    if output_prefix != 'vcf_genotype':
        args.extend(['-o', output_prefix])
        
    if samples != 'all':
        args.extend(['-s', samples])
        
    if biallelic_only:
        args.append('--biallelic-only')
        
    if split_by_chromosome.lower() not in ['n', 'no']:
        args.extend(['-e', split_by_chromosome])
        
    if output_type != 'txt':
        args.extend(['-t', output_type])
        
    # The default './' is tricky. If user provides it explicitly, we pass it.
    # Otherwise, argparse handles the default. For simplicity, we pass it if it's not the default.
    if output_dir != './':
        args.extend(['--output-dir', output_dir])

    # 保存并恢复sys.argv 💾 | Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args
    
    try:
        # 调用原始的main函数 🚀 | Call original main function
        vcf_genotype_main()
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
    vcf_genotype()