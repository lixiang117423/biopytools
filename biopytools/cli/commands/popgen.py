"""
🧬 群体遗传分析命令 | Population Genetics Analysis Command
"""

import click
import sys
# In your actual project, you would use a relative import like this:
# from ...popgen.main import main as popgen_main

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
        parser.add_argument('-o', '--output', default='./popgen_output')
        parser.add_argument('-g', '--groups')
        parser.add_argument('--all', action='store_true', default=True)
        parser.add_argument('--fst', action='store_true')
        parser.add_argument('--pi', action='store_true')
        parser.add_argument('--theta-w', action='store_true')
        parser.add_argument('--tajima-d', action='store_true')
        parser.add_argument('--ibd', action='store_true')
        parser.add_argument('--ld', action='store_true')
        parser.add_argument('--ne', action='store_true')
        parser.add_argument('-w', '--windows', nargs='+', type=int, default=[10000, 100000, 500000])
        parser.add_argument('--overlap', type=float, default=0.9)
        parser.add_argument('-m', '--maf', type=float, default=0.01)
        parser.add_argument('-M', '--missing', type=float, default=0.1)
        parser.add_argument('-H', '--hwe', type=float, default=1e-6)
        parser.add_argument('--min-dp', type=int, default=10)
        parser.add_argument('--max-dp', type=int, default=100)
        parser.add_argument('-f', '--format', choices=['txt', 'csv', 'tsv', 'json'], default='txt')
        parser.add_argument('-t', '--threads', type=int, default=4)
        parser.add_argument('--vcftools-path', default='vcftools')
        parser.add_argument('--plink-path', default='plink')
        parser.add_argument('--bcftools-path', default='bcftools')
        parser.add_argument('--smcpp-path', default='smc++')
        try:
            args = parser.parse_args()
            print(f"Argparse would have parsed arguments as: {args}")
        except Exception as e:
            print(f"Argparse simulation failed: {e}")
        
        print("--- ✅ Analysis finished (simulated) ---")
    return main_placeholder
popgen_main = get_original_main_for_demo()
# END: Placeholder

@click.command(context_settings=dict(help_option_names=['-h', '--help']), short_help = "群体遗传多样性参数计算工具")
# --- Required arguments ---
@click.option('--vcf', '-v',
              required=True,
              type=click.Path(exists=True, dir_okay=False, resolve_path=True),
              help='🧬 输入VCF文件路径 | Input VCF file path.')
# --- Optional arguments ---
@click.option('--output', '-o',
              default='./popgen_output', show_default=True,
              type=click.Path(file_okay=False, resolve_path=True),
              help='📂 输出目录 | Output directory.')
@click.option('--groups', '-g',
              type=click.Path(exists=True, dir_okay=False, resolve_path=True),
              help='👥 分组信息文件 | Group information file.')
# --- Analysis selection parameters ---
@click.option('--all', 'run_all', is_flag=True, default=True, show_default=True, help='🔬 计算所有参数 (默认开启) | Calculate all parameters.')
@click.option('--fst', is_flag=True, help='📈 计算Fst | Calculate Fst.')
@click.option('--pi', is_flag=True, help='📈 计算π | Calculate π.')
@click.option('--theta-w', is_flag=True, help='📈 计算θw | Calculate θw.')
@click.option('--tajima-d', is_flag=True, help='📈 计算Tajima\'s D | Calculate Tajima\'s D.')
@click.option('--ibd', is_flag=True, help='🔗 计算IBD | Calculate IBD.')
@click.option('--ld', is_flag=True, help='🔗 计算LD | Calculate LD.')
@click.option('--ne', is_flag=True, help='📈 计算有效群体大小 | Calculate effective population size.')
# --- Sliding window parameters ---
@click.option('--windows', '-w',
              multiple=True, type=int, default=[10000, 100000, 500000], show_default=True,
              help='🪟 滑动窗口大小(bp) | Sliding window sizes (bp).')
@click.option('--overlap',
              type=float, default=0.9, show_default=True,
              help='🔄 窗口重叠率 | Window overlap rate.')
# --- Quality control parameters ---
@click.option('--maf', '-m',
              type=float, default=0.01, show_default=True,
              help='🗑️ MAF阈值 | MAF threshold.')
@click.option('--missing', '-M',
              type=float, default=0.1, show_default=True,
              help='🗑️ 缺失率阈值 | Missing rate threshold.')
@click.option('--hwe', '-H',
              type=float, default=1e-6, show_default=True,
              help='🗑️ HWE p值阈值 | HWE p-value threshold.')
@click.option('--min-dp',
              type=int, default=10, show_default=True,
              help='✅ 最小测序深度 | Minimum depth.')
@click.option('--max-dp',
              type=int, default=100, show_default=True,
              help='🚫 最大测序深度 | Maximum depth.')
# --- Output and resources ---
@click.option('--format', '-f',
              type=click.Choice(['txt', 'csv', 'tsv', 'json'], case_sensitive=False),
              default='txt', show_default=True,
              help='📄 输出文件格式 | Output file format.')
@click.option('--threads', '-t',
              type=int, default=4, show_default=True,
              help='🚀 线程数 | Number of threads.')
# --- Tool paths ---
@click.option('--vcftools-path', default='vcftools', show_default=True, help='🔧 VCFtools路径 | VCFtools path.')
@click.option('--plink-path', default='plink', show_default=True, help='🔧 PLINK路径 | PLINK path.')
@click.option('--bcftools-path', default='bcftools', show_default=True, help='🔧 BCFtools路径 | BCFtools path.')
@click.option('--smcpp-path', default='smc++', show_default=True, help='🔧 SMC++路径 | SMC++ path.')
def popgen(**kwargs):
    """
    群体遗传多样性参数计算工具.

    一个集成的工具，用于执行一系列标准的群体遗传学分析，
    包括遗传多样性、群体分化、连锁不平衡等。
    
    🌟 示例 | Examples:
    
    \b
    # 🎯 运行所有分析 (默认行为)
    biopytools popgen -v variants.vcf.gz -g groups.txt -o results
    
    \b
    # 🔬 只运行Fst和多样性分析
    biopytools popgen -v data.vcf -g groups.txt -o out --no-all --fst --pi
    
    \b
    # 🚀 使用更多线程并调整QC
    biopytools popgen -v data.vcf -o out --maf 0.05 -t 16
    """
    
    # 构建参数列表以传递给原始的main函数 🔄 | Build argument list for original main function
    args = ['biopytools', 'popgen']
    
    # 特殊处理 --all 标志的逻辑
    # click v8+ has `allow_from_autoenv`, but for compatibility, we handle it manually.
    # The default for the 'run_all' flag is True. click passes it as True unless --no-all is used.
    # The original argparse script has a bug where --all is always True.
    # We will correctly represent the user's intent.
    specific_analyses_chosen = any([
        kwargs.get('fst'), kwargs.get('pi'), kwargs.get('theta_w'),
        kwargs.get('tajima_d'), kwargs.get('ibd'), kwargs.get('ld'), kwargs.get('ne')
    ])
    
    # If specific analyses are chosen, we must assume the user wants to override 'all'.
    # The original argparse script doesn't support --no-all, so we achieve this by NOT passing --all.
    # But since its default is True, this won't work without fixing the argparse script.
    # The following code correctly translates the click command, but relies on a fixed argparse script.
    
    if kwargs.get('run_all') and not specific_analyses_chosen:
        args.append('--all')
        
    # 遍历所有参数
    for key, value in kwargs.items():
        if value is None or key == 'run_all': # Skip 'run_all' as we handled it
            continue

        param_name = '--' + key.replace('_', '-')
        
        # 处理布尔标志
        if isinstance(value, bool) and value:
            args.append(param_name)
        # 处理多值参数 (windows)
        elif isinstance(value, tuple) and value:
            default_val = popgen.params_by_name[key].default
            if sorted(value) != sorted(default_val):
                args.append(param_name)
                args.extend(map(str, value))
        # 处理非布尔、非默认值的常规参数
        elif not isinstance(value, (bool, tuple)):
            default_val = popgen.params_by_name[key].default
            if value != default_val:
                args.append(param_name)
                args.append(str(value))
    
    # 补上必需参数
    args.extend(['-v', kwargs['vcf']])

    # 保存并恢复sys.argv 💾 | Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args
    
    try:
        # 调用原始的main函数 🚀 | Call original main function
        popgen_main()
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
    popgen()