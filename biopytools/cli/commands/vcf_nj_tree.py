"""
🧬 VCF系统发育分析命令 | VCF Phylogenetic Analysis Command
"""

import click
import sys
# In your actual project, you would use a relative import like this:
# from ...vcf_phylo.main import main as vcf_phylo_main

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
        parser.add_argument("-i", "--input", dest="vcf_file")
        parser.add_argument("-d", "--distance-matrix")
        parser.add_argument("-o", "--output", dest="output_prefix", default="phylo_analysis")
        parser.add_argument("-t", "--tree-output")
        parser.add_argument("--vcf2dis-path", default="VCF2Dis")
        parser.add_argument("-w", "--working-dir", default=".")
        parser.add_argument("--skip-vcf2dis", action="store_true")
        try:
            args = parser.parse_args()
            # Simulate validation logic
            if args.skip_vcf2dis and not args.distance_matrix:
                print("Error: --distance-matrix must be specified when using --skip-vcf2dis")
                sys.exit(1)
            if not args.skip_vcf2dis and not args.vcf_file:
                print("Error: Input VCF file must be specified")
                sys.exit(1)
            print(f"Argparse would have parsed arguments as: {args}")
        except Exception as e:
            print(f"Argparse simulation failed: {e}")
        
        print("--- ✅ Analysis finished (simulated) ---")
    return main_placeholder
vcf_phylo_main = get_original_main_for_demo()
# END: Placeholder

@click.command(context_settings=dict(help_option_names=['-h', '--help']))
# --- Input files ---
@click.option('--vcf-file', '-i',
              type=click.Path(exists=True, dir_okay=False, resolve_path=True),
              help='🧬 输入VCF文件路径 | Input VCF file path.')
@click.option('--distance-matrix', '-d',
              type=click.Path(exists=True, dir_okay=False, resolve_path=True),
              help='📊 已有的距离矩阵文件 (用于跳过VCF2Dis) | Existing distance matrix file.')
# --- Output files ---
@click.option('--output-prefix', '-o',
              default='phylo_analysis', show_default=True,
              help='📝 输出文件前缀 | Output file prefix.')
@click.option('--tree-output', '-t',
              type=click.Path(dir_okay=False, resolve_path=True),
              help='🌳 系统发育树输出文件路径 (默认: <prefix>.nwk) | Phylogenetic tree output file path.')
# --- Tool settings ---
@click.option('--vcf2dis-path',
              default='VCF2Dis', show_default=True,
              help='⚙️ VCF2Dis程序路径 | VCF2Dis program path.')
@click.option('--working-dir', '-w',
              default='.', show_default=True,
              type=click.Path(file_okay=False, resolve_path=True),
              help='📂 工作目录 | Working directory.')
# --- Behavior control ---
@click.option('--skip-vcf2dis',
              is_flag=True,
              help='⏩ 跳过VCF2Dis步骤，直接从距离矩阵构建树 | Skip VCF2Dis step.')
def vcf_nj_tree(vcf_file, distance_matrix, output_prefix, tree_output, vcf2dis_path, working_dir, skip_vcf2dis):
    """
    VCF文件构建NJ系统发育树.

    使用VCF2Dis和邻接法(NJ)从VCF文件构建系统发育树。
    
    🌟 示例 | Examples:
    
    \b
    # 🎯 基本用法 (VCF -> 树)
    biopytools vcf-phylo -i variants.vcf -o my_tree
    
    \b
    # ⏩ 从已有的距离矩阵构建树
    biopytools vcf-phylo --distance-matrix dist.mat -o my_tree --skip-vcf2dis
        
    \b
    # 📂 指定所有路径
    biopytools vcf-phylo -i data.vcf -o results/tree_prefix -t results/final_tree.nwk
    """
    
    # 验证参数逻辑 | Validate argument logic
    if skip_vcf2dis and not distance_matrix:
        raise click.UsageError("❌ 使用 --skip-vcf2dis 时必须提供 --distance-matrix | --distance-matrix must be provided with --skip-vcf2dis.")
    if not skip_vcf2dis and not vcf_file:
        raise click.UsageError("❌ 必须提供 --vcf-file (除非使用 --skip-vcf2dis) | --vcf-file must be provided (unless using --skip-vcf2dis).")

    # 构建参数列表以传递给原始的main函数 🔄 | Build argument list for original main function
    args = ['biopytools', 'vcf-phylo']
    
    # 输入参数 | Input parameters
    if vcf_file:
        args.extend(['-i', vcf_file])
    if distance_matrix:
        args.extend(['-d', distance_matrix])
        
    # 输出参数 | Output parameters
    if output_prefix != 'phylo_analysis':
        args.extend(['-o', output_prefix])
    else: # 总是传递，让argparse处理默认值
        args.extend(['-o', output_prefix])
        
    if tree_output:
        args.extend(['-t', tree_output])
        
    # 工具参数 | Tool parameters
    if vcf2dis_path != 'VCF2Dis':
        args.extend(['--vcf2dis-path', vcf2dis_path])
        
    if working_dir != '.':
        args.extend(['-w', working_dir])
        
    # 行为参数 | Behavior parameters
    if skip_vcf2dis:
        args.append('--skip-vcf2dis')
        
    # 保存并恢复sys.argv 💾 | Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args
    
    try:
        # 调用原始的main函数 🚀 | Call original main function
        vcf_phylo_main()
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
    vcf_nj_tree()