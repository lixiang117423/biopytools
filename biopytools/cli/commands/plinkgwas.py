"""
🧬 PLINK GWAS分析命令 | PLINK GWAS Analysis Command
"""

import click
import sys
# In your actual project, you would use a relative import like this:
# from ...plink_gwas.main import main as plink_gwas_main

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
        parser.add_argument("-v", "--vcf-file", required=True)
        parser.add_argument("-p", "--phenotype-file", required=True)
        parser.add_argument("-t", "--trait-type", choices=["qualitative", "quantitative"], default="qualitative")
        parser.add_argument('-m', '--genetic-model', choices=['additive', 'dominant', 'recessive', 'all'], default='additive')
        parser.add_argument("-o", "--output-dir", default="gwas_results")
        parser.add_argument("--mind", type=float, default=0.05)
        parser.add_argument("--geno", type=float, default=0.05)
        parser.add_argument("--maf", type=float, default=0.01)
        parser.add_argument("--hwe", type=float, default=1e-6)
        parser.add_argument("--ld-window-size", type=int, default=50)
        parser.add_argument("--ld-step-size", type=int, default=5)
        parser.add_argument("--ld-r2-threshold", type=float, default=0.2)
        parser.add_argument("--pca-components", type=int, default=10)
        parser.add_argument("--pca-use", type=int, default=5)
        parser.add_argument("--correction-method", choices=["bonferroni", "suggestive", "fdr", "all"], default="all")
        parser.add_argument("--bonferroni-alpha", type=float, default=0.05)
        parser.add_argument("--suggestive-threshold", type=float, default=1e-5)
        parser.add_argument("--fdr-alpha", type=float, default=0.05)
        parser.add_argument("--threads", type=int, default=1)
        try:
            args = parser.parse_args()
            print(f"Argparse would have parsed arguments as: {args}")
        except Exception as e:
            print(f"Argparse simulation failed: {e}")
        
        print("--- ✅ Analysis finished (simulated) ---")
    return main_placeholder
plink_gwas_main = get_original_main_for_demo()
# END: Placeholder

@click.group(invoke_without_command=True, context_settings=dict(help_option_names=['-h', '--help']), short_help = "完整的PLINK GWAS分析流程")
@click.option('--vcf-file', '-v',
              required=True,
              type=click.Path(exists=True, dir_okay=False, resolve_path=True),
              help='🧬 VCF文件路径 (支持.gz) | Input VCF file path (supports .gz).')
@click.option('--phenotype-file', '-p',
              required=True,
              type=click.Path(exists=True, dir_okay=False, resolve_path=True),
              help='📋 表型文件路径 | Phenotype file path.')
@click.option('--trait-type', '-t',
              type=click.Choice(['qualitative', 'quantitative'], case_sensitive=False),
              default='qualitative', show_default=True,
              help='📊 表型类型 | Trait type.')
@click.option('--genetic-model', '-m',
              type=click.Choice(['additive', 'dominant', 'recessive', 'all'], case_sensitive=False),
              default='additive', show_default=True,
              help='🔬 遗传模型 | Genetic model.')
@click.option('--output-dir', '-o',
              default='gwas_results', show_default=True,
              type=click.Path(file_okay=False, resolve_path=True),
              help='📂 输出目录 | Output directory.')
@click.option('--mind',
              type=float, default=0.05, show_default=True,
              help='🗑️ 个体缺失率阈值 | Individual missing rate threshold.')
@click.option('--geno',
              type=float, default=0.05, show_default=True,
              help='🗑️ SNP缺失率阈值 | SNP missing rate threshold.')
@click.option('--maf',
              type=float, default=0.01, show_default=True,
              help='🧬 最小等位基因频率 | Minor allele frequency threshold.')
@click.option('--hwe',
              type=float, default=1e-6, show_default=True,
              help='⚖️ Hardy-Weinberg平衡P值阈值 | HWE p-value threshold.')
@click.option('--ld-window-size',
              type=int, default=50, show_default=True,
              help='🔗 LD剪枝窗口大小(kb) | LD pruning window size (kb).')
@click.option('--ld-step-size',
              type=int, default=5, show_default=True,
              help='🔗 LD剪枝步长(SNP数) | LD pruning step size (SNPs).')
@click.option('--ld-r2-threshold',
              type=float, default=0.2, show_default=True,
              help='🔗 LD剪枝r²阈值 | LD pruning r² threshold.')
@click.option('--pca-components',
              type=int, default=10, show_default=True,
              help='📊 计算的主成分数量 | Number of PCs to compute.')
@click.option('--pca-use',
              type=int, default=5, show_default=True,
              help='📊 关联分析中使用的主成分数量 | Number of PCs to use in analysis.')
@click.option('--correction-method',
              type=click.Choice(['bonferroni', 'suggestive', 'fdr', 'all'], case_sensitive=False),
              default='all', show_default=True,
              help='📈 显著性校正方法 | Significance correction method.')
@click.option('--bonferroni-alpha',
              type=float, default=0.05, show_default=True,
              help='📈 Bonferroni校正alpha水平 | Alpha for Bonferroni correction.')
@click.option('--suggestive-threshold',
              type=float, default=1e-5, show_default=True,
              help='📈 提示性关联阈值 | Suggestive association threshold.')
@click.option('--fdr-alpha',
              type=float, default=0.05, show_default=True,
              help='📈 FDR校正q值阈值 | Q-value threshold for FDR correction.')
@click.option('--threads',
              type=int, default=1, show_default=True,
              help='🚀 使用的线程数 | Number of threads to use.')
@click.pass_context
def plinkgwas(ctx, **kwargs):
    """
    完整的PLINK GWAS分析流程.

    支持质量性状和数量性状，多种遗传模型和显著性校正方法。
    
    🌟 示例 | Examples:
    
    \b
    # 🎯 基本质量性状分析 (加性模型)
    biopytools plink-gwas -v data.vcf.gz -p pheno.txt -t qualitative -o results
    
    \b
    # 🔬 使用显性模型
    biopytools plink-gwas -v data.vcf.gz -p pheno.txt -t qualitative -m dominant
    
    \b
    # 📊 数量性状分析并调整QC参数
    biopytools plink-gwas -v data.vcf.gz -p quant_pheno.txt -t quantitative \\
        --maf 0.05 --hwe 1e-5
    """
    if ctx.invoked_subcommand is None:
        # 构建参数列表以传递给原始的main函数 🔄 | Build argument list for original main function
        args = ['biopytools', 'plink-gwas']
        
        # 遍历所有click传递的参数
        for key, value in kwargs.items():
            if value is None:
                continue
                
            param_name = '--' + key.replace('_', '-')
            default_val = ctx.command.params_by_name[key].default
            
            # 只有当值不等于默认值时才添加 (除了必需参数)
            if value != default_val or ctx.command.params_by_name[key].required:
                if isinstance(value, bool) and value:
                    args.append(param_name)
                elif not isinstance(value, bool):
                    args.append(param_name)
                    args.append(str(value))

        # 保存并恢复sys.argv 💾 | Save and restore sys.argv
        original_argv = sys.argv
        sys.argv = args
        
        try:
            # 调用原始的main函数 🚀 | Call original main function
            plink_gwas_main()
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
    plinkgwas()