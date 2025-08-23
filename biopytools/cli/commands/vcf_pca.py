"""
🧬 VCF PCA分析命令 | VCF PCA Analysis Command
"""

import click
import sys
# In your actual project, you would use a relative import like this:
# from ...vcf_pca.main import main as vcf_pca_main

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
        parser.add_argument('-o', '--output', default='./pca_output')
        parser.add_argument('-s', '--sample-info')
        parser.add_argument('-c', '--components', type=int, default=10)
        parser.add_argument('-m', '--maf', type=float, default=0.05)
        parser.add_argument('--missing', type=float, default=0.1)
        parser.add_argument('--hwe', type=float, default=1e-6)
        parser.add_argument('--skip-qc', action='store_true')
        parser.add_argument('-p', '--plot', action='store_true')
        parser.add_argument('-g', '--group-column')
        parser.add_argument('--plink-path', default='plink')
        parser.add_argument('--bcftools-path', default='bcftools')
        try:
            args = parser.parse_args()
            print(f"Argparse would have parsed arguments as: {args}")
        except Exception as e:
            print(f"Argparse simulation failed: {e}")
        
        print("--- ✅ Analysis finished (simulated) ---")
    return main_placeholder
vcf_pca_main = get_original_main_for_demo()
# END: Placeholder

@click.command(context_settings=dict(help_option_names=['-h', '--help']))
# --- Required arguments ---
@click.option('--vcf', '-v',
              required=True,
              type=click.Path(exists=True, dir_okay=False, resolve_path=True),
              help='🧬 输入VCF文件路径 | Input VCF file path.')
# --- Optional arguments ---
@click.option('--output', '-o', 'output_dir',
              default='./pca_output', show_default=True,
              type=click.Path(file_okay=False, resolve_path=True),
              help='📂 输出目录 | Output directory.')
@click.option('--sample-info', '-s', 'sample_info_file',
              type=click.Path(exists=True, dir_okay=False, resolve_path=True),
              help='👥 样本信息文件 (用于分组绘图) | Sample information file (for plot grouping).')
# --- PCA parameters ---
@click.option('--components', '-c',
              type=int, default=10, show_default=True,
              help='📊 主成分数量 | Number of principal components.')
# --- Quality control parameters ---
@click.option('--maf', '-m',
              type=float, default=0.05, show_default=True,
              help='🗑️ 最小等位基因频率阈值 | MAF threshold.')
@click.option('--missing',
              type=float, default=0.1, show_default=True,
              help='🗑️ 最大缺失率阈值 | Maximum missing rate threshold.')
@click.option('--hwe',
              type=float, default=1e-6, show_default=True,
              help='🗑️ HWE p值阈值 | HWE p-value threshold.')
@click.option('--skip-qc',
              is_flag=True,
              help='⏩ 跳过质量控制过滤 | Skip quality control filtering.')
# --- Visualization parameters ---
@click.option('--plot', '-p',
              is_flag=True,
              help='📈 生成PCA可视化图表 | Generate PCA visualization plots.')
@click.option('--group-column', '-g',
              help='🏷️ 样本信息文件中用于分组的列名 | Column name for grouping in sample info file.')
# --- Tool paths ---
@click.option('--plink-path',
              default='plink', show_default=True,
              help='⚙️ PLINK软件路径 | PLINK software path.')
@click.option('--bcftools-path',
              default='bcftools', show_default=True,
              help='⚙️ BCFtools软件路径 | BCFtools software path.')
def vcf_pca(vcf, output_dir, sample_info_file, components, maf, missing, hwe, skip_qc, plot, group_column, plink_path, bcftools_path):
    """
    VCF PCA分析工具.

    使用PLINK对VCF文件进行主成分分析(PCA)，并可选地生成可视化图表，
    用于探究群体结构。
    
    🌟 示例 | Examples:
    
    \b
    # 🎯 基本PCA分析
    biopytools vcf-pca -v variants.vcf.gz -o pca_results
    
    \b
    # 📈 生成图表并按群体分组
    biopytools vcf-pca -v data.vcf -o results -p -s samples.tsv -g Population
        
    \b
    # ⏩ 跳过QC，直接对已有VCF进行PCA
    biopytools vcf-pca -v clean_data.vcf -o out --skip-qc -p
    """
    
    # 构建参数列表以传递给原始的main函数 🔄 | Build argument list for original main function
    args = ['biopytools', 'vcf-pca']
    
    # 必需参数 📋 | Required parameters
    args.extend(['-v', vcf])
    
    # 可选参数 ⚙️ | Optional parameters
    if output_dir != './pca_output':
        args.extend(['-o', output_dir])
    else: # 总是传递输出目录
        args.extend(['-o', output_dir])
        
    if sample_info_file:
        args.extend(['-s', sample_info_file])
        
    if components != 10:
        args.extend(['-c', str(components)])
        
    if maf != 0.05:
        args.extend(['-m', str(maf)])
        
    if missing != 0.1:
        args.extend(['--missing', str(missing)])
        
    if hwe != 1e-6:
        args.extend(['--hwe', str(hwe)])
        
    if skip_qc:
        args.append('--skip-qc')
        
    if plot:
        args.append('-p')
        
    if group_column:
        args.extend(['-g', group_column])
        
    if plink_path != 'plink':
        args.extend(['--plink-path', plink_path])
        
    if bcftools_path != 'bcftools':
        args.extend(['--bcftools-path', bcftools_path])
        
    # 保存并恢复sys.argv 💾 | Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args
    
    try:
        # 调用原始的main函数 🚀 | Call original main function
        vcf_pca_main()
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
    vcf_pca()