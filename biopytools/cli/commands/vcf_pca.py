"""
VCF PCA分析命令|VCF PCA Analysis Command
"""

import click
import sys


def _lazy_import_vcf_pca_main():
    """延迟加载vcf_pca主函数|Lazy load vcf_pca main function"""
    try:
        from ...vcf_pca.main import main as vcf_pca_main
        return vcf_pca_main
    except ImportError as e:
        click.echo(f"导入错误|Import error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否是帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


@click.command(short_help='VCF PCA主成分分析|VCF PCA principal component analysis',
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--input', '-i',
              required=True,
              type=click.Path(exists=True),
              help='输入VCF文件|Input VCF file')
@click.option('--output', '-o',
              default='./pca_output',
              type=click.Path(),
              show_default=True,
              help='输出目录|Output directory')
@click.option('--sample-info', '-s',
              type=click.Path(exists=True),
              help='样本信息文件|Sample information file')
@click.option('--components', '-c',
              default=10,
              type=int,
              show_default=True,
              help='主成分数量|Number of principal components')
@click.option('--maf',
              default=0.05,
              type=float,
              show_default=True,
              help='最小等位基因频率阈值|Minor allele frequency threshold')
@click.option('--missing',
              default=0.1,
              type=float,
              show_default=True,
              help='最大缺失率阈值|Maximum missing rate threshold')
@click.option('--hwe',
              default=1e-6,
              type=float,
              show_default=True,
              help='Hardy-Weinberg平衡p值阈值|Hardy-Weinberg equilibrium p-value threshold')
@click.option('--skip-qc',
              is_flag=True,
              help='跳过质量控制过滤|Skip quality control filtering')
@click.option('--plot', '-P',
              is_flag=True,
              help='生成PCA可视化图表|Generate PCA visualization plots')
@click.option('--group-column', '-g',
              type=str,
              help='分组列名|Column name for grouping')
@click.option('--plink-path',
              default='plink',
              type=str,
              show_default=True,
              help='PLINK软件路径|PLINK software path')
@click.option('--bcftools-path',
              default='bcftools',
              type=str,
              show_default=True,
              help='BCFtools软件路径|BCFtools software path')
def vcf_pca(input, output, sample_info, components, maf, missing, hwe, skip_qc,
            plot, group_column, plink_path, bcftools_path):
    """
    VCF PCA主成分分析工具|VCF PCA Principal Component Analysis Tool

    示例|Examples: biopytools vcf-pca -i variants.vcf -o pca_results
    """

    # 构建参数列表传递给原始main函数|Build argument list for original main function
    args = ['vcf_pca.py']

    # 必需参数|Required parameters
    args.extend(['--vcf', input])

    # 可选参数（只在非默认值时添加）| Optional parameters (add only when non-default)
    if output != './pca_output':
        args.extend(['-o', output])

    if sample_info:
        args.extend(['-s', sample_info])

    if components != 10:
        args.extend(['-c', str(components)])

    if maf != 0.05:
        args.extend(['--maf', str(maf)])

    if missing != 0.1:
        args.extend(['--missing', str(missing)])

    if hwe != 1e-6:
        args.extend(['--hwe', str(hwe)])

    if group_column:
        args.extend(['-g', group_column])

    if plink_path != 'plink':
        args.extend(['--plink-path', plink_path])

    if bcftools_path != 'bcftools':
        args.extend(['--bcftools-path', bcftools_path])

    # 标志参数|Flag parameters
    if skip_qc:
        args.append('--skip-qc')

    if plot:
        args.append('--plot')

    # 保存并恢复sys.argv|Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        # 延迟加载|Lazy loading
        vcf_pca_main_func = _lazy_import_vcf_pca_main()

        # 调用原始的main函数|Call original main function
        vcf_pca_main_func()
    except SystemExit as e:
        # 处理程序正常退出|Handle normal program exit
        if e.code != 0:
            sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n用户中断|User interrupted", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"运行时错误|Runtime error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
