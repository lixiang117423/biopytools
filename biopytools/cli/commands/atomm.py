"""ATOMM CLI命令|ATOMM CLI Command"""

import click
import sys


def _lazy_import_atomm_main():
    """延迟加载ATOMM主函数|Lazy load ATOMM main function"""
    try:
        from ...atomm.main import main as atomm_main
        return atomm_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


@click.command(
    short_help='ATOMM双物种混合效应模型关联分析|ATOMM two-organism mixed model association',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--host-vcf', required=True,
              help='宿主VCF文件(支持.gz)|Host VCF file (.gz supported)')
@click.option('--pathogen-vcf', required=True,
              help='病原VCF文件(支持.gz)|Pathogen VCF file (.gz supported)')
@click.option('--phenotype-matrix', required=True,
              help='交叉感染表型矩阵(行=宿主,列=病原)|Cross-infection phenotype matrix (rows=hosts, cols=pathogens)')
@click.option('-o', '--output-dir', default='./output', show_default=True,
              help='输出目录|Output directory')
@click.option('--maf', default=0.05, type=float, show_default=True,
              help='MAF过滤阈值|MAF filter threshold')
@click.option('--encoding', default='auto', type=click.Choice(['auto', 'haploid', 'dosage']), show_default=True,
              help='基因型编码方式|Genotype encoding mode')
@click.option('--convert-maf', default=0.05, type=float, show_default=True,
              help='VCF转换时MAF阈值|MAF threshold for VCF conversion')
@click.option('--missing-value', default='NA', show_default=True,
              help='表型缺失值标记|Missing value marker in phenotype matrix')
@click.option('--host-snp-range', type=int, nargs=2,
              help='宿主边际检验SNP范围|Host marginal test SNP range')
@click.option('--pathogen-snp-range', type=int, nargs=2,
              help='病原边际检验SNP范围|Pathogen marginal test SNP range')
@click.option('--interaction-host-range', type=int, nargs=2,
              help='交互检验宿主SNP范围|Interaction test host SNP range')
@click.option('--interaction-pathogen-range', type=int, nargs=2,
              help='交互检验病原SNP范围|Interaction test pathogen SNP range')
@click.option('--tol', default=1e-6, type=float, show_default=True,
              help='优化容忍度|Optimizer tolerance')
@click.option('--maxiter', default=10000, type=int, show_default=True,
              help='最大迭代次数|Max iterations')
def atomm(host_vcf, pathogen_vcf, phenotype_matrix, output_dir,
           maf, encoding, convert_maf, missing_value,
           host_snp_range, pathogen_snp_range,
           interaction_host_range, interaction_pathogen_range,
           tol, maxiter):
    """ATOMM: 双物种混合效应模型关联分析

    基于Wang et al. 2018 PNAS的双物种混合效应模型
    Two-organism mixed-effect model based on Wang et al. 2018 PNAS
    """

    atomm_main = _lazy_import_atomm_main()

    args = ['atomm.py', '--host-vcf', host_vcf, '--pathogen-vcf', pathogen_vcf,
            '--phenotype-matrix', phenotype_matrix,
            '-o', output_dir, '--maf', str(maf),
            '--encoding', encoding, '--convert-maf', str(convert_maf),
            '--missing-value', missing_value,
            '--tol', str(tol), '--maxiter', str(maxiter)]

    if host_snp_range:
        args.extend(['--host-snp-range', str(host_snp_range[0]), str(host_snp_range[1])])
    if pathogen_snp_range:
        args.extend(['--pathogen-snp-range', str(pathogen_snp_range[0]), str(pathogen_snp_range[1])])
    if interaction_host_range:
        args.extend(['--interaction-host-range', str(interaction_host_range[0]), str(interaction_host_range[1])])
    if interaction_pathogen_range:
        args.extend(['--interaction-pathogen-range', str(interaction_pathogen_range[0]), str(interaction_pathogen_range[1])])

    original_argv = sys.argv
    sys.argv = args

    try:
        atomm_main()
    except SystemExit as e:
        if e.code != 0:
            sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n用户中断操作|Analysis interrupted by user", err=True)
        sys.exit(130)
    except Exception as e:
        click.echo(f"分析执行失败|Analysis execution failed: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
