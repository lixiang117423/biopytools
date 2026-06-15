"""
PLINK GWAS分析命令|PLINK GWAS Analysis Command
"""

import click
import sys
import os


def _lazy_import_plinkgwas_main():
    """延迟加载plink gwas主函数|Lazy load plink gwas main function"""
    try:
        from ...plinkgwas.main import main as plinkgwas_main
        return plinkgwas_main
    except ImportError as e:
        click.echo(f"导入错误|Import error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_file_exists(file_path):
    """验证文件存在(仅在非帮助模式)|Validate file existence (only in non-help mode)"""
    if _is_help_request():
        return file_path
    if not os.path.exists(file_path):
        raise click.BadParameter(f"文件不存在|File does not exist: {file_path}")
    return file_path


@click.command(
    short_help='PLINK GWAS分析工具|PLINK GWAS analysis tool',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--vcf', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='VCF文件路径|Path to VCF file')
@click.option('--phenotype', '-p',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='表型文件路径|Path to phenotype file')
@click.option('--trait-type', '-T',
              type=click.Choice(['qualitative', 'quantitative']),
              default='qualitative',
              show_default=True,
              help='表型类型|Trait type')
@click.option('--genetic-model', '-m',
              type=click.Choice(['additive', 'dominant', 'recessive', 'all']),
              default='additive',
              show_default=True,
              help='遗传模型|Genetic model')
@click.option('--output-dir', '-o',
              default='gwas_results',
              show_default=True,
              help='输出目录|Output directory')
@click.option('--no-strat-corr',
              is_flag=True,
              help='禁用群体结构校正|Disable population stratification correction')
@click.option('--mind',
              type=float,
              show_default=True,
              help='个体缺失率阈值|Individual missing rate threshold')
@click.option('--geno',
              type=float,
              default=0.05,
              show_default=True,
              help='SNP缺失率阈值|SNP missing rate threshold')
@click.option('--maf',
              type=float,
              default=0.05,
              show_default=True,
              help='最小等位基因频率|Minor allele frequency')
@click.option('--hwe',
              type=float,
              show_default=True,
              help='Hardy-Weinberg平衡P值阈值|HWE p-value threshold')
@click.option('--ld-window-size',
              type=int,
              default=50,
              show_default=True,
              help='LD窗口大小|LD window size in kb')
@click.option('--ld-step-size',
              type=int,
              default=5,
              show_default=True,
              help='LD步长大小|LD step size in SNPs')
@click.option('--ld-r2-threshold',
              type=float,
              default=0.2,
              show_default=True,
              help='LD r²阈值|LD r² threshold')
@click.option('--pca-components',
              type=int,
              default=10,
              show_default=True,
              help='主成分数量|Number of PCA components')
@click.option('--pca-use',
              type=int,
              default=5,
              show_default=True,
              help='关联分析中使用的主成分数量|Number of PCs to use in association')
@click.option('--correction-method',
              type=click.Choice(['bonferroni', 'suggestive', 'fdr', 'all']),
              default='all',
              show_default=True,
              help='显著性校正方法|Significance correction method')
@click.option('--bonferroni-alpha',
              type=float,
              default=0.05,
              show_default=True,
              help='Bonferroni校正alpha水平|Bonferroni alpha level')
@click.option('--suggestive-threshold',
              type=float,
              default=1e-5,
              show_default=True,
              help='提示性关联阈值|Suggestive threshold')
@click.option('--fdr-alpha',
              type=float,
              default=0.05,
              show_default=True,
              help='FDR校正q值阈值|FDR q-value threshold')
@click.option('--threads', '-t',
              type=int,
              default=12,
              show_default=True,
              help='线程数|Number of threads')
@click.option('--verbose', '-v',
              count=True,
              help='详细输出模式|Verbose mode (-v: INFO, -vv: DEBUG)')
@click.option('--quiet',
              is_flag=True,
              help='静默模式|Quiet mode (only ERROR)')
@click.option('--log-file',
              type=click.Path(),
              help='日志文件路径|Log file path')
@click.option('--force', '-f',
              is_flag=True,
              help='强制覆盖已存在的输出目录|Force overwrite existing output directory')
@click.option('--dry-run',
              is_flag=True,
              help='模拟运行不实际执行分析|Dry run without actual analysis')
def plinkgwas(vcf, phenotype, trait_type, genetic_model, output_dir,
              no_strat_corr, mind, geno, maf, hwe, ld_window_size, ld_step_size,
              ld_r2_threshold, pca_components, pca_use, correction_method,
              bonferroni_alpha, suggestive_threshold, fdr_alpha, threads,
              verbose, quiet, log_file, force, dry_run):
    """
    PLINK GWAS分析工具|PLINK GWAS Analysis Tool

    基于PLINK进行全基因组关联分析|Perform genome-wide association study using PLINK

    示例|Examples: biopytools plink-gwas -i data.vcf.gz -p pheno.txt -o results
    """

    # 延迟加载|Lazy load
    plinkgwas_main = _lazy_import_plinkgwas_main()

    # 构建主函数参数列表|Build argument list for main function
    args = ['plinkgwas.py']

    # 必需参数|Required parameters
    args.extend(['-i', vcf])
    args.extend(['-p', phenotype])

    # 可选参数(仅在非默认时添加)|Optional parameters (add only when non-default)
    if trait_type != 'qualitative':
        args.extend(['-T', trait_type])

    if genetic_model != 'additive':
        args.extend(['-m', genetic_model])

    if output_dir != 'gwas_results':
        args.extend(['-o', output_dir])

    if no_strat_corr:
        args.append('--no-strat-corr')

    if mind is not None:
        args.extend(['--mind', str(mind)])

    if geno != 0.05:
        args.extend(['--geno', str(geno)])

    if maf != 0.05:
        args.extend(['--maf', str(maf)])

    if hwe is not None:
        args.extend(['--hwe', str(hwe)])

    if ld_window_size != 50:
        args.extend(['--ld-window-size', str(ld_window_size)])

    if ld_step_size != 5:
        args.extend(['--ld-step-size', str(ld_step_size)])

    if ld_r2_threshold != 0.2:
        args.extend(['--ld-r2-threshold', str(ld_r2_threshold)])

    if pca_components != 10:
        args.extend(['--pca-components', str(pca_components)])

    if pca_use != 5:
        args.extend(['--pca-use', str(pca_use)])

    if correction_method != 'all':
        args.extend(['--correction-method', correction_method])

    if bonferroni_alpha != 0.05:
        args.extend(['--bonferroni-alpha', str(bonferroni_alpha)])

    if suggestive_threshold != 1e-5:
        args.extend(['--suggestive-threshold', str(suggestive_threshold)])

    if fdr_alpha != 0.05:
        args.extend(['--fdr-alpha', str(fdr_alpha)])

    if threads != 1:
        args.extend(['-t', str(threads)])

    if verbose > 0:
        args.extend(['-' + 'v' * verbose])

    if quiet:
        args.append('--quiet')

    if log_file:
        args.extend(['--log-file', log_file])

    if force:
        args.append('-f')

    if dry_run:
        args.append('--dry-run')

    # 保存并恢复sys.argv|Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        # 调用主函数|Call main function
        plinkgwas_main()
    except SystemExit as e:
        # 处理系统退出|Handle system exit
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
