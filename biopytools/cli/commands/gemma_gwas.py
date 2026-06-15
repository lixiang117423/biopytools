"""
GEMMA GWAS批量分析命令|GEMMA GWAS Batch Analysis Command
"""

import click
import sys
import os


def _lazy_import_gemma_gwas_main():
    """延迟加载gemma_gwas主函数|Lazy load gemma_gwas main function"""
    try:
        from ...gemma_gwas.main import main as gemma_gwas_main
        return gemma_gwas_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_file_exists(file_path):
    """验证文件存在(仅在非帮助模式)|Validate file exists (only in non-help mode)"""
    if not _is_help_request() and file_path and not os.path.exists(file_path):
        raise click.BadParameter(f"文件不存在|File does not exist: {file_path}")
    return file_path


@click.command(
    short_help='GEMMA GWAS批量分析|GEMMA GWAS Batch Analysis',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--input', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='VCF格式基因型文件|VCF format genotype file')
@click.option('--pheno', '-p',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='表型文件|Phenotype file (first column: sample ID, has header)')
@click.option('--output-dir', '-o',
              default='gemma_results',
              show_default=True,
              type=click.Path(),
              help='输出目录|Output directory')
@click.option('--n-pca',
              type=int,
              default=10,
              show_default=True,
              help='PCA主成分数|Number of PCA components')
@click.option('--threads',
              type=int,
              default=12,
              show_default=True,
              help='线程数|Number of threads')
@click.option('--gemma',
              default='~/.local/bin/gemma',
              show_default=True,
              help='GEMMA程序路径|GEMMA program path')
@click.option('--maf',
              type=float,
              help='最小等位基因频率阈值|Minor allele frequency threshold')
@click.option('--geno',
              type=float,
              help='SNP缺失率阈值|SNP missing rate threshold')
@click.option('--mind',
              type=float,
              help='样本缺失率阈值|Sample missing rate threshold')
@click.option('--hwe',
              type=float,
              help='Hardy-Weinberg p值阈值|Hardy-Weinberg p-value threshold')
@click.option('--no-qc',
              is_flag=True,
              help='跳过PLINK质控|Skip PLINK quality control')
@click.option('--lmm',
              type=int,
              help='LMM检验方法: 1=Wald, 2=LRT, 3=Score, 4=全部|LMM test method: 1=Wald, 2=LRT, 3=Score, 4=all')
@click.option('--gk',
              type=int,
              help='亲缘关系矩阵方法: 1=中心化, 2=标准化|Kinship matrix method: 1=centered, 2=standardized')
@click.option('--miss-gemma',
              type=float,
              help='GEMMA缺失率阈值|GEMMA missing rate threshold')
@click.option('--maf-gemma',
              type=float,
              help='GEMMA MAF阈值|GEMMA MAF threshold')
@click.option('--notsnp',
              is_flag=True,
              help='不输出每个SNP的估计值(更快)|Do not output estimated values for each SNP (faster)')
@click.option('--verbose', '-v',
              count=True,
              help='详细模式(-v: INFO, -vv: DEBUG)|Verbose mode (-v: INFO, -vv: DEBUG)')
@click.option('--quiet',
              is_flag=True,
              help='静默模式(仅ERROR)|Quiet mode (only ERROR)')
@click.option('--log-file',
              type=click.Path(),
              help='日志文件路径|Log file path')
def gemma_gwas(input, pheno, output_dir, n_pca, threads, gemma,
               maf, geno, mind, hwe, no_qc,
               lmm, gk, miss_gemma, maf_gemma, notsnp,
               verbose, quiet, log_file):
    """
    GEMMA GWAS批量分析|GEMMA GWAS Batch Analysis

    使用GEMMA进行全基因组关联分析，包括PCA、亲缘关系矩阵和LMM检验|Perform GWAS using GEMMA with PCA, kinship matrix, and LMM test

    示例|Examples: biopytools gemma-gwas -i genotype.vcf.gz -p phenotype.txt
    """

    # 延迟加载|Lazy load: import only when actually called
    gemma_gwas_main = _lazy_import_gemma_gwas_main()

    # 构建主函数的参数列表|Build argument list for main function
    args = ['gemma_gwas.py']

    # 必需参数|Required parameters
    args.extend(['--input', input])
    args.extend(['-p', pheno])

    # 可选参数|Optional parameters (add only when non-default)
    if output_dir != 'gemma_results':
        args.extend(['-o', output_dir])

    if n_pca != 10:
        args.extend(['--n-pca', str(n_pca)])

    if threads != 12:
        args.extend(['--threads', str(threads)])

    if gemma != '~/.local/bin/gemma':
        args.extend(['--gemma', gemma])

    # PLINK质控参数|PLINK QC parameters
    if maf is not None:
        args.extend(['--maf', str(maf)])
    if geno is not None:
        args.extend(['--geno', str(geno)])
    if mind is not None:
        args.extend(['--mind', str(mind)])
    if hwe is not None:
        args.extend(['--hwe', str(hwe)])

    if no_qc:
        args.append('--no-qc')

    # GEMMA参数|GEMMA parameters
    if lmm is not None:
        args.extend(['--lmm', str(lmm)])
    if gk is not None:
        args.extend(['--gk', str(gk)])
    if miss_gemma is not None:
        args.extend(['--miss-gemma', str(miss_gemma)])
    if maf_gemma is not None:
        args.extend(['--maf-gemma', str(maf_gemma)])
    if notsnp:
        args.append('--notsnp')

    # 日志选项|Logging options
    if verbose:
        args.extend(['-v'] * verbose)
    if quiet:
        args.append('--quiet')
    if log_file:
        args.extend(['--log-file', log_file])

    # 保存并恢复sys.argv|Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        # 调用原始主函数|Call original main function
        gemma_gwas_main()
    except SystemExit as e:
        # 处理正常程序退出|Handle normal program exit
        sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n用户中断|User interrupted", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
