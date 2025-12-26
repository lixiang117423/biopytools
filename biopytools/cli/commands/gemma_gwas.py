"""
GEMMA GWAS Batch Analysis Command
GEMMA全基因组关联分析批量分析命令
"""

import click
import sys
import os


def _lazy_import_gemma_gwas_main():
    """懒加载gemma_gwas main函数"""
    try:
        from ...gemma_gwas.main import main as gemma_gwas_main
        return gemma_gwas_main
    except ImportError as e:
        click.echo(f"Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否是帮助请求"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_file_exists(file_path):
    """验证文件是否存在（仅在非帮助模式下）"""
    if not _is_help_request() and not os.path.exists(file_path):
        raise click.BadParameter(f"File does not exist: {file_path}")
    return file_path


@click.command(short_help="GEMMA GWAS批量分析工具",
               context_settings=dict(help_option_names=['-h', '--help'],
                                    max_content_width=120))
@click.option('--vcf',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='VCF format genotype file (can be compressed)')
@click.option('--pheno', '-p',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='Phenotype file (first column: sample ID, has header)')
@click.option('--outdir', '-o',
              default='gemma_results',
              type=click.Path(),
              help='Output directory (default: gemma_results)')
@click.option('--n-pca',
              type=int,
              default=10,
              help='Number of PCA components (default: 10)')
@click.option('--threads',
              type=int,
              default=12,
              help='Number of threads (default: 12)')
@click.option('--gemma',
              default='/share/org/YZWL/yzwl_lixg/.local/bin/gemma',
              help='GEMMA program path')
@click.option('--maf',
              type=float,
              help='Minor allele frequency threshold (default: 0.05)')
@click.option('--geno',
              type=float,
              help='SNP missing rate threshold (default: 0.1)')
@click.option('--mind',
              type=float,
              help='Sample missing rate threshold (default: 0.1)')
@click.option('--hwe',
              type=float,
              help='Hardy-Weinberg p-value threshold (default: 1e-6)')
@click.option('--no-qc',
              is_flag=True,
              help='Skip PLINK quality control')
@click.option('--lmm',
              type=int,
              help='LMM test method: 1=Wald, 2=LRT, 3=Score, 4=all (default: 4)')
@click.option('--gk',
              type=int,
              help='Kinship matrix method: 1=centered, 2=standardized (default: 1)')
@click.option('--miss-gemma',
              type=float,
              help='GEMMA missing rate threshold (default: 0.05)')
@click.option('--maf-gemma',
              type=float,
              help='GEMMA MAF threshold (default: 0.01)')
@click.option('--notsnp',
              is_flag=True,
              help='Do not output estimated values for each SNP (faster)')
@click.option('--verbose', '-v',
              count=True,
              help='Verbose mode (-v: INFO, -vv: DEBUG)')
@click.option('--quiet',
              is_flag=True,
              help='Quiet mode (only ERROR)')
@click.option('--log-file',
              type=click.Path(),
              help='Log file path')
def gemma_gwas(vcf, pheno, outdir, n_pca, threads, gemma,
               maf, geno, mind, hwe, no_qc,
               lmm, gk, miss_gemma, maf_gemma, notsnp,
               verbose, quiet, log_file):
    """
    GEMMA GWAS批量分析工具

    自动化GWAS分析流程，支持多表型批量分析，包含质控、PCA、kinship计算
    和线性混合模型关联分析。

    示例 | Examples:

    \b
    # 基本使用
    biopytools gemma-gwas \\
        --vcf genotype.vcf.gz \\
        --pheno phenotype.txt

    \b
    # 指定输出目录和PCA数量
    biopytools gemma-gwas \\
        --vcf genotype.vcf.gz \\
        --pheno phenotype.txt \\
        -o results \\
        --n-pca 15

    \b
    # 完整参数示例
    biopytools gemma-gwas \\
        --vcf genotype.vcf.gz \\
        --pheno phenotype.txt \\
        -o results \\
        --n-pca 10 \\
        --maf 0.05 \\
        --geno 0.05 \\
        --lmm 4 \\
        --notsnp

    \b
    # 跳过质控
    biopytools gemma-gwas \\
        --vcf genotype.vcf.gz \\
        --pheno phenotype.txt \\
        --no-qc
    """

    # 懒加载：只有在实际调用时才导入模块
    gemma_gwas_main = _lazy_import_gemma_gwas_main()

    # 构建参数列表传递给原始main函数
    args = ['gemma_gwas.py']

    # 必需参数
    args.extend(['--vcf', vcf])
    args.extend(['-p', pheno])

    # 可选参数（只在非默认值时添加）
    if outdir != 'gemma_results':
        args.extend(['-o', outdir])

    if n_pca != 10:
        args.extend(['--n-pca', str(n_pca)])

    if threads != 12:
        args.extend(['--threads', str(threads)])

    if gemma != '/share/org/YZWL/yzwl_lixg/.local/bin/gemma':
        args.extend(['--gemma', gemma])

    # PLINK质控参数
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

    # GEMMA参数
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

    # 日志参数
    if verbose:
        args.extend(['-v'] * verbose)
    if quiet:
        args.append('--quiet')
    if log_file:
        args.extend(['--log-file', log_file])

    # 保存并恢复sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        # 调用原始的main函数
        gemma_gwas_main()
    except SystemExit as e:
        # 处理程序正常退出
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
