"""
GCTB全基因组复杂性状贝叶斯分析|GCTB Genome-wide Complex Trait Bayesian Analysis
"""

import click
import sys
import os


def _lazy_import_gctb_main():
    """延迟加载gctb主函数|Lazy load gctb main function"""
    try:
        from ...gctb.main import main as gctb_main
        return gctb_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_file_exists(file_path):
    """验证文件存在(仅在非帮助模式)|Validate file existence (only in non-help mode)"""
    if not _is_help_request() and file_path and not os.path.exists(file_path):
        raise click.BadParameter(f"文件不存在|File does not exist: {file_path}")
    return file_path


@click.command(
    short_help='GCTB全基因组复杂性状贝叶斯分析工具|GCTB Bayesian analysis tool for complex traits',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('-i', '--vcf-file',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='VCF变异文件|VCF variant file')
@click.option('-p', '--pheno-file',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='表型文件或GWAS汇总统计文件|Phenotype file or GWAS summary statistics file')
@click.option('-o', '--output-dir',
              default='./gctb_output',
              show_default=True,
              help='输出目录|Output directory')
@click.option('--gctb-path',
              default='~/miniforge3/envs/gctb/bin/gctb',
              show_default=True,
              help='GCTB软件路径|GCTB software path')
@click.option('--plink-path',
              default='~/miniforge3/envs/Population_genetics/bin/plink',
              show_default=True,
              help='PLINK软件路径|PLINK software path')
@click.option('--maf-threshold',
              type=float,
              default=0.01,
              show_default=True,
              help='MAF阈值|MAF threshold')
@click.option('--miss-threshold',
              type=float,
              default=0.1,
              show_default=True,
              help='缺失率阈值|Missing rate threshold')
@click.option('--bayes-type',
              type=click.Choice(['S', 'R', 'C']),
              default='S',
              show_default=True,
              help='贝叶斯模型类型|Bayesian model type')
@click.option('--analysis-mode',
              type=click.Choice(['individual', 'summary']),
              default='individual',
              show_default=True,
              help='分析模式|Analysis mode')
@click.option('--ld-matrix-type',
              type=click.Choice(['sparse', 'block', 'eigen']),
              default='sparse',
              show_default=True,
              help='LD矩阵类型|LD matrix type')
@click.option('-t', '--threads',
              type=int,
              default=12,
              show_default=True,
              help='线程数|Number of threads')
@click.option('--batch/--no-batch',
              default=True,
              show_default=False,
              help='批量处理多个表型（默认开启）|Batch process multiple phenotypes (enabled by default)')
@click.option('--step',
              type=click.Choice(['convert', 'qc', 'freq', 'ld', 'analysis']),
              help='只运行指定步骤|Run only specified step')
def gctb(vcf_file, pheno_file, output_dir, gctb_path, plink_path,
         maf_threshold, miss_threshold, bayes_type, analysis_mode,
         ld_matrix_type, threads, batch, step):
    """
    GCTB全基因组复杂性状贝叶斯分析工具|GCTB Bayesian Analysis Tool for Complex Traits

    从VCF文件开始，自动完成数据转换、质控、LD矩阵构建和贝叶斯分析
    Automate data conversion, QC, LD matrix construction, and Bayesian analysis from VCF

    默认启用批量表型处理（自动检测并拆分多个表型）|Batch phenotype processing is enabled by default (auto-detect and split multiple phenotypes)

    示例|Example: biopytools gctb -i variants.vcf -p phenotype.txt -o results
    禁用批量处理|Disable batch: biopytools gctb -i variants.vcf -p phenotype.txt -o results --no-batch
    """

    # 延迟加载|Lazy loading
    gctb_main = _lazy_import_gctb_main()

    # 构建参数列表|Build argument list
    args = ['gctb.py']

    # 必需参数|Required parameters
    args.extend(['-i', vcf_file])
    args.extend(['-p', pheno_file])

    # 可选参数|Optional parameters
    if output_dir != './gctb_output':
        args.extend(['-o', output_dir])

    if gctb_path != '~/miniforge3/envs/gctb/bin/gctb':
        args.extend(['--gctb-path', gctb_path])

    if plink_path != '~/miniforge3/envs/Population_genetics/bin/plink':
        args.extend(['--plink-path', plink_path])

    if maf_threshold != 0.01:
        args.extend(['--maf-threshold', str(maf_threshold)])

    if miss_threshold != 0.1:
        args.extend(['--miss-threshold', str(miss_threshold)])

    if bayes_type != 'S':
        args.extend(['--bayes-type', bayes_type])

    if analysis_mode != 'individual':
        args.extend(['--analysis-mode', analysis_mode])

    if ld_matrix_type != 'sparse':
        args.extend(['--ld-matrix-type', ld_matrix_type])

    if threads != 12:
        args.extend(['--threads', str(threads)])

    if not batch:
        args.extend(['--no-batch'])

    if step:
        args.extend(['--step', step])

    # 执行主程序|Execute main program
    original_argv = sys.argv
    sys.argv = args

    try:
        gctb_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
