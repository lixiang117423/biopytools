"""
TASSEL GWAS命令|TASSEL GWAS Analysis Command
"""

import click
import sys
import os


def _lazy_import_tassel_gwas_main():
    """延迟加载tassel_gwas主函数|Lazy load tassel_gwas main function"""
    try:
        from ...tassel_gwas.main import main as tassel_gwas_main
        return tassel_gwas_main
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


def _validate_directory_exists(dir_path):
    """验证目录存在(如不存在则创建)|Validate directory exists (create if not exists)"""
    if not _is_help_request() and dir_path:
        if not os.path.exists(dir_path):
            try:
                os.makedirs(dir_path, exist_ok=True)
            except OSError:
                raise click.BadParameter(f"无法创建输出目录|Cannot create output directory: {dir_path}")
    return dir_path


@click.command(short_help="TASSEL GWAS分析工具|TASSEL GWAS Analysis Tool",
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--vcf', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='VCF基因型文件|VCF genotype file')
@click.option('--pheno', '-p',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='表型数据文件|Phenotype data file')
@click.option('--output', '-o',
              required=True,
              callback=lambda ctx, param, value: _validate_directory_exists(value) if value else None,
              help='输出目录|Output directory')
@click.option('--model', '-m',
              type=click.Choice(['GLM', 'MLM', 'BOTH']),
              default='MLM',
              show_default=True,
              help='GWAS模型|GWAS model (GLM/MLM/BOTH)')
@click.option('--memory',
              default='100g',
              show_default=True,
              help='Java最大内存|Java maximum memory (e.g., 100g, 200g)')
@click.option('--threads', '-t',
              type=int,
              default=12,
              show_default=True,
              help='并行线程数|Number of parallel threads')
@click.option('--maf',
              type=float,
              help='最小等位基因频率过滤|Minimum allele frequency filter')
@click.option('--miss',
              type=float,
              help='最大缺失率过滤|Maximum missing rate filter')
@click.option('--pca-components',
              type=int,
              default=5,
              show_default=True,
              help='PCA主成分数(MLM协变量)|Number of PCA components for MLM covariates')
@click.option('--q-matrix', '-q',
              type=click.Path(exists=True),
              help='群体结构Q矩阵文件|Population structure Q matrix file')
@click.option('--kinship', '-k',
              type=click.Path(exists=True),
              help='亲缘关系矩阵文件|Kinship matrix file')
@click.option('--tassel-path',
              type=click.Path(exists=True),
              help='TASSEL安装路径|TASSEL installation path')
@click.option('--skip-sort',
              is_flag=True,
              help='跳过VCF排序|Skip VCF sorting')
@click.option('--keep-temp',
              is_flag=True,
              help='保留临时文件|Keep temporary files')
@click.option('--parallel',
              is_flag=True,
              help='并行处理多个性状|Parallel process multiple traits')
@click.option('--workers',
              type=int,
              default=4,
              show_default=True,
              help='并行工作进程数|Number of parallel workers')
@click.option('--log-level',
              type=click.Choice(['DEBUG', 'INFO', 'WARN', 'ERROR']),
              default='INFO',
              show_default=True,
              help='日志级别|Log level')
def tassel_gwas(vcf, pheno, output, model, memory, threads, maf, miss, pca_components, q_matrix,
                kinship, tassel_path, skip_sort, keep_temp, parallel, workers, log_level):
    """
    TASSEL GWAS分析工具|TASSEL GWAS Analysis Tool

    使用TASSEL进行全基因组关联分析(GWAS)|Perform GWAS using TASSEL

    示例|Examples: biopytools tassel-gwas -i input.vcf.gz -p traits.txt -o results
    """

    # 延迟加载|Lazy load: import only when actually called
    tassel_gwas_main = _lazy_import_tassel_gwas_main()

    # 构建主函数的参数列表|Build argument list for main function
    args = ['tassel_gwas.py']

    # 必需参数|Required parameters
    args.extend(['--vcf', vcf])
    args.extend(['--pheno', pheno])
    args.extend(['--output', output])

    # 可选参数|Optional parameters
    args.extend(['--model', model])
    args.extend(['--memory', memory])
    args.extend(['--threads', str(threads)])
    args.extend(['--log-level', log_level])

    # 仅在非默认值时添加|Add only when non-default
    if maf is not None:
        args.extend(['--maf', str(maf)])

    if miss is not None:
        args.extend(['--miss', str(miss)])

    if pca_components != 5:
        args.extend(['--pca-components', str(pca_components)])

    if q_matrix:
        args.extend(['--q-matrix', q_matrix])

    if kinship:
        args.extend(['--kinship', kinship])

    if tassel_path:
        args.extend(['--tassel-path', tassel_path])

    if skip_sort:
        args.append('--skip-sort')

    if keep_temp:
        args.append('--keep-temp')

    if parallel:
        args.append('--parallel')

    if workers != 4:
        args.extend(['--workers', str(workers)])

    # 保存并恢复sys.argv|Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        # 调用原始主函数|Call original main function
        tassel_gwas_main()
    except SystemExit as e:
        # 处理正常程序退出|Handle normal program exit
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"运行时错误|Runtime Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
