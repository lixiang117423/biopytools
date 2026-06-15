"""
Pixy群体遗传学统计命令|Pixy Population Genetics Statistics Command
"""

import click
import sys
import os


def _lazy_import_pixy_main():
    """延迟加载pixy主函数|Lazy load pixy main function"""
    try:
        from ...pixy.main import main as pixy_main
        return pixy_main
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
    short_help='Pixy群体遗传学统计工具|Pixy population genetics statistics tool',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('-i', '--vcf-file',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='VCF文件路径（需用bgzip压缩并建立tabix索引）|VCF file path (must be bgzip-compressed and tabix-indexed)')
@click.option('-p', '--pop-file',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='群体文件路径（两列：样本ID 群体名）|Population file path (two columns: sample_id population_name)')
@click.option('-o', '--output-dir',
              default='./pixy_output',
              show_default=True,
              type=click.Path(),
              help='输出目录|Output directory')
@click.option('--stats',
              default='pi,dxy,fst,watterson,tajima',
              show_default=True,
              help='要计算的统计量，逗号分隔|Statistics to calculate, comma-separated')
@click.option('--calc-pi',
              is_flag=True,
              help='计算pi（核苷酸多样性）|Calculate pi (nucleotide diversity)')
@click.option('--calc-dxy',
              is_flag=True,
              help='计算dxy（群体间核苷酸差异）|Calculate dxy (nucleotide divergence)')
@click.option('--calc-fst',
              is_flag=True,
              help='计算fst（遗传分化系数）|Calculate fst (genetic differentiation)')
@click.option('--calc-watterson-theta',
              is_flag=True,
              help='计算Watterson\'s theta|Calculate Watterson\'s theta')
@click.option('--calc-tajima-d',
              is_flag=True,
              help='计算Tajima\'s D|Calculate Tajima\'s D')
@click.option('-w', '--window-size',
              type=int,
              default=None,
              help='窗口大小bp（不设置则全基因组计算）|Window size in bp (null for genome-wide)')
@click.option('-b', '--bed-file',
              default=None,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='BED文件定义窗口|BED file defining windows')
@click.option('-s', '--sites-file',
              default=None,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='位点文件（只计算特定位点）|Sites file (calculate only specific sites)')
@click.option('--min-samples',
              type=int,
              default=0,
              show_default=True,
              help='每个群体最小样本数|Minimum samples per population')
@click.option('--max-missing',
              type=float,
              default=1.0,
              show_default=True,
              help='最大缺失率|Maximum missing rate')
@click.option('--min-maf',
              type=float,
              default=0.0,
              show_default=True,
              help='最小等位基因频率|Minor allele frequency')
@click.option('--zscore-window',
              type=int,
              default=None,
              help='Z-score过滤窗口大小|Z-score filtering window size')
@click.option('-c', '--chromosomes',
              default=None,
              help='指定染色体列表，逗号分隔|List of chromosomes, comma-separated')
@click.option('--pixy-path',
              default='pixy',
              show_default=True,
              help='pixy可执行文件路径|pixy executable path')
@click.option('--conda-env',
              default='~/miniforge3/envs/pixy_v.2.0.0',
              show_default=True,
              help='conda环境路径|conda environment path')
@click.option('-t', '--threads',
              type=int,
              default=12,
              show_default=True,
              help='线程数|Number of threads')
@click.option('--bypass-invariant-check',
              is_flag=True,
              help='强制绕过不变位点检查（默认自动检测VCF并自动绕过）|Force bypass invariant sites check (default: auto-detect VCF and bypass if needed)')
@click.option('--keep-intermediate',
              is_flag=True,
              help='保留中间文件|Keep intermediate files')
@click.option('-v', '--verbose',
              is_flag=True,
              help='详细输出|Verbose output')
def pixy(vcf_file, pop_file, output_dir, stats, calc_pi, calc_dxy, calc_fst,
         calc_watterson_theta, calc_tajima_d, window_size, bed_file, sites_file,
         min_samples, max_missing, min_maf, zscore_window, chromosomes,
         pixy_path, conda_env, threads, bypass_invariant_check, keep_intermediate, verbose):
    """
    Pixy群体遗传学统计工具|Pixy Population Genetics Statistics Tool

    计算群体遗传学统计量（pi、dxy、fst、Watterson's theta、Tajima's D）|Calculate population genetics statistics (pi, dxy, fst, Watterson's theta, Tajima's D)

    示例|Example: biopytools pixy -i variants.vcf.gz -p populations.txt -o pixy_output -w 100000

    注意|Note: pixy要求必须指定窗口大小(-w)、BED文件(-b)或位点文件(-s)|pixy requires window_size (-w), bed_file (-b), or sites_file (-s)
    """

    # 延迟加载|Lazy loading
    pixy_main = _lazy_import_pixy_main()

    # 构建参数列表|Build argument list
    args = ['pixy.py']

    # 必需参数|Required parameters
    args.extend(['-i', vcf_file])
    args.extend(['-p', pop_file])
    args.extend(['-o', output_dir])

    # 统计量选择|Statistics selection
    if stats != 'pi,dxy,fst,watterson,tajima':
        args.extend(['--stats', stats])

    if calc_pi:
        args.append('--calc-pi')

    if calc_dxy:
        args.append('--calc-dxy')

    if calc_fst:
        args.append('--calc-fst')

    if calc_watterson_theta:
        args.append('--calc-watterson-theta')

    if calc_tajima_d:
        args.append('--calc-tajima-d')

    # 窗口参数|Window parameters
    if window_size:
        args.extend(['-w', str(window_size)])

    if bed_file:
        args.extend(['-b', bed_file])

    if sites_file:
        args.extend(['-s', sites_file])

    # 质控参数|QC parameters
    if min_samples != 0:
        args.extend(['--min-samples', str(min_samples)])

    if max_missing != 1.0:
        args.extend(['--max-missing', str(max_missing)])

    if min_maf != 0.0:
        args.extend(['--min-maf', str(min_maf)])

    if zscore_window:
        args.extend(['--zscore-window', str(zscore_window)])

    # 染色体参数|Chromosome parameters
    if chromosomes:
        args.extend(['-c', chromosomes])

    # 工具路径|Tool paths
    if pixy_path != 'pixy':
        args.extend(['--pixy-path', pixy_path])

    if conda_env != '~/miniforge3/envs/pixy_v.2.0.0':
        args.extend(['--conda-env', conda_env])

    # 其他参数|Other parameters
    if threads != 12:
        args.extend(['-t', str(threads)])

    if bypass_invariant_check:
        args.append('--bypass-invariant-check')

    if keep_intermediate:
        args.append('--keep-intermediate')

    if verbose:
        args.append('--verbose')

    # 执行主程序|Execute main program
    original_argv = sys.argv
    sys.argv = args

    try:
        pixy_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
