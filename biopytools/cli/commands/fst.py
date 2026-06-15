"""
Fst变异注释|Fst Calculation Command

"""

import click
import sys
import os


def _lazy_import_fst_main():
    """延迟加载fst主函数|Lazy load fst main function"""
    try:
        from ...fst.main import main as fst_main
        return fst_main
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
    short_help='Fst计算工具|Fst calculation tool',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('-i', '--vcf-file',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='VCF文件路径|VCF file path')
@click.option('-p', '--pop-file',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='群体文件路径（样本ID + 群体标签）|Population file path (sample ID + population label)')
@click.option('-o', '--output-dir',
              default='./fst_output',
              show_default=True,
              type=click.Path(),
              help='输出目录|Output directory')
@click.option('--plink-path',
              default=None,
              help='PLINK软件路径|PLINK software path (default: auto-detect)')
@click.option('--maf',
              type=float,
              default=0.05,
              show_default=True,
              help='最小等位基因频率阈值|Minor allele frequency threshold')
@click.option('--geno',
              type=float,
              default=0.1,
              show_default=True,
              help='位点缺失率阈值|Genotype missing rate threshold')
@click.option('--mind',
              type=float,
              default=0.1,
              show_default=True,
              help='样本缺失率阈值|Sample missing rate threshold')
@click.option('--hwe',
              type=float,
              default=1e-6,
              show_default=True,
              help='Hardy-Weinberg平衡p值阈值|Hardy-Weinberg equilibrium p-value threshold')
@click.option('--enable-qc',
              is_flag=True,
              help='启用质控过滤（默认禁用）|Enable quality control filtering')
@click.option('--no-keep-intermediate',
              is_flag=True,
              help='不保留中间文件|Do not keep intermediate files')
@click.option('--enable-bootstrap',
              is_flag=True,
              help='启用bootstrap抽样|Enable bootstrap sampling')
@click.option('--bootstrap-iterations',
              type=int,
              default=100,
              show_default=True,
              help='Bootstrap迭代次数|Bootstrap iterations')
@click.option('--min-samples',
              type=int,
              default=10,
              show_default=True,
              help='最小样本数阈值|Minimum sample count threshold')
@click.option('--exclude-pops',
              default=None,
              help='手动排除群体（逗号分隔）|Manually exclude populations (comma-separated)')
@click.option('--no-ld-prune',
              is_flag=True,
              help='禁用LD pruning|Disable LD pruning')
@click.option('--ld-window',
              type=int,
              default=50,
              show_default=True,
              help='LD pruning窗口大小|LD pruning window size')
@click.option('--ld-step',
              type=int,
              default=10,
              show_default=True,
              help='LD pruning步长|LD pruning step size')
@click.option('--ld-r2',
              type=float,
              default=0.2,
              show_default=True,
              help='LD pruning R2阈值|LD pruning R2 threshold')
@click.option('--thin',
              type=float,
              default=None,
              help='SNP抽稀比例|SNP thinning ratio')
@click.option('-t', '--threads',
              type=int,
              default=12,
              show_default=True,
              help='线程数|Number of threads')
def fst(vcf_file, pop_file, output_dir, plink_path, maf, geno, mind, hwe,
        enable_qc, no_keep_intermediate, enable_bootstrap, bootstrap_iterations,
        min_samples, exclude_pops, no_ld_prune, ld_window, ld_step,
        ld_r2, thin, threads):
    """
    Fst计算工具|Fst Calculation Tool

    计算群体间遗传分化系数(Fst)，支持VCF输入和灵活的质量控制|Calculate genetic differentiation coefficient (Fst) between populations, supporting VCF input and flexible quality control

    示例|Example: biopytools fst -i variants.vcf -p population.txt -o fst_output
    """

    # 延迟加载|Lazy loading
    fst_main = _lazy_import_fst_main()

    # 构建参数列表|Build argument list
    args = ['fst.py']

    # 必需参数|Required parameters
    args.extend(['-i', vcf_file])
    args.extend(['-p', pop_file])
    args.extend(['-o', output_dir])

    # 可选参数|Optional parameters
    if plink_path:
        args.extend(['--plink-path', plink_path])

    if maf != 0.05:
        args.extend(['--maf', str(maf)])

    if geno != 0.1:
        args.extend(['--geno', str(geno)])

    if mind != 0.1:
        args.extend(['--mind', str(mind)])

    if hwe != 1e-6:
        args.extend(['--hwe', str(hwe)])

    if enable_qc:
        args.append('--enable-qc')

    if no_keep_intermediate:
        args.append('--no-keep-intermediate')

    if enable_bootstrap:
        args.append('--enable-bootstrap')

    if bootstrap_iterations != 100:
        args.extend(['--bootstrap-iterations', str(bootstrap_iterations)])

    if min_samples != 10:
        args.extend(['--min-samples', str(min_samples)])

    if exclude_pops:
        args.extend(['--exclude-pops', exclude_pops])

    if no_ld_prune:
        args.append('--no-ld-prune')

    if ld_window != 50:
        args.extend(['--ld-window', str(ld_window)])

    if ld_step != 10:
        args.extend(['--ld-step', str(ld_step)])

    if ld_r2 != 0.2:
        args.extend(['--ld-r2', str(ld_r2)])

    if thin is not None:
        args.extend(['--thin', str(thin)])

    if threads != 12:
        args.extend(['--threads', str(threads)])

    # 执行主程序|Execute main program
    original_argv = sys.argv
    sys.argv = args

    try:
        fst_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
