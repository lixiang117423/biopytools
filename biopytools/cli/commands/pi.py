"""
核苷酸多样性计算命令|Nucleotide Diversity Calculation Command
"""

import click
import sys
import os


def _lazy_import_main():
    """延迟加载主函数|Lazy load main function"""
    try:
        from ...pi.main import main as pi_main
        return pi_main
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
    short_help='核苷酸多样性计算工具|Nucleotide diversity calculation tool',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('-i', '--vcf-file',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='VCF文件路径（bgzip压缩+tabix索引）|VCF file path (bgzip-compressed + tabix-indexed)')
@click.option('-p', '--pop-file',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='群体文件路径（样本ID 群体名）|Population file path (sample_id population_name)')
@click.option('-g', '--genome',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='参考基因组fasta文件路径（需有.fai索引）|Reference genome fasta path (requires .fai index)')
@click.option('-o', '--output-dir',
              default='./pi_output',
              show_default=True,
              type=click.Path(),
              help='输出目录|Output directory')
@click.option('-w', '--window-size',
              type=int,
              default=None,
              help='窗口大小bp（不设置则全基因组计算）|Window size in bp (omit for genome-wide)')
@click.option('--window-step',
              type=int,
              default=None,
              help='窗口步长bp（默认等于窗口大小，无重叠）|Window step in bp (default=window_size, no overlap)')
@click.option('--default-window-size',
              type=int,
              default=100000,
              show_default=True,
              help='自动滑窗默认窗口大小bp（全基因组模式时同步运行）|Default windowed window size in bp')
@click.option('--default-window-step',
              type=int,
              default=100000,
              show_default=True,
              help='自动滑窗默认步长bp|Default windowed step in bp')
@click.option('--maf',
              type=float,
              default=0.0,
              show_default=True,
              help='最小等位基因频率|Minor allele frequency')
@click.option('--max-missing',
              type=float,
              default=1.0,
              show_default=True,
              help='最大缺失率|Maximum missing rate')
@click.option('-t', '--threads',
              type=int,
              default=12,
              show_default=True,
              help='线程数|Number of threads')
@click.option('--keep-intermediate',
              is_flag=True,
              help='保留中间文件|Keep intermediate files')
@click.option('--vcftools-path',
              default=None,
              help='vcftools路径|vcftools path')
def pi(vcf_file, pop_file, genome, output_dir, window_size, window_step,
       default_window_size, default_window_step, maf, max_missing,
       threads, keep_intermediate, vcftools_path):
    """
    核苷酸多样性(pi)计算工具|Nucleotide Diversity (pi) Calculation Tool

    使用vcftools计算群体内核苷酸多样性，输出汇总表
    Calculate within-population nucleotide diversity using vcftools, output summary table

    示例|Examples: biopytools pi -i variants.vcf.gz -p populations.txt -g reference.fasta -o pi_output
    """
    pi_main = _lazy_import_main()

    # 构建参数列表|Build argument list
    args = ['pi.py']

    # 必需参数|Required parameters
    args.extend(['-i', vcf_file])
    args.extend(['-p', pop_file])
    args.extend(['-g', genome])
    if output_dir != './pi_output':
        args.extend(['-o', output_dir])

    # 窗口参数|Window parameters
    if window_size is not None:
        args.extend(['-w', str(window_size)])

    if window_step is not None:
        args.extend(['--window-step', str(window_step)])

    # 默认滑窗参数|Default windowed parameters
    if default_window_size != 100000:
        args.extend(['--default-window-size', str(default_window_size)])

    if default_window_step != 100000:
        args.extend(['--default-window-step', str(default_window_step)])

    # 质控参数|QC parameters
    if maf != 0.0:
        args.extend(['--maf', str(maf)])

    if max_missing != 1.0:
        args.extend(['--max-missing', str(max_missing)])

    # 其他参数|Other parameters
    if threads != 12:
        args.extend(['-t', str(threads)])

    if keep_intermediate:
        args.append('--keep-intermediate')

    if vcftools_path:
        args.extend(['--vcftools-path', vcftools_path])

    # 执行主程序|Execute main program
    original_argv = sys.argv
    sys.argv = args

    try:
        pi_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
