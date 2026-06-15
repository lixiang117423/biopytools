"""
RNA-Bloom转录组从头组装命令|RNA-Bloom De Novo Transcriptome Assembly Command
"""

import click
import sys
import os


def _lazy_import_rnabloom_main():
    """延迟加载rnabloom主函数|Lazy load rnabloom main function"""
    try:
        from ...rnabloom.main import main as rnabloom_main
        return rnabloom_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_file_exists(file_path):
    """验证文件存在(仅在非帮助模式)|Validate file existence (only in non-help mode)"""
    if file_path and not _is_help_request() and not os.path.exists(file_path):
        raise click.BadParameter(f"文件不存在|File does not exist: {file_path}")
    return file_path


@click.command(
    short_help='RNA-Bloom转录组从头组装工具|RNA-Bloom de novo transcriptome assembly tool',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
# 输入文件|Input files
@click.option('--left', '-1',
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='左端reads文件|Left reads file (paired-end)')
@click.option('--right', '-2',
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='右端reads文件|Right reads file (paired-end)')
@click.option('--sef',
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='单端正向reads文件|Single-end forward reads file')
@click.option('--ser',
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='单端反向reads文件|Single-end reverse reads file')
@click.option('--long',
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='长reads文件|Long reads file (ONT/PacBio)')
@click.option('--cell-list',
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='单细胞列表文件|Single-cell list file (pooled assembly)')
@click.option('-o', '--output-dir',
              required=True,
              type=click.Path(),
              help='输出目录|Output directory')
# 可选参数|Optional parameters
@click.option('-t', '--threads',
              default=12,
              show_default=True,
              type=int,
              help='线程数|Number of threads')
@click.option('--rnabloom-path',
              default='rnabloom',
              show_default=True,
              help='RNA-Bloom工具路径|RNA-Bloom tool path')
# Bloom filter配置|Bloom filter configuration
@click.option('--mem',
              type=float,
              help='Bloom filter总大小(GB)|Total Bloom filter size in GB')
@click.option('--fpr',
              type=float,
              help='假阳性率|False positive rate (0-1)')
@click.option('--nk', '--num-kmers',
              type=int,
              help='唯一k-mer数量|Number of unique kmers')
# 数据类型配置|Data type configuration
@click.option('--stranded',
              is_flag=True,
              help='链特异性数据|Strand-specific data')
@click.option('--revcomp-left',
              is_flag=True,
              help='反向互补左端reads|Reverse-complement left reads')
@click.option('--revcomp-right',
              is_flag=True,
              help='反向互补右端reads|Reverse-complement right reads')
@click.option('--pacbio',
              is_flag=True,
              help='PacBio数据|PacBio data (default: ONT)')
# 参考引导|Reference-guided
@click.option('--ref', '--reference',
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='参考转录本文件|Reference transcript file')
# 输出选项|Output options
@click.option('--min-length',
              default=200,
              show_default=True,
              type=int,
              help='最小转录本长度|Minimum transcript length')
@click.option('--uracil',
              is_flag=True,
              help='输出尿嘧啶(U)而非胸腺嘧啶(T)|Write uracil (U) instead of thymine (T)')
@click.option('--no-nr',
              is_flag=True,
              help='不导出去冗余转录本|Do not export non-redundant transcripts')
# 步骤控制|Step control
@click.option('--stage',
              type=click.Choice(['1', '2', '3']),
              help='停止阶段|Stop at stage (1-3)')
def rnabloom(left, right, sef, ser, long, cell_list, output_dir, threads,
             rnabloom_path, mem, fpr, nk, stranded, revcomp_left, revcomp_right,
             pacbio, ref, min_length, uracil, no_nr, stage):
    """
    RNA-Bloom转录组从头组装工具|RNA-Bloom De Novo Transcriptome Assembly Tool

    使用RNA-Bloom进行无参考转录组组装|Perform de novo transcriptome assembly using RNA-Bloom

    示例|Example: biopytools rnabloom --left reads_1.fq --right reads_2.fq -o ./results
    """

    # 延迟加载|Lazy loading
    rnabloom_main = _lazy_import_rnabloom_main()

    # 构建参数列表|Build argument list
    args = ['rnabloom.py']

    # 输入文件|Input files
    if left:
        args.extend(['--left', left])
    if right:
        args.extend(['--right', right])
    if sef:
        args.extend(['--sef', sef])
    if ser:
        args.extend(['--ser', ser])
    if long:
        args.extend(['--long', long])
    if cell_list:
        args.extend(['--cell-list', cell_list])

    # 输出|Output
    args.extend(['-o', output_dir])

    # 可选参数|Optional parameters
    if threads != 12:
        args.extend(['-t', str(threads)])

    if rnabloom_path != 'rnabloom':
        args.extend(['--rnabloom-path', rnabloom_path])

    # Bloom filter|Bloom filter
    if mem:
        args.extend(['--mem', str(mem)])
    if fpr:
        args.extend(['--fpr', str(fpr)])
    if nk:
        args.extend(['--nk', str(nk)])

    # 数据类型|Data type
    if stranded:
        args.append('--stranded')
    if revcomp_left:
        args.append('--revcomp-left')
    if revcomp_right:
        args.append('--revcomp-right')
    if pacbio:
        args.append('--pacbio')

    # 参考|Reference
    if ref:
        args.extend(['--ref', ref])

    # 输出选项|Output options
    if min_length != 200:
        args.extend(['--min-length', str(min_length)])
    if uracil:
        args.append('--uracil')
    if no_nr:
        args.append('--no-nr')

    # 步骤控制|Step control
    if stage:
        args.extend(['--stage', stage])

    # 执行主程序|Execute main program
    original_argv = sys.argv
    sys.argv = args

    try:
        exit_code = rnabloom_main()
        sys.exit(exit_code)
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
