"""
CentIER着丝粒鉴定命令|CentIER Centromere Identification Command
"""

import click
import sys
import os


def _lazy_import_centier_main():
    """延迟加载centier主函数|Lazy load centier main function"""
    try:
        from ...centier.main import main as centier_main
        return centier_main
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
    short_help='CentIER着丝粒鉴定工具|CentIER centromere identification tool',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--input', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value),
              help='基因组FASTA文件|Genome FASTA file path')
@click.option('--output-dir', '-o',
              default='./centier_output',
              show_default=True,
              help='输出目录|Output directory')
@click.option('--centier-path',
              default='~/software/CentIER/CentIER-main',
              show_default=True,
              help='CentIER软件路径|CentIER software path')
@click.option('--gff',
              default=None,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='GFF/GTF注释文件|GFF/GTF annotation file (optional)')
@click.option('--matrix1',
              default=None,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='Hi-C矩阵文件(100000分辨率)|Hi-C matrix file at 100000 resolution (optional)')
@click.option('--matrix2',
              default=None,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='Hi-C矩阵文件(200000分辨率)|Hi-C matrix file at 200000 resolution (optional)')
@click.option('--bed1',
              default=None,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='Hi-C BED文件(对应matrix1)|Hi-C BED file for matrix1 (optional)')
@click.option('--bed2',
              default=None,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='Hi-C BED文件(对应matrix2)|Hi-C BED file for matrix2 (optional)')
@click.option('--threads', '-t',
              type=int,
              default=12,
              show_default=True,
              help='线程数|Number of threads')
@click.option('--kmer-size', '-k',
              type=int,
              default=21,
              show_default=True,
              help='K-mer大小|K-mer size')
@click.option('--center-tolerance', '-c',
              type=int,
              default=15,
              show_default=True,
              help='中心容差|Center tolerance')
@click.option('--step-len',
              type=int,
              default=10000,
              show_default=True,
              help='步长|Step length')
@click.option('--mul-cents',
              is_flag=True,
              help='保留所有潜在的着丝粒区域|Retain all potential centromeric regions')
@click.option('--mingap',
              type=int,
              default=2,
              show_default=True,
              help='最小Gap值|Minimum gap value n*100000')
@click.option('--signal-threshold',
              type=float,
              default=0.7,
              show_default=True,
              help='信号阈值|Signal threshold')
@click.option('--step', '-s',
              type=click.Choice(['1', '2', '3', '4', '5', '6']),
              help='运行指定步骤|Run only specified step (1-6)')
@click.option('--skip-dependency-check',
              is_flag=True,
              help='跳过依赖检查|Skip dependency check')
@click.option('--summary',
              is_flag=True,
              help='输出分析结果摘要|Output analysis result summary')
def centier(input, output_dir, centier_path, gff, matrix1, matrix2, bed1, bed2,
            threads, kmer_size, center_tolerance, step_len, mul_cents, mingap, signal_threshold,
            step, skip_dependency_check, summary):
    """
    CentIER着丝粒鉴定工具|CentIER Centromere Identification Tool

    用于T2T基因组组装的着丝粒识别和注释|Identify and annotate centromeric regions in T2T-assembled genomes

    示例|Example: biopytools centier -i genome.fa -o output_dir/
    """

    # 延迟加载|Lazy loading
    centier_main = _lazy_import_centier_main()

    # 构建参数列表|Build argument list
    args = ['centier.py']

    # 必需参数|Required parameters
    args.extend(['-i', input])
    args.extend(['-o', output_dir])

    # 可选路径参数|Optional path parameters
    if centier_path != '~/software/CentIER/CentIER-main':
        args.extend(['--centier-path', centier_path])

    if gff:
        args.extend(['--gff', gff])

    # Hi-C数据|Hi-C data
    if matrix1:
        args.extend(['--matrix1', matrix1])
    if matrix2:
        args.extend(['--matrix2', matrix2])
    if bed1:
        args.extend(['--bed1', bed1])
    if bed2:
        args.extend(['--bed2', bed2])

    # 参数|Parameters
    if threads != 12:
        args.extend(['--threads', str(threads)])
    if kmer_size != 21:
        args.extend(['--kmer-size', str(kmer_size)])
    if center_tolerance != 15:
        args.extend(['--center-tolerance', str(center_tolerance)])
    if step_len != 10000:
        args.extend(['--step-len', str(step_len)])
    if mingap != 2:
        args.extend(['--mingap', str(mingap)])
    if signal_threshold != 0.7:
        args.extend(['--signal-threshold', str(signal_threshold)])

    # 布尔选项|Boolean options
    if mul_cents:
        args.append('--mul-cents')
    if skip_dependency_check:
        args.append('--skip-dependency-check')
    if summary:
        args.append('--summary')

    # 步骤控制|Step control
    if step:
        args.extend(['--step', step])

    # 执行主程序|Execute main program
    original_argv = sys.argv
    sys.argv = args

    try:
        centier_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
