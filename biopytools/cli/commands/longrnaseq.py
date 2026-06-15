"""
三代转录组比对命令|Long RNA-seq Alignment Command
"""

import click
import sys
import os


def _lazy_import_longrnaseq_main():
    """延迟加载longrnaseq主函数|Lazy load longrnaseq main function"""
    try:
        from ...longrnaseq.main import main as longrnaseq_main
        return longrnaseq_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_path_exists(path):
    """验证路径存在(仅在非帮助模式)|Validate path exists (only in non-help mode)"""
    if not _is_help_request() and path and not os.path.exists(path):
        raise click.BadParameter(f"路径不存在|Path does not exist: {path}")
    # 允许文件夹路径|Allow directory path
    return path


@click.command(short_help="三代转录组比对|Long RNA-seq Alignment",
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('-i', '--input-file',
              required=True,
              callback=lambda ctx, param, value: _validate_path_exists(value),
              help='输入文件或文件夹（BAM/FASTQ）|Input file or directory (BAM/FASTQ)')
@click.option('-r', '--ref-genome',
              required=True,
              callback=lambda ctx, param, value: _validate_path_exists(value),
              help='参考基因组文件|Reference genome file')
@click.option('-o', '--output-dir',
              required=True,
              help='输出目录|Output directory')
@click.option('-s', '--sample-name',
              required=False,
              default=None,
              help='样本名称(可选，默认从输入文件名提取)|Sample name (optional, auto-extracted from input filename)')
@click.option('-t', '--threads',
              type=int,
              default=12,
              show_default=True,
              help='线程数|Number of threads')
@click.option('--max-intron',
              type=int,
              default=100000,
              show_default=True,
              help='最大intron长度|Maximum intron length')
@click.option('--min-mapq',
              type=int,
              default=20,
              show_default=True,
              help='最小mapping quality|Minimum mapping quality')
@click.option('--no-secondary',
              is_flag=True,
              default=False,
              help='不输出次优比对|Do not output secondary alignments')
@click.option('--minimap2-path',
              default='minimap2',
              show_default=True,
              help='minimap2可执行文件路径|minimap2 executable path')
@click.option('--samtools-path',
              default='samtools',
              show_default=True,
              help='samtools可执行文件路径|samtools executable path')
def longrnaseq(input_file, ref_genome, output_dir, sample_name, threads,
               max_intron, min_mapq, no_secondary, minimap2_path, samtools_path):
    """
    三代转录组比对工具|Long RNA-seq Alignment Tool

    使用minimap2将三代转录组数据比对到参考基因组
    Align long-read RNA-seq data to reference genome using minimap2

    示例|Examples: biopytools longrnaseq -i input.bam -r genome.fa -o output_dir
    """

    # 延迟加载|Lazy load: import only when actually called
    longrnaseq_main_func = _lazy_import_longrnaseq_main()

    # 构建主函数的参数列表|Build argument list for main function
    args = ['longrnaseq.py']

    # 必需参数|Required parameters
    args.extend(['-i', input_file])
    args.extend(['-r', ref_genome])
    args.extend(['-o', output_dir])
    args.extend(['-s', sample_name])

    # 可选参数|Optional parameters
    args.extend(['-t', str(threads)])
    args.extend(['--max-intron', str(max_intron)])
    args.extend(['--min-mapq', str(min_mapq)])
    args.extend(['--minimap2-path', minimap2_path])
    args.extend(['--samtools-path', samtools_path])

    if no_secondary:
        args.append('--no-secondary')

    # 保存并恢复sys.argv|Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        # 调用原始主函数|Call original main function
        longrnaseq_main_func()
    except SystemExit as e:
        # 处理正常程序退出|Handle normal program exit
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
