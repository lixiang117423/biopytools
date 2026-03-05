"""
基因组分析命令|Genome Analysis Command
"""

import click
import sys
import os


def _lazy_import_genome_analysis_main():
    """延迟加载genome_analysis主函数|Lazy load genome_analysis main function"""
    try:
        from ...genome_analysis.main import main as genome_analysis_main
        return genome_analysis_main
    except ImportError as e:
        click.echo(f"Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_input(input_path):
    """验证输入路径(文件或目录)|Validate input path (file or directory)"""
    if _is_help_request():
        return input_path
    if not os.path.exists(input_path):
        raise click.BadParameter(f"Path does not exist: {input_path}")

    # 如果是目录，直接返回
    if os.path.isdir(input_path):
        return input_path

    # 如果是文件，检查是否为FASTQ文件
    if os.path.isfile(input_path):
        fastq_extensions = ['.fastq', '.fq', '.fastq.gz', '.fq.gz']
        if not any(input_path.endswith(ext) for ext in fastq_extensions):
            raise click.BadParameter(
                f"File must be a FASTQ file (*.fastq, *.fq, *.fastq.gz, *.fq.gz): {input_path}"
            )
        return input_path

    raise click.BadParameter(f"Path is neither a file nor a directory: {input_path}")


@click.command(
    short_help='基因组分析|Genome Analysis',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--input', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_input(value) if value else None,
              help='输入FASTQ文件或目录|Input FASTQ file or directory')
@click.option('--output-dir', '-o',
              required=True,
              type=click.Path(),
              help='输出目录|Output directory')
@click.option('--read-length', '-l',
              type=int,
              default=150,
              show_default=True,
              help='测序读长|Read length')
@click.option('--kmer-size', '-k',
              type=int,
              default=21,
              show_default=True,
              help='K-mer大小|K-mer size')
@click.option('--threads', '-t',
              type=int,
              default=12,
              show_default=True,
              help='线程数|Number of threads')
@click.option('--hash-size', '-s',
              default='10G',
              show_default=True,
              help='Jellyfish哈希表大小|Jellyfish hash size')
@click.option('--max-kmer-cov', '-c',
              type=int,
              default=1000,
              show_default=True,
              help='最大k-mer覆盖度|Max k-mer coverage')
@click.option('--genomescope-r',
              default='~/software/scripts/genomescope.R',
              show_default=True,
              help='GenomeScope R脚本路径|GenomeScope R script path')
@click.option('--skip-smudgeplot',
              is_flag=True,
              default=False,
              help='跳过Smudgeplot倍性分析|Skip Smudgeplot ploidy analysis')
@click.option('--fastk-table',
              default='',
              help='FastK表文件路径|FastK table file path')
@click.option('--fastk-memory',
              default='100G',
              show_default=True,
              help='FastK内存大小|FastK memory size')
@click.option('--read1-suffix',
              default='*_1.clean.fq.gz',
              show_default=True,
              help='Read1文件后缀模式|Read1 file suffix pattern')
def genome_analysis(input, output_dir, read_length, kmer_size, threads,
                   hash_size, max_kmer_cov, genomescope_r, skip_smudgeplot,
                   fastk_table, fastk_memory, read1_suffix):
    """
    基因组分析工具|Genome Analysis Tool

    使用Jellyfish和GenomeScope进行k-mer分析和基因组特征评估|Perform k-mer analysis and genome characterization using Jellyfish and GenomeScope

    示例|Example: biopytools genome-analysis -i fastq_dir -o output_dir -l 150

    """

    # 延迟加载|Lazy loading: import only when actually called
    genome_analysis_main = _lazy_import_genome_analysis_main()

    # 构建参数列表|Build argument list
    args = ['genome_analysis.py']

    # 必需参数|Required parameters
    args.extend(['-i', input])
    args.extend(['-o', output_dir])
    args.extend(['-l', str(read_length)])

    # 可选参数|Optional parameters
    if kmer_size != 21:
        args.extend(['-k', str(kmer_size)])

    if threads != 64:
        args.extend(['-t', str(threads)])

    if hash_size != '10G':
        args.extend(['-s', hash_size])

    if max_kmer_cov != 1000:
        args.extend(['-c', str(max_kmer_cov)])

    if genomescope_r != '~/software/scripts/genomescope.R':
        args.extend(['--genomescope-r', genomescope_r])

    if skip_smudgeplot:
        args.append('--skip-smudgeplot')

    if fastk_table:
        args.extend(['--fastk-table', fastk_table])

    if fastk_memory != '16G':
        args.extend(['--fastk-memory', fastk_memory])

    if read1_suffix != '*_1.clean.fq.gz':
        args.extend(['--read1-suffix', read1_suffix])

    # 保存和恢复sys.argv|Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        # 调用原始主函数|Call original main function
        genome_analysis_main()
    except SystemExit as e:
        # 处理正常程序退出|Handle normal program exit
        if e.code != 0:
            sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n用户中断|Interrupted by user", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"执行失败|Execution failed: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
