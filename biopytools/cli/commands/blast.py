"""
BLAST序列比对命令|BLAST Alignment Command
多序列文件与目标序列的BLAST比对分析完整流程|Complete pipeline for BLAST alignment analysis
"""

import click
import sys
import os


def _lazy_import_blast_main():
    """懒加载blast main函数|Lazy load blast main function"""
    try:
        from ...blast.main import main as blast_main
        return blast_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否是帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_file_exists(file_path):
    """验证文件是否存在（仅在非帮助模式下）| Validate file existence (only in non-help mode)"""
    if not _is_help_request() and file_path and not os.path.exists(file_path):
        raise click.BadParameter(f"文件不存在|File does not exist: {file_path}")
    return file_path


@click.command(short_help="BLAST序列比对分析",
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--input', '-i',
              type=click.Path(exists=True),
              help='输入文件或目录路径|Input file or directory path')
@click.option('--sample-map-file', '-s',
              type=click.Path(exists=True),
              help='样品映射文件，格式: 文件路径<TAB>样品名|Sample mapping file, format: filepath<TAB>samplename')
@click.option('--reference', '-r',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              type=click.Path(exists=True),
              help='目标基因序列文件|Target gene sequence file')
@click.option('--output', '-o',
              default='./blast_output',
              type=click.Path(),
              help='输出目录路径 (默认: ./blast_output)|Output directory (default: ./blast_output)')
@click.option('--blast-type',
              default='blastn',
              type=click.Choice(['blastn', 'blastp', 'blastx', 'tblastn', 'tblastx']),
              help='BLAST程序类型 (默认: blastn)|BLAST program type (default: blastn)')
@click.option('--evalue', '-e',
              default=1e-5,
              type=float,
              help='E-value阈值 (默认: 1e-5)|E-value threshold (default: 1e-5)')
@click.option('--max-target-seqs',
              default=10,
              type=int,
              help='最大目标序列数 (默认: 10)|Maximum target sequences (default: 10)')
@click.option('--word-size',
              default=11,
              type=int,
              help='词大小，仅适用于blastn/tblastx (默认: 11)|Word size for blastn/tblastx only (default: 11)')
@click.option('--threads', '-j',
              default=4,
              type=int,
              help='线程数 (默认: 4)|Number of threads (default: 4)')
@click.option('--input-suffix',
              default='*.fa',
              type=str,
              help='输入文件后缀模式 (默认: *.fa)|Input file suffix pattern (default: *.fa)')
@click.option('--target-db-type',
              default='nucl',
              type=click.Choice(['nucl', 'prot']),
              help='目标数据库类型 (默认: nucl)|Target database type (default: nucl)')
@click.option('--min-identity',
              default=70.0,
              type=float,
              help='最小序列相似度% (默认: 70.0)|Minimum sequence identity%% (default: 70.0)')
@click.option('--min-coverage',
              default=50.0,
              type=float,
              help='最小覆盖度% (默认: 50.0)|Minimum coverage%% (default: 50.0)')
@click.option('--high-quality-evalue',
              default=1e-10,
              type=float,
              help='高质量比对E-value阈值 (默认: 1e-10)|High quality alignment E-value threshold (default: 1e-10)')
@click.option('--auto-detect-samples',
              is_flag=True,
              default=True,
              help='自动检测样品名称 (仅在使用-i时)|Auto-detect sample names (only when using -i)')
@click.option('--sample-name-pattern',
              default=r'([^/]+?)(?:\.fa|\.fasta|\.fna)?$',
              type=str,
              help='样品名提取正则表达式|Sample name extraction regex')
@click.option('--makeblastdb-path',
              default='makeblastdb',
              type=str,
              help='makeblastdb程序路径 (默认: makeblastdb)|makeblastdb program path (default: makeblastdb)')
@click.option('--blastn-path',
              default='blastn',
              type=str,
              help='blastn程序路径 (默认: blastn)|blastn program path (default: blastn)')
@click.option('--blastp-path',
              default='blastp',
              type=str,
              help='blastp程序路径 (默认: blastp)|blastp program path (default: blastp)')
@click.option('--blastx-path',
              default='blastx',
              type=str,
              help='blastx程序路径 (默认: blastx)|blastx program path (default: blastx)')
@click.option('--tblastn-path',
              default='tblastn',
              type=str,
              help='tblastn程序路径 (默认: tblastn)|tblastn program path (default: tblastn)')
@click.option('--tblastx-path',
              default='tblastx',
              type=str,
              help='tblastx程序路径 (默认: tblastx)|tblastx program path (default: tblastx)')
@click.option('--alignment-output',
              default='both',
              type=click.Choice(['none', 'text', 'html', 'both']),
              help='比对可视化输出格式 (默认: both)|Alignment visualization output format (default: both)')
@click.option('--alignment-width',
              default=80,
              type=int,
              help='比对每行显示的字符数 (默认: 80)|Characters per line in alignment display (default: 80)')
@click.option('--alignment-min-identity',
              default=0.0,
              type=float,
              help='比对可视化最小相似度% (默认: 0.0)|Minimum identity for alignment visualization%% (default: 0.0)')
@click.option('--alignment-min-coverage',
              default=0.0,
              type=float,
              help='比对可视化最小覆盖度% (默认: 0.0)|Minimum coverage for alignment visualization%% (default: 0.0)')
@click.option('--alignment-max-per-sample',
              default=100,
              type=int,
              help='每个样品最多显示的比对数 (默认: 100)|Maximum alignments to display per sample (default: 100)')
@click.option('--html-theme',
              default='modern',
              type=click.Choice(['modern', 'classic', 'dark']),
              help='HTML主题样式 (默认: modern)|HTML theme style (default: modern)')
@click.option('--verbose', '-v',
              count=True,
              help='详细输出模式 (多次使用增加详细程度)|Verbose output mode (use multiple times for more verbosity)')
@click.option('--quiet',
              is_flag=True,
              help='静默模式，仅输出错误信息|Quiet mode, only output error messages')
@click.option('--log-level',
              type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL']),
              default='INFO',
              help='日志级别 (默认: INFO)|Log level (default: INFO)')
@click.option('--log-file',
              type=click.Path(),
              help='日志文件路径|Log file path')
@click.option('--force', '-f',
              is_flag=True,
              help='强制覆盖已存在的结果|Force overwrite existing results')
@click.option('--dry-run',
              is_flag=True,
              help='测试模式，不执行实际命令|Test mode, do not execute actual commands')
@click.option('--version', '-V',
              is_flag=True,
              help='显示版本信息|Show version information')
def blast(version, input, sample_map_file, reference, output, blast_type, evalue,
          max_target_seqs, word_size, threads, input_suffix, target_db_type,
          min_identity, min_coverage, high_quality_evalue, auto_detect_samples,
          sample_name_pattern, makeblastdb_path, blastn_path, blastp_path,
          blastx_path, tblastn_path, tblastx_path, alignment_output,
          alignment_width, alignment_min_identity, alignment_min_coverage,
          alignment_max_per_sample, html_theme, verbose, quiet,
          log_level, log_file, force, dry_run):
    """
    BLAST序列比对分析工具|BLAST Sequence Alignment Analysis Tool

    基于BLAST+的多序列比对分析，支持批量序列与目标序列的比对分析
    支持5种BLAST程序类型，自动样品检测，结果可视化和统计

    使用示例|Examples:

    \b
    # 基本用法 - 分析目录中的所有序列
    biopytools blast -i sequences/ -r nlr_genes.fa -o blast_results

    \b
    # 使用自定义E-value和相似度阈值
    biopytools blast -i sequences/ -r targets.fa -o results -e 1e-10 --min-identity 80

    \b
    # 使用样品映射文件
    biopytools blast -s sample_map.txt -r nlr_genes.fa -o results

    \b
    # 使用tblastx进行翻译比对
    biopytools blast -i proteins.fa -r genome.fa -o results --blast-type tblastx

    \b
    # 测试模式
    biopytools blast -i sequences/ -r targets.fa -o results --dry-run
    """

    # 处理版本信息|Handle version information
    if version:
        from ...blast import __version__
        click.echo(f"BLAST Analysis Pipeline Version: {__version__}")
        sys.exit(0)

    # 懒加载：只有在实际调用时才导入模块|Lazy loading: import only when actually called
    blast_main_func = _lazy_import_blast_main()

    # 构建参数列表传递给原始main函数|Build argument list for original main function
    args = ['blast.py']

    # 必需参数|Required parameters
    args.extend(['-r', reference])

    # 可选参数（只在非默认值时添加）| Optional parameters (add only when non-default)
    if input:
        args.extend(['--input', input])

    if sample_map_file:
        args.extend(['--sample-map-file', sample_map_file])

    if output != './blast_output':
        args.extend(['-o', output])

    if blast_type != 'blastn':
        args.extend(['--blast-type', blast_type])

    if evalue != 1e-5:
        args.extend(['--evalue', str(evalue)])

    if max_target_seqs != 10:
        args.extend(['--max-target-seqs', str(max_target_seqs)])

    if word_size != 11:
        args.extend(['--word-size', str(word_size)])

    if threads != 4:
        args.extend(['--threads', str(threads)])

    if input_suffix != '*.fa':
        args.extend(['--input-suffix', input_suffix])

    if target_db_type != 'nucl':
        args.extend(['--target-db-type', target_db_type])

    if min_identity != 70.0:
        args.extend(['--min-identity', str(min_identity)])

    if min_coverage != 50.0:
        args.extend(['--min-coverage', str(min_coverage)])

    if high_quality_evalue != 1e-10:
        args.extend(['--high-quality-evalue', str(high_quality_evalue)])

    if not auto_detect_samples:
        args.append('--no-auto-detect-samples')

    if sample_name_pattern != r'([^/]+?)(?:\.fa|\.fasta|\.fna)?$':
        args.extend(['--sample-name-pattern', sample_name_pattern])

    if makeblastdb_path != 'makeblastdb':
        args.extend(['--makeblastdb-path', makeblastdb_path])

    if blastn_path != 'blastn':
        args.extend(['--blastn-path', blastn_path])

    if blastp_path != 'blastp':
        args.extend(['--blastp-path', blastp_path])

    if blastx_path != 'blastx':
        args.extend(['--blastx-path', blastx_path])

    if tblastn_path != 'tblastn':
        args.extend(['--tblastn-path', tblastn_path])

    if tblastx_path != 'tblastx':
        args.extend(['--tblastx-path', tblastx_path])

    if alignment_output != 'both':
        args.extend(['--alignment-output', alignment_output])

    if alignment_width != 80:
        args.extend(['--alignment-width', str(alignment_width)])

    if alignment_min_identity != 0.0:
        args.extend(['--alignment-min-identity', str(alignment_min_identity)])

    if alignment_min_coverage != 0.0:
        args.extend(['--alignment-min-coverage', str(alignment_min_coverage)])

    if alignment_max_per_sample != 100:
        args.extend(['--alignment-max-per-sample', str(alignment_max_per_sample)])

    if html_theme != 'modern':
        args.extend(['--html-theme', html_theme])

    # 处理verbose (count)|Handle verbose (count)
    if verbose > 0:
        args.extend(['--verbose'] * verbose)

    if quiet:
        args.append('--quiet')

    if log_level != 'INFO':
        args.extend(['--log-level', log_level])

    if log_file:
        args.extend(['--log-file', log_file])

    if force:
        args.append('--force')

    if dry_run:
        args.append('--dry-run')

    # 保存并恢复sys.argv|Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        # 调用原始的main函数|Call original main function
        blast_main_func()
    except SystemExit as e:
        # 处理程序正常退出|Handle normal program exit
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
