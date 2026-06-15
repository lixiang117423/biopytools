"""
KMC k-mer分析命令|KMC K-mer Analysis Command
"""

import click
import sys
import os


def _lazy_import_kmc_main():
    """延迟加载kmc主函数|Lazy load kmc main function"""
    try:
        from ...kmc.main import main as kmc_main
        return kmc_main
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


# 通用参数选项|Common option decorators
def input_files_option(f):
    """输入文件选项装饰器|Input files option decorator"""
    return click.option('--input', '-i',
                       multiple=True,
                       required=True,
                       help='输入文件(FASTQ/FASTA)，可多次使用|Input files (FASTQ/FASTA), can be used multiple times')(f)


def sample_names_option(f):
    """样本名称选项装饰器|Sample names option decorator"""
    return click.option('--sample-names', '-n',
                       multiple=True,
                       help='样本名称，可多次使用(默认使用文件名)|Sample names, can be used multiple times (default: use filename)')(f)


def kmer_size_option(f):
    """k-mer大小选项装饰器|k-mer size option decorator"""
    return click.option('--kmer-size', '-k',
                       type=int,
                       default=21,
                       show_default=True,
                       help='k-mer大小|k-mer size')(f)


def min_count_option(f):
    """最小计数选项装饰器|Minimum count option decorator"""
    return click.option('--min-count',
                       type=int,
                       default=2,
                       show_default=True,
                       help='最小计数阈值|Minimum count threshold')(f)


def max_count_option(f):
    """最大计数选项装饰器|Maximum count option decorator"""
    return click.option('--max-count',
                       type=int,
                       help='最大计数阈值|Maximum count threshold')(f)


def output_dir_option(f):
    """输出目录选项装饰器|Output directory option decorator"""
    return click.option('--output-dir', '-o',
                       default='./kmc_output',
                       show_default=True,
                       type=click.Path(),
                       help='输出目录|Output directory')(f)


def kmc_path_option(f):
    """KMC路径选项装饰器|KMC path option decorator"""
    return click.option('--kmc-path',
                       default='~/miniforge3/envs/kmc_v.3.2.4/bin',
                       show_default=True,
                       help='KMC软件路径|KMC software path')(f)


def tmp_dir_option(f):
    """临时目录选项装饰器|Temporary directory option decorator"""
    return click.option('--tmp-dir',
                       default='./kmc_tmp',
                       show_default=True,
                       type=click.Path(),
                       help='临时文件目录|Temporary directory')(f)


def threads_option(f):
    """线程数选项装饰器|Threads option decorator"""
    return click.option('--threads', '-t',
                       type=int,
                       default=12,
                       show_default=True,
                       help='线程数|Number of threads')(f)


def memory_limit_option(f):
    """内存限制选项装饰器|Memory limit option decorator"""
    return click.option('--memory-limit',
                       help='内存限制(如: 12G)|Memory limit (e.g., 12G)')(f)


@click.group(
    short_help='KMC k-mer统计和分析工具|KMC k-mer counting and analysis tool',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
def kmc():
    """
    KMC k-mer分析工具|KMC K-mer Analysis Tool

    使用KMC进行k-mer统计、丰度矩阵构建和查询|Use KMC for k-mer counting, abundance matrix building and query

    示例|Example: biopytools kmc count -i sample1.fq -i sample2.fq -k 21 -o kmc_output
    """
    pass


@kmc.command(short_help='统计k-mer|Count k-mers')
@click.option('--input-dir', '-d',
              help='输入目录(自动识别双末端测序)|Input directory (auto-detect paired-end)')
@click.option('--input', '-i',
              multiple=True,
              help='输入文件(FASTQ/FASTA)，可多次使用|Input files (FASTQ/FASTA), can be used multiple times')
@click.option('--sample-names', '-n',
              multiple=True,
              help='样本名称，可多次使用(默认使用文件名)|Sample names, can be used multiple times (default: use filename)')
@click.option('--read1-suffix',
              default='_1.clean.fq.gz',
              show_default=True,
              help='Read1文件后缀|Read1 file suffix')
@click.option('--read2-suffix',
              default='_2.clean.fq.gz',
              show_default=True,
              help='Read2文件后缀|Read2 file suffix')
@click.option('--single-end',
              is_flag=True,
              help='单末端测序模式|Single-end sequencing mode')
@kmer_size_option
@min_count_option
@max_count_option
@output_dir_option
@tmp_dir_option
@kmc_path_option
@threads_option
@memory_limit_option
def count(input_dir, input, sample_names, read1_suffix, read2_suffix, single_end,
          kmer_size, min_count, max_count, output_dir, tmp_dir, kmc_path, threads, memory_limit):
    """
    统计k-mer|Count k-mers

    为每个样本建立KMC数据库，统计k-mer丰度|Build KMC database for each sample and count k-mer abundance

    示例|Example:
      # 目录模式(自动识别双末端)|Directory mode (auto-detect paired-end):
      biopytools kmc count -d raw_data/ -o kmc_output

      # 文件列表模式|File list mode:
      biopytools kmc count -i sample1.fq -i sample2.fq -k 21 -o output

      # 单末端模式|Single-end mode:
      biopytools kmc count -d raw_data/ --single-end -o output
    """
    kmc_main = _lazy_import_kmc_main()

    args = _build_args('count', {
        'input_dir': input_dir,
        'input': input,
        'sample_names': sample_names,
        'read1_suffix': read1_suffix,
        'read2_suffix': read2_suffix,
        'single_end': single_end,
        'kmer_size': kmer_size,
        'min_count': min_count,
        'max_count': max_count,
        'output_dir': output_dir,
        'tmp_dir': tmp_dir,
        'kmc_path': kmc_path,
        'threads': threads,
        'memory_limit': memory_limit
    })

    _execute_main(kmc_main, args)


@kmc.command(short_help='构建丰度矩阵|Build abundance matrix')
@click.option('--input-dir', '-i',
              help='包含kmc_databases的目录(即count步骤的-o参数)|Directory containing kmc_databases (i.e., -o from count step)')
@output_dir_option
@kmc_path_option
@threads_option
@click.option('--max-memory',
              type=int,
              default=500,
              show_default=True,
              help='最大内存使用量(GB)|Maximum memory usage (GB)')
@click.option('--matrix-format',
              type=click.Choice(['hdf5', 'tsv', 'sqlite']),
              default='hdf5',
              show_default=True,
              help='矩阵存储格式|Matrix storage format')
@click.option('--dense-matrix',
              is_flag=True,
              help='使用密集矩阵(默认稀疏)|Use dense matrix (default: sparse)')
@click.option('--no-export',
              is_flag=True,
              help='不导出TSV文件(默认自动导出)|Do not export TSV files (auto export by default)')
@click.option('--keep-dump',
              is_flag=True,
              help='保留dump文件到dump_files目录(默认保留)|Keep dump files in dump_files directory (default: True)')
@click.option('--no-keep-dump',
              is_flag=True,
              help='不保留dump文件(节省空间)|Do not keep dump files (save space)')
def matrix(input_dir, output_dir, kmc_path, threads, max_memory, matrix_format, dense_matrix, no_export, keep_dump, no_keep_dump):
    """
    构建k-mer丰度矩阵|Build k-mer abundance matrix

    基于已存在的KMC数据库自动构建跨样本丰度矩阵|Auto-build cross-sample abundance matrix from existing KMC databases

    注意：需要先运行count命令建立数据库|Note: Run count command first to build databases

    注意：
      - 如果指定-i/--input-dir，使用该目录作为输入目录|If -i/--input-dir specified, use it as input directory
      - 如果未指定-i，使用-o/--output-dir（即扫描output_dir/kmc_databases/）|If -i not specified, use -o (scan output_dir/kmc_databases/)

    示例|Example:
      # 方式1：使用-i指定输入目录|Method 1: Use -i to specify input directory
      biopytools kmc matrix -i kmc_database -o kmc_database

      # 方式2：不指定-i，使用-o作为输入和输出目录|Method 2: Don't specify -i, use -o for both input and output
      biopytools kmc matrix -o kmc_database
    """
    kmc_main = _lazy_import_kmc_main()

    args = _build_args('matrix', {
        'input_dir': input_dir,
        'output_dir': output_dir,
        'kmc_path': kmc_path,
        'threads': threads,
        'max_memory': max_memory,
        'matrix_format': matrix_format,
        'dense_matrix': dense_matrix,
        'no_export': no_export,
        'keep_dump': not no_keep_dump  # 默认True，--no-keep-dump时为False
    })

    _execute_main(kmc_main, args)


@kmc.command(short_help='查询k-mer|Query k-mer')
@click.option('--input-fasta', '-f',
              required=True,
              help='查询的FASTA文件(支持百万级批量查询)|FASTA file for query (supports millions)')
@click.option('--output-file', '-o',
              help='输出文件(可选，默认输出到屏幕)|Output file (optional, default to screen)')
@click.option('--db-dir', '--output-dir', 'output_dir',
              default='./kmc_output',
              show_default=True,
              help='KMC数据库目录(包含abundance_matrix.h5等文件)|KMC database directory (contains abundance_matrix.h5, etc.)')
@click.option('--threads', '-t',
              type=int,
              default=12,
              show_default=True,
              help='线程数|Number of threads')
def query(input_fasta, output_file, output_dir, threads):
    """
    查询k-mer丰度|Query k-mer abundance

    从FASTA文件中查询所有序列在各样本中的丰度|Query all sequences in FASTA across samples

    示例|Example: biopytools kmc query -f query.fasta --db-dir output -o results.tsv
    """
    kmc_main = _lazy_import_kmc_main()

    args = _build_args('query', {
        'query_fasta': input_fasta,
        'output_file': output_file,
        'output_dir': output_dir,
        'threads': threads
    })

    _execute_main(kmc_main, args)


@kmc.command(short_help='添加新样本|Add new samples')
@click.option('--input-dir', '-d',
              help='输入目录(自动识别双末端测序)|Input directory (auto-detect paired-end)')
@click.option('--input', '-i',
              multiple=True,
              help='新样本文件(可选，不指定则仅更新矩阵)|New sample files (optional, if not specified, only update matrix)')
@click.option('--sample-names', '-n',
              multiple=True,
              help='新样本名称(可选)|New sample names (optional)')
@click.option('--read1-suffix',
              default='_1.clean.fq.gz',
              show_default=True,
              help='Read1文件后缀|Read1 file suffix')
@click.option('--read2-suffix',
              default='_2.clean.fq.gz',
              show_default=True,
              help='Read2文件后缀|Read2 file suffix')
@click.option('--single-end',
              is_flag=True,
              help='单末端测序模式|Single-end sequencing mode')
@output_dir_option
@tmp_dir_option
@kmc_path_option
@threads_option
@memory_limit_option
@click.option('--max-memory',
              type=int,
              default=500,
              show_default=True,
              help='最大内存使用量(GB)|Maximum memory usage (GB)')
@click.option('--matrix-format',
              type=click.Choice(['hdf5', 'tsv', 'sqlite']),
              default='hdf5',
              show_default=True,
              help='矩阵存储格式|Matrix storage format')
@click.option('--dense-matrix',
              is_flag=True,
              help='使用密集矩阵(默认稀疏)|Use dense matrix (default: sparse)')
def add(input_dir, input, sample_names, read1_suffix, read2_suffix, single_end,
        output_dir, tmp_dir, kmc_path, threads, memory_limit, max_memory, matrix_format, dense_matrix):
    """
    添加新样本到矩阵|Add new samples to matrix

    如果指定了新样本文件，会先为它们建立数据库，然后更新矩阵|If new sample files specified, build databases for them first, then update matrix
    如果未指定，仅使用现有数据库更新矩阵|If not specified, only update matrix using existing databases

    示例|Example: biopytools kmc add -d new_samples/ -o output
    """
    kmc_main = _lazy_import_kmc_main()

    args = _build_args('add', {
        'input_dir': input_dir,
        'input': input,
        'sample_names': sample_names,
        'read1_suffix': read1_suffix,
        'read2_suffix': read2_suffix,
        'single_end': single_end,
        'output_dir': output_dir,
        'tmp_dir': tmp_dir,
        'kmc_path': kmc_path,
        'threads': threads,
        'memory_limit': memory_limit,
        'max_memory': max_memory,
        'matrix_format': matrix_format,
        'dense_matrix': dense_matrix
    })

    _execute_main(kmc_main, args)


@kmc.command(short_help='导出矩阵为TSV|Export matrix to TSV')
@click.option('--input-dir', '-i',
              help='包含矩阵文件的目录(即count/matrix步骤的-o参数)|Directory containing matrix files (i.e., -o from count/matrix step)')
@click.option('--output-file', '-o',
              default='abundance_matrix.tsv',
              show_default=True,
              help='输出TSV文件路径|Output TSV file path')
@click.option('--format',
              type=click.Choice(['full', 'sparse']),
              default='sparse',
              show_default=True,
              help='输出格式|Output format:\n'
                   'full: 完整矩阵(所有k-mer x 所有样本)|full matrix (all k-mers x all samples)\n'
                   'sparse: 稀疏格式(只包含非零值)|sparse format (non-zero values only)')
@click.option('--min-abundance',
              type=int,
              default=1,
              show_default=True,
              help='最小丰度阈值|Minimum abundance threshold')
def export(input_dir, output_file, format, min_abundance):
    """
    导出HDF5矩阵为TSV格式|Export HDF5 matrix to TSV format

    将构建好的k-mer丰度矩阵从HDF5格式导出为TSV文本格式|Export k-mer abundance matrix from HDF5 format to TSV text format

    注意：需要先运行matrix命令构建矩阵|Note: Run matrix command first to build matrix

    示例|Example:
      # 导出为稀疏格式（推荐，适合大规模数据）
      biopytools kmc export -i output -o output/abundance_sparse.tsv

      # 导出为完整矩阵（适合小规模数据）
      biopytools kmc export -i output -o output/abundance_full.tsv --format full

      # 只导出丰度>=5的k-mer
      biopytools kmc export -i output -o output/abundance_filtered.tsv --min-abundance 5
    """
    kmc_main = _lazy_import_kmc_main()

    args = _build_args('export', {
        'input_dir': input_dir,
        'output_file': output_file,
        'format': format,
        'min_abundance': min_abundance
    })

    _execute_main(kmc_main, args)


def _build_args(mode, params):
    """构建参数列表|Build argument list"""
    args = ['kmc.py', '-m', mode]

    # 输入目录或文件|Input directory or files
    if 'input_dir' in params and params['input_dir']:
        args.extend(['--input-dir', params['input_dir']])

    if 'input' in params and params['input']:
        for input_file in params['input']:
            args.extend(['-i', input_file])

    # 样本名称|Sample names
    if 'sample_names' in params and params['sample_names']:
        for sample_name in params['sample_names']:
            args.extend(['-n', sample_name])

    # 双末端测序参数|Paired-end parameters
    if 'read1_suffix' in params and params['read1_suffix'] != '_1.clean.fq.gz':
        args.extend(['--read1-suffix', params['read1_suffix']])

    if 'read2_suffix' in params and params['read2_suffix'] != '_2.clean.fq.gz':
        args.extend(['--read2-suffix', params['read2_suffix']])

    if 'single_end' in params and params['single_end']:
        args.append('--single-end')

    # k-mer大小|k-mer size
    if 'kmer_size' in params and params['kmer_size'] != 21:
        args.extend(['-k', str(params['kmer_size'])])

    # 最小计数|Minimum count
    if 'min_count' in params and params['min_count'] != 2:
        args.extend(['--min-count', str(params['min_count'])])

    # 最大计数|Maximum count
    if 'max_count' in params and params['max_count'] is not None:
        args.extend(['--max-count', str(params['max_count'])])

    # KMC路径|KMC path
    if 'kmc_path' in params and params['kmc_path'] != '~/miniforge3/envs/kmc_v.3.2.4/bin':
        args.extend(['--kmc-path', params['kmc_path']])

    # 输出目录/数据库目录|Output directory / Database directory
    if 'output_dir' in params and params['output_dir'] != './kmc_output':
        args.extend(['--output-dir', params['output_dir']])

    # 临时目录|Temporary directory
    if 'tmp_dir' in params and params['tmp_dir'] != './kmc_tmp':
        args.extend(['--tmp-dir', params['tmp_dir']])

    # 线程数|Threads
    if 'threads' in params and params['threads'] != 12:
        args.extend(['-t', str(params['threads'])])

    # 内存限制|Memory limit
    if 'memory_limit' in params and params['memory_limit']:
        args.extend(['--memory-limit', params['memory_limit']])

    # 最大内存|Maximum memory
    if 'max_memory' in params and params['max_memory'] != 500:
        args.extend(['--max-memory', str(params['max_memory'])])

    # 矩阵格式|Matrix format
    if 'matrix_format' in params and params['matrix_format'] != 'hdf5':
        args.extend(['--matrix-format', params['matrix_format']])

    # 密集矩阵|Dense matrix
    if 'dense_matrix' in params and params['dense_matrix']:
        args.append('--dense-matrix')

    # 不导出TSV|No export TSV
    if 'no_export' in params and params['no_export']:
        args.append('--no-export')

    # 保留dump文件|Keep dump files
    if 'keep_dump' in params and not params['keep_dump']:
        args.append('--no-keep-dump')

    # FASTA批量查询|FASTA batch query
    if 'query_fasta' in params and params['query_fasta']:
        args.extend(['-f', params['query_fasta']])

    # 输出文件|Output file
    if 'output_file' in params and params['output_file']:
        args.extend(['--output-file', params['output_file']])

    # Export格式|Export format
    if 'format' in params and params['format'] != 'sparse':
        args.extend(['--format', params['format']])

    # 最小丰度|Minimum abundance
    if 'min_abundance' in params and params['min_abundance'] != 1:
        args.extend(['--min-abundance', str(params['min_abundance'])])

    return args


def _execute_main(kmc_main, args):
    """执行主程序|Execute main program"""
    original_argv = sys.argv
    sys.argv = args

    try:
        kmc_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
