"""|Sequence Extraction Tool Command
--help
"""

import click
import sys
import os


def _lazy_import_sequence_main():
    """main|Lazy load sequence extraction main function"""
    try:
        from ...parse_seq.main import main as sequence_main
        return sequence_main
    except ImportError as e:
        click.echo(f"Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_file_exists(file_path):
    """| Validate file existence (only in non-help mode)"""
    if not _is_help_request() and not os.path.exists(file_path):
        raise click.BadParameter(f"File does not exist: {file_path}")
    return file_path


@click.command(short_help='序列提取工具|Sequence extraction tool',
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--sequence-file', '-s',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='输入序列文件(FASTA格式)|Input sequence file (FASTA format)')
@click.option('--regions-file', '-r',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='区域文件(BED格式:染色体 起始 终止)|Regions file (BED format: chromosome start end)')
@click.option('--output-file', '-o',
              required=True,
              type=click.Path(),
              help='输出序列文件|Output sequence file')
@click.option('--type', '--sequence-type',
              type=click.Choice(['dna', 'protein']),
              default='dna',
              help='序列类型|Sequence type (default: dna)')
@click.option('--threads', '-t',
              type=int,
              default=12,
              help='线程数|Number of threads (default: 12)')
@click.option('--merge-output',
              is_flag=True,
              default=True,
              help='合并输出到一个文件|Merge output to one file')
@click.option('--separate-output',
              is_flag=True,
              help='分别输出到多个文件|Output to separate files')
@click.option('--no-headers',
              is_flag=True,
              help='不包含区域信息在序列名中|Do not include region info in sequence names')
@click.option('--reverse-complement',
              is_flag=True,
              help='反向互补DNA序列|Reverse complement DNA sequences')
@click.option('--translate',
              is_flag=True,
              help='将DNA翻译为蛋白质|Translate DNA to protein')
@click.option('--line-width',
              type=int,
              default=80,
              help='FASTA序列每行字符数|Characters per line in FASTA sequence (default: 80)')
@click.option('--verbose', '-v',
              is_flag=True,
              default=True,
              help='详细输出|Verbose output')
@click.option('--quiet', '-q',
              is_flag=True,
              help='安静模式|Quiet mode')
@click.option('--samtools-path',
              default='samtools',
              help='samtools程序路径|samtools program path')
def parse_seq(sequence_file, regions_file, output_file, type, threads, merge_output,
                      separate_output, no_headers, reverse_complement, translate, line_width,
                      verbose, quiet, samtools_path):
    """
    序列提取工具：从FASTA文件中提取指定区域序列|Sequence extraction tool: Extract specified regions from FASTA files

    示例|Examples: biopytools parse-seq -s genome.fasta -r regions.bed -o extracted.fasta
    """
    
    # Lazy loading: import only when actually called
    sequence_main = _lazy_import_sequence_main()
    
    # main|Build argument list for original main function
    args = ['parse_seq.py']
    
    # Required parameters
    args.extend(['-s', sequence_file])
    args.extend(['-r', regions_file])
    args.extend(['-o', output_file])
    
    #Optional parameters (add only when non-default)
    if type != 'dna':
        args.extend(['--type', type])
    
    if threads != 12:
        args.extend(['-t', str(threads)])
    
    if line_width != 80:
        args.extend(['--line-width', str(line_width)])
    
    if samtools_path != 'samtools':
        args.extend(['--samtools-path', samtools_path])
    
    #Handle mutually exclusive parameters
    if separate_output:
        args.append('--separate-output')
    elif not merge_output:
        # merge_outputseparate_output
        pass
    
    if quiet:
        args.append('--quiet')
    elif not verbose:
        # verbosequiet
        pass
    
    # Boolean options
    if no_headers:
        args.append('--no-headers')
    
    if reverse_complement:
        args.append('--reverse-complement')
    
    if translate:
        args.append('--translate')
    
    # sys.argv|Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args
    
    try:
        # main|Call original main function
        sequence_main()
    except SystemExit as e:
        # Handle normal program exit
        if e.code != 0:
            sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n操作被用户中断|Operation interrupted by user", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"序列提取失败|Sequence extraction failed: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv