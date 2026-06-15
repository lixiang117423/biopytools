"""候选基因RNA-seq转录验证命令|Candidate Gene RNA-seq Transcriptional Validation Command"""

import click
import sys
import os


def _lazy_import_main():
    """延迟加载主函数|Lazy load main function"""
    try:
        from ...gene_rnaseq_check.main import main as module_main
        return module_main
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
    return path


@click.command(
    short_help='候选基因RNA-seq转录验证|Candidate gene RNA-seq transcriptional validation',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120),
)
@click.option('-g', '--genome', required=True,
              callback=lambda ctx, param, value: _validate_path_exists(value) if value else None,
              help='基因组FASTA文件|Genome FASTA file')
@click.option('-a', '--annotation', required=True,
              callback=lambda ctx, param, value: _validate_path_exists(value) if value else None,
              help='GFF3注释文件|GFF3 annotation file')
@click.option('-e', '--gene-list', required=True,
              callback=lambda ctx, param, value: _validate_path_exists(value) if value else None,
              help='目标基因ID列表文件|Target gene ID list file')
@click.option('-r', '--reads-dir', required=True,
              callback=lambda ctx, param, value: _validate_path_exists(value) if value else None,
              help='RNA-seq reads目录|RNA-seq reads directory')
@click.option('-o', '--output-dir', default='./gene_rnaseq_check_output',
              help='输出目录|Output directory')
@click.option('-t', '--threads', default=12,
              help='线程数|Threads')
@click.option('--reads-pattern', default='',
              help='FASTQ命名模式|FASTQ naming pattern')
@click.option('--steps', default='all',
              help='执行步骤(逗号分隔或all)|Steps (comma-separated or all)')
@click.option('--sample-timeout', default=21600,
              help='单样本超时时间(秒)|Sample timeout in seconds')
@click.option('--flanking-window', default=500,
              help='上下游分析窗口(bp)|Flanking analysis window (bp)')
@click.option('--junction-tolerance', default=5,
              help='Junction容差(bp)|Junction tolerance (bp)')
@click.option('--strandness-confidence', default=70.0,
              help='链特异性判定置信度(%%)|Strandness confidence threshold (%%)')
@click.option('-v', '--verbose', is_flag=True,
              help='详细模式|Verbose mode')
@click.option('--quiet', is_flag=True,
              help='静默模式|Quiet mode')
@click.option('--force', is_flag=True,
              help='强制重新运行|Force re-run')
def gene_rnaseq_check(genome, annotation, gene_list, reads_dir, output_dir,
                      threads, reads_pattern, steps, sample_timeout,
                      flanking_window, junction_tolerance, strandness_confidence,
                      verbose, quiet, force):
    """候选基因RNA-seq转录验证|Candidate gene RNA-seq transcriptional validation

    示例|Examples: biopytools gene-rnaseq-check -g genome.fa -a anno.gff -e genes.txt -r ./reads/ -o out
    """
    module_main = _lazy_import_main()

    args = ['gene_rnaseq_check.py']
    args.extend(['-g', genome])
    args.extend(['-a', annotation])
    args.extend(['-e', gene_list])
    args.extend(['-r', reads_dir])
    if output_dir != './gene_rnaseq_check_output':
        args.extend(['-o', output_dir])
    if threads != 12:
        args.extend(['-t', str(threads)])
    if reads_pattern:
        args.extend(['--reads-pattern', reads_pattern])
    if steps != 'all':
        args.extend(['--steps', steps])
    if sample_timeout != 21600:
        args.extend(['--sample-timeout', str(sample_timeout)])
    if flanking_window != 500:
        args.extend(['--flanking-window', str(flanking_window)])
    if junction_tolerance != 5:
        args.extend(['--junction-tolerance', str(junction_tolerance)])
    if strandness_confidence != 70.0:
        args.extend(['--strandness-confidence', str(strandness_confidence)])
    if verbose:
        args.append('--verbose')
    if quiet:
        args.append('--quiet')
    if force:
        args.append('--force')

    original_argv = sys.argv
    sys.argv = args
    try:
        module_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
