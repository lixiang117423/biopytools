"""
EviAnn基因组注释CLI包装器|EviAnn Genome Annotation CLI Wrapper
"""

import click
import sys
import os


def _lazy_import_eviann_main():
    """延迟加载eviann主函数|Lazy load eviann main function"""
    try:
        from ...eviann.main import main as eviann_main
        return eviann_main
    except ImportError as e:
        def error_func():
            click.echo(f"导入错误|Import Error: {e}", err=True)
            sys.exit(1)
        return error_func


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_path_exists(path):
    """验证路径存在|Validate path exists"""
    if not _is_help_request() and path and not os.path.exists(path):
        raise click.BadParameter(f"路径不存在|Path does not exist: {path}")
    return path


@click.command(short_help='EviAnn基因组注释|EviAnn Genome Annotation',
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('-g', '--genome',
              required=True,
              type=click.Path(exists=True),
              help='基因组FASTA文件|Genome FASTA file (required)')
@click.option('--short-reads',
              type=click.Path(exists=True),
              help='二代转录组数据（文件或目录）|Short-read RNA-seq data (file or directory)')
@click.option('--long-reads',
              type=click.Path(exists=True),
              help='三代转录组数据（文件或目录）|Long-read RNA-seq data (file or directory)')
@click.option('-e', '--transcripts',
              type=click.Path(exists=True),
              help='转录本FASTA文件|Transcripts FASTA file')
@click.option('-p', '--proteins',
              type=click.Path(exists=True),
              help='蛋白质FASTA文件|Proteins FASTA file')
@click.option('-s', '--uniprot',
              type=click.Path(exists=True),
              help='UniProt-SwissProt FASTA|UniProt-SwissProt FASTA')
@click.option('-t', '--threads',
              default=12,
              type=int,
              show_default=True,
              help='线程数|Number of threads')
@click.option('-m', '--max-intron',
              type=int,
              help='最大内含子长度|Maximum intron length (default: auto)')
@click.option('-d', '--ploidy',
              default=2,
              type=int,
              show_default=True,
              help='基因组倍性|Genome ploidy')
@click.option('-c', '--cds-gff',
              type=click.Path(exists=True),
              help='现有CDS的GFF文件|GFF file with existing CDS')
@click.option('--lncrna-tpm',
              default=1.0,
              type=float,
              show_default=True,
              help='lncRNA最小TPM|Minimum TPM for lncRNA')
@click.option('--partial',
              is_flag=True,
              default=False,
              help='包含部分CDS|Include partial CDS')
@click.option('--functional',
              is_flag=True,
              default=False,
              help='执行功能注释|Perform functional annotation')
@click.option('--mito-contigs',
              type=click.Path(exists=True),
              help='线粒体contig列表|File with mitochondrial contigs')
@click.option('--extra-gff',
              type=click.Path(exists=True),
              help='额外的GFF特征|Extra features from GFF')
@click.option('--debug',
              is_flag=True,
              default=False,
              help='调试模式|Debug mode')
@click.option('--verbose',
              is_flag=True,
              default=False,
              help='详细输出|Verbose output')
def eviann(genome, short_reads, long_reads, transcripts, proteins, uniprot, threads, max_intron,
           ploidy, cds_gff, lncrna_tpm, partial, functional, mito_contigs, extra_gff,
           debug, verbose):
    """
    EviAnn基因组注释流程|EviAnn Genome Annotation Pipeline

    基于RNA-seq和/或蛋白质比对进行真核生物基因组注释
    Evidence-based eukaryotic genome annotation using RNA-seq and/or protein alignments

    示例|Example: biopytools eviann -g genome.fa --long-reads longreads.fq.gz -p proteins.fa -t 12
    """

    # 延迟加载|Lazy loading
    eviann_main = _lazy_import_eviann_main()

    # 构建参数列表|Build argument list
    args = ['eviann.py']

    # 必需参数|Required parameters
    args.extend(['-g', genome])

    # 可选参数|Optional parameters
    if short_reads:
        args.extend(['--short-reads', short_reads])

    if long_reads:
        args.extend(['--long-reads', long_reads])

    if transcripts:
        args.extend(['-e', transcripts])

    if proteins:
        args.extend(['-p', proteins])

    if uniprot:
        args.extend(['-s', uniprot])

    if threads != 1:
        args.extend(['-t', str(threads)])

    if max_intron:
        args.extend(['-m', str(max_intron)])

    if ploidy != 2:
        args.extend(['-d', str(ploidy)])

    if cds_gff:
        args.extend(['-c', cds_gff])

    if lncrna_tpm != 1.0:
        args.extend(['--lncrna-tpm', str(lncrna_tpm)])

    if partial:
        args.append('--partial')

    if functional:
        args.append('--functional')

    if mito_contigs:
        args.extend(['--mito-contigs', mito_contigs])

    if extra_gff:
        args.extend(['--extra-gff', extra_gff])

    if debug:
        args.append('--debug')

    if verbose:
        args.append('--verbose')

    # 保存并恢复sys.argv|Save and restore sys.argv
    original_argv = sys.argv

    try:
        # 调用主函数|Call main function
        sys.argv = args
        eviann_main()

    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
