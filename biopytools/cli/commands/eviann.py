"""
EviAnn基因组注释CLI包装器|EviAnn Genome Annotation CLI Wrapper
"""

import click
import os
import sys


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


def _validate_rnaseq_data(value):
    """校验逗号分隔的转录组数据路径|Validate comma-separated data paths"""
    if not value or _is_help_request():
        return value
    for entry in value.split(','):
        entry = entry.strip()
        if entry and not os.path.exists(os.path.expanduser(entry)):
            raise click.BadParameter(
                f"路径不存在|Path does not exist: {entry}")
    return value


@click.command(short_help='EviAnn基因组注释|EviAnn Genome Annotation',
               context_settings=dict(help_option_names=['-h', '--help'],
                                     max_content_width=120))
@click.option('-g', '--genome',
              required=True,
              type=click.Path(exists=True),
              help='基因组FASTA文件|Genome FASTA file (required)')
@click.option('-o', '--output-dir',
              required=True,
              help='输出目录|Output directory (required)')
@click.option('--rnaseq-data',
              default=None,
              callback=lambda ctx, param, value: _validate_rnaseq_data(value),
              help='转录组数据文件或目录(逗号分隔多个),自动识别二代/三代|'
                   'RNA-seq file(s) or dir(s), comma-separated')
@click.option('--sample-sheet',
              default=None,
              type=click.Path(exists=True),
              help='样本清单TSV|Sample sheet TSV')
@click.option('-r', '--rnaseq',
              default=None,
              type=click.Path(exists=True),
              help='EviAnn原生-r描述文件(透传)|EviAnn native -r file')
@click.option('-e', '--transcripts',
              default=None,
              type=click.Path(exists=True),
              help='近缘物种转录本FASTA|Transcripts FASTA')
@click.option('-p', '--proteins',
              default=None,
              type=click.Path(exists=True),
              help='近缘物种蛋白质FASTA|Proteins FASTA')
@click.option('-s', '--uniprot',
              default=None,
              type=click.Path(exists=True),
              help='UniProt-SwissProt FASTA|UniProt-SwissProt FASTA')
@click.option('-t', '--threads',
              default=12,
              type=int,
              show_default=True,
              help='线程数|Number of threads')
@click.option('-m', '--max-intron',
              default=None,
              type=int,
              help='最大内含子长度|Maximum intron length (default: auto)')
@click.option('-d', '--ploidy',
              default=2,
              type=int,
              show_default=True,
              help='基因组倍性|Genome ploidy')
@click.option('-c', '--cds-gff',
              default=None,
              type=click.Path(exists=True),
              help='含现有CDS的GFF|GFF with existing CDS')
@click.option('--lncrna-tpm',
              default=1.0,
              type=float,
              show_default=True,
              help='lncRNA最小TPM|Minimum TPM for lncRNA')
@click.option('--min-prot',
              default=None,
              type=int,
              help='无同源证据时ab initio ORF最小蛋白长度(aa)|'
                   'Min protein length for ab initio ORF')
@click.option('--partial',
              is_flag=True,
              default=False,
              help='包含部分CDS|Include partial CDS')
@click.option('-f', '--functional',
              is_flag=True,
              default=False,
              help='执行功能注释|Perform functional annotation')
@click.option('--mito-contigs',
              default=None,
              type=click.Path(exists=True),
              help='线粒体contig列表文件|File with mitochondrial contigs')
@click.option('--extra-gff',
              default=None,
              type=click.Path(exists=True),
              help='额外GFF特征|Extra features from external GFF')
@click.option('--debug',
              is_flag=True,
              default=False,
              help='保留中间文件|Keep intermediate files')
@click.option('--verbose',
              is_flag=True,
              default=False,
              help='详细输出|Verbose output')
def eviann(genome, output_dir, rnaseq_data, sample_sheet, rnaseq,
           transcripts, proteins, uniprot, threads, max_intron, ploidy,
           cds_gff, lncrna_tpm, min_prot, partial, functional, mito_contigs,
           extra_gff, debug, verbose):
    """
    EviAnn基因组注释流程|EviAnn Genome Annotation Pipeline

    基于RNA-seq和/或蛋白质比对的证据驱动真核基因组注释,自动识别
    二代(双端配对)/三代数据并生成EviAnn输入描述文件
    |Evidence-based eukaryotic genome annotation; auto-classifies
    short-read pairs and long-read data into EviAnn input

    示例|Example: biopytools eviann -g genome.fa --rnaseq-data ./rna_data/ -p proteins.fa -t 12 -o out/
    """

    # 延迟加载|Lazy loading
    eviann_main = _lazy_import_eviann_main()

    # 构建参数列表|Build argument list
    args = ['eviann.py']
    args.extend(['-g', genome, '-o', output_dir])

    # 数据输入|Data inputs
    if rnaseq_data:
        args.extend(['--rnaseq-data', rnaseq_data])
    if sample_sheet:
        args.extend(['--sample-sheet', sample_sheet])
    if rnaseq:
        args.extend(['-r', rnaseq])
    if transcripts:
        args.extend(['-e', transcripts])
    if proteins:
        args.extend(['-p', proteins])

    # 可选参数(默认值显式透传)|Optional args (explicit passthrough)
    if uniprot:
        args.extend(['-s', uniprot])
    args.extend(['-t', str(threads)])
    if max_intron:
        args.extend(['-m', str(max_intron)])
    args.extend(['-d', str(ploidy)])
    if cds_gff:
        args.extend(['-c', cds_gff])
    args.extend(['--lncrna-tpm', str(lncrna_tpm)])
    if min_prot:
        args.extend(['--min-prot', str(min_prot)])
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
