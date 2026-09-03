"""全长转录本分析CLI包装器|Full-length RNA analysis CLI wrapper"""

import sys

import click


def _lazy_import_main():
    """延迟加载主函数|Lazy load main function"""
    try:
        from ...rna_iso.main import main as rna_iso_main
        return rna_iso_main
    except ImportError as e:
        error_msg = str(e)

        def error_func():
            click.echo(f"导入错误|Import error: {error_msg}", err=True)
            sys.exit(1)

        return error_func


@click.command(short_help='全长转录本分析(IsoSeq3+IsoQuant)|Full-length transcript analysis',
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('-i', '--reads',
              multiple=True,
              required=True,
              type=click.Path(exists=True),
              help='输入文件(可多个):subreads.bam/ccs.bam/fasta/fastq(±gz)|Input files (multiple allowed)')
@click.option('--data-type',
              type=click.Choice(['pacbio', 'ont']),
              default=None,
              help='reads文件(fasta/fastq)时必填;BAM输入自动嗅探|Required for fasta/fastq inputs; auto-sniffed for BAM')
@click.option('-g', '--reference',
              type=click.Path(exists=True),
              help='参考基因组FASTA(isoquant引擎必填)|Reference genome FASTA (required for isoquant engine)')
@click.option('--genedb',
              type=click.Path(exists=True),
              help='参考注释GTF/GFF(可选,提高精确度)|Reference annotation GTF/GFF (optional)')
@click.option('--engine',
              type=click.Choice(['isoquant', 'isoseq3', 'both']),
              default='isoquant',
              show_default=True,
              help='转录本重建引擎|Transcript engine')
@click.option('--primers',
              type=click.Path(exists=True),
              help='引物fasta(默认内置Clontech SMARTer)|Primer fasta (default: built-in Clontech SMARTer)')
@click.option('--min-passes',
              type=int,
              default=1,
              show_default=True,
              help='ccs最小pass数(Iso-Seq官方推荐1)|ccs min passes')
@click.option('-t', '--threads',
              type=int,
              default=12,
              show_default=True,
              help='线程数|Number of threads')
@click.option('-p', '--prefix',
              default='rna_sample',
              show_default=True,
              help='样本前缀|Sample prefix')
@click.option('-o', '--output-dir',
              default='./rna_iso_output',
              show_default=True,
              type=click.Path(),
              help='输出目录|Output directory')
def rna_iso(reads, data_type, reference, genedb, engine, primers,
            min_passes, threads, prefix, output_dir):
    """全长转录本分析流程|Full-length transcript analysis pipeline

    融合IsoSeq3(ccs+refine+cluster2)与IsoQuant(参考引导重建+定量),
    三种输入形态一条命令到全长转录本
    |Fuses IsoSeq3 and IsoQuant; one command from PacBio subreads/CCS or ONT
    reads to full-length transcripts.

    示例|Examples: biopytools rna-iso -i sample.fq.gz --data-type ont -g genome.fa -o output/
    """

    rna_iso_main = _lazy_import_main()

    args = ['rna_iso.py']

    # 必需参数|Required parameters
    args.extend(['--reads'] + list(reads))

    # 可选参数(仅显式传入时透传)|Optional (forward only when set)
    if data_type is not None:
        args.extend(['--data-type', data_type])
    if reference is not None:
        args.extend(['--reference', reference])
    if genedb is not None:
        args.extend(['--genedb', genedb])
    if engine != 'isoquant':
        args.extend(['--engine', engine])
    if primers is not None:
        args.extend(['--primers', primers])
    if min_passes != 1:
        args.extend(['--min-passes', str(min_passes)])
    if threads != 12:
        args.extend(['--threads', str(threads)])
    if prefix != 'rna_sample':
        args.extend(['--prefix', prefix])
    if output_dir != './rna_iso_output':
        args.extend(['--output-dir', output_dir])

    original_argv = sys.argv
    try:
        sys.argv = args
        rna_iso_main()
    except SystemExit as e:
        if e.code != 0:
            sys.exit(e.code)
    finally:
        sys.argv = original_argv
