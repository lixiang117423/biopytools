"""
ps-gene-anno: BRAKER 后效应子查漏补缺|Post-BRAKER effector gap-filling
"""

import click
import sys
import os


def _lazy_import_main():
    """延迟加载主函数|Lazy load main"""
    try:
        from ...ps_gene_anno import main as psga_module
        return psga_module.main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    return any(a in {'-h', '--help'} for a in sys.argv)


def _validate_file(path):
    if not _is_help_request() and path and not os.path.exists(path):
        raise click.BadParameter(f"文件不存在|File not found: {path}")
    return path


@click.command(
    short_help='BRAKER后效应子查漏补缺(miniprot证据补回多拷贝)|Post-BRAKER effector gap-filling',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('-g', '--genome', required=True,
              callback=lambda c, p, v: _validate_file(v) if v else None,
              help='未mask原始基因组|Unmasked raw genome')
@click.option('-b', '--braker-gff3', required=True,
              callback=lambda c, p, v: _validate_file(v) if v else None,
              help='BRAKER输出GFF3|BRAKER output GFF3')
@click.option('-p', '--prot-seq', required=True,
              callback=lambda c, p, v: _validate_file(v) if v else None,
              help='近缘蛋白|Protein evidence')
@click.option('-o', '--output-dir', required=True, help='输出目录|Output directory')
@click.option('--rnaseq-bam', help='RNA-seq BAM(逗号分隔)|RNA-seq BAMs')
@click.option('--isoseq-bam', help='三代BAM|Long-read BAM')
@click.option('--repeat-out', help='RepeatMasker .out(真TE区排除)|RepeatMasker out')
@click.option('--prefix', help='输出前缀(默认genome stem)|Output prefix')
@click.option('-t', '--threads', type=int, default=12, show_default=True,
              help='线程数|Threads')
@click.option('--split-min-copy-coverage', type=float, default=80, show_default=True,
              help='保守合并判据:完整拷贝覆盖率%|Conservative split copy coverage')
@click.option('--no-split', is_flag=True, help='关闭合并拆分|Disable merged-gene split')
def ps_gene_anno(genome, braker_gff3, prot_seq, output_dir,
                 rnaseq_bam, isoseq_bam, repeat_out, prefix,
                 threads, split_min_copy_coverage, no_split):
    """
    BRAKER后效应子查漏补缺|Post-BRAKER effector gap-filling

    示例|Example: biopytools ps-gene-anno -g genome.fa -b braker.gff3 -p prot.fa -o out/
    """
    psga_main = _lazy_import_main()
    args = ['ps_gene_anno.py', '-g', genome, '-b', braker_gff3,
            '-p', prot_seq, '-o', output_dir]
    if rnaseq_bam:
        args.extend(['--rnaseq-bam', rnaseq_bam])
    if isoseq_bam:
        args.extend(['--isoseq-bam', isoseq_bam])
    if repeat_out:
        args.extend(['--repeat-out', repeat_out])
    if prefix:
        args.extend(['--prefix', prefix])
    if threads != 12:
        args.extend(['-t', str(threads)])
    if split_min_copy_coverage != 80:
        args.extend(['--split-min-copy-coverage', str(split_min_copy_coverage)])
    if no_split:
        args.append('--no-split')

    original_argv = sys.argv
    sys.argv = args
    try:
        psga_main()
    except SystemExit as e:
        sys.exit(e.code)
    finally:
        sys.argv = original_argv
