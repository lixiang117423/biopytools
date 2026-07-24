"""rnaseq2vcf CLI 包装(lazy)|Click lazy wrapper (VCF-only, no annotation)"""

import os
import sys
import click


def _lazy_import_main():
    """延迟加载主函数|Lazy load main"""
    try:
        from ...rnaseq2vcf.main import main as rnaseq2vcf_main
        return rnaseq2vcf_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    return any(a in ('-h', '--help') for a in sys.argv)


def _validate_exists(path):
    if not _is_help_request() and path and not os.path.exists(path):
        raise click.BadParameter(f"路径不存在|Path does not exist: {path}")
    return path


@click.command(short_help="转录组变异检测(到VCF)|RNA-seq variant calling (to VCF)",
               context_settings={'show_default': True})
@click.option('-g', '--genome', required=True, callback=lambda c, p, v: _validate_exists(v),
              help='参考基因组 FASTA|Reference genome FASTA')
@click.option('--gff3', callback=lambda c, p, v: _validate_exists(v),
              help='参考 GFF3(可选,HISAT2 剪接位点;不提供则 HISAT2 de novo 发现 junction)|Reference GFF3 (optional, HISAT2 splice sites)')
@click.option('-i', '--input', callback=lambda c, p, v: _validate_exists(v),
              help='原始 FASTQ 目录(跑 fastp)|Raw FASTQ dir (runs fastp)')
@click.option('--clean-fastq-dir', callback=lambda c, p, v: _validate_exists(v),
              help='已清洗 FASTQ 目录(跳过 QC)|Clean FASTQ dir (skip QC)')
@click.option('-o', '--output-dir', default='.', help='输出目录|Output dir')
@click.option('-t', '--threads', type=int, default=12, help='线程数|Threads')
@click.option('--min-conf', type=int, default=20, help='HaplotypeCaller 最小置信度|min confidence')
@click.option('-s', '--step', type=click.IntRange(0, 4),
              help='0=仅建索引|index only;省略=全流程|omit for full pipeline')
@click.option('--no-checkpoint', is_flag=True, help='关闭断点续传|Disable checkpoint')
@click.option('-f', '--force', is_flag=True, help='忽略断点重跑|Force rerun')
@click.option('--dry-run', is_flag=True, help='只打印命令|Dry run')
@click.option('--skip-qc', is_flag=True, help='跳过 fastp|Skip QC')
@click.option('--log-file', help='日志文件|Log file')
@click.option('--log-level', default='INFO', help='日志级别|Log level')
def rnaseq2vcf(genome, gff3, input, clean_fastq_dir, output_dir, threads,
               min_conf, step, no_checkpoint, force, dry_run, skip_qc,
               log_file, log_level):
    """转录组变异检测全流程(到 VCF;ANNOVAR 注释请手动运行 biopytools annovar)|
    RNA-seq variant calling pipeline (to VCF; run `biopytools annovar` manually for annotation)

    示例|Examples: biopytools rnaseq2vcf -g genome.fa --gff3 anno.gff3 -i reads/ -o out/
    """
    main = _lazy_import_main()
    args = ['rnaseq2vcf.py', '-g', genome, '--gff3', gff3, '-o', output_dir,
            '-t', str(threads), '--min-conf', str(min_conf), '--log-level', log_level]
    if input:
        args += ['-i', input]
    if clean_fastq_dir:
        args += ['--clean-fastq-dir', clean_fastq_dir]
    if step is not None:
        args += ['-s', str(step)]
    if no_checkpoint:
        args.append('--no-checkpoint')
    if force:
        args.append('-f')
    if dry_run:
        args.append('--dry-run')
    if skip_qc:
        args.append('--skip-qc')
    if log_file:
        args += ['--log-file', log_file]

    original = sys.argv
    sys.argv = args
    try:
        main()
    except SystemExit as e:
        sys.exit(e.code)
    finally:
        sys.argv = original
