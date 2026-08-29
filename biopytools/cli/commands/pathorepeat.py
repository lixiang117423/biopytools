"""pathorepeat 病原菌重复序列注释命令|pathorepeat pathogen repeat command"""

import os
import sys

import click


def _lazy_import_main():
    """延迟加载主函数|Lazy load main function"""
    try:
        from ...pathorepeat.main import main as pathorepeat_main
        return pathorepeat_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    return any(arg in {'-h', '--help'} for arg in sys.argv)


def _validate_path_exists(path):
    """验证路径存在(非帮助模式)|Validate path exists (non-help mode)"""
    if not _is_help_request() and path and not os.path.exists(path):
        raise click.BadParameter(f"路径不存在|Path does not exist: {path}")
    return path


@click.command(
    short_help='病原菌重复序列注释(RepeatModeler2+RepeatMasker+TEsorter)|'
               'Pathogen repeat annotation (RepeatModeler2+RepeatMasker+TEsorter)',
    context_settings=dict(help_option_names=['-h', '--help'],
                          max_content_width=120))
@click.option('--input', '-i', default=None, required=True,
              callback=lambda ctx, param, value: _validate_path_exists(value),
              type=click.Path(),
              help='基因组FASTA或文件夹(批量)|Genome FASTA or directory (batch)')
@click.option('--output-dir', '-o', default='./pathorepeat_output',
              show_default=True, type=click.Path(),
              help='输出目录|Output directory')
@click.option('--threads', '-t', default=12, show_default=True, type=int,
              help='线程数|Thread count')
@click.option('--masking-mode', default='xsmall', show_default=True,
              type=click.Choice(['xsmall', 'soft', 'hard', 'x']),
              help='屏蔽模式(xsmall=小写软屏蔽,病原菌默认)|Masking mode '
                   '(xsmall=lowercase soft mask, pathogen default)')
@click.option('--ltr-struct/--no-ltr-struct', default=True, show_default=True,
              help='RepeatModeler -LTRStruct(默认开)|-LTRStruct (default on)')
@click.option('--tesorter-db', default='rexdb', show_default=True,
              type=click.Choice(['gydb', 'rexdb', 'rexdb-plant', 'rexdb-metazoa',
                                 'rexdb-v3', 'rexdb-plantv3', 'rexdb-metazoav3',
                                 'rexdb-pnas', 'rexdb-line', 'sine']),
              help='TEsorter数据库(REXdb植物/动物为主,卵菌可试gydb)|TEsorter db')
@click.option('--db-hmm', default=None,
              help='自定义TEsorter HMM文件(优先于--tesorter-db)|Custom HMM file')
@click.option('--famdb-dir', default=None,
              callback=lambda ctx, param, value: _validate_path_exists(value),
              type=click.Path(),
              help='Dfam famdb数据目录(注入FAMDB_DIR启用RM2自带分类;不设则分类失败'
                   '自动降级)|Dfam famdb dir (injected as FAMDB_DIR; auto-degrades '
                   'if unset)')
@click.option('--effector-bed', default=None,
              callback=lambda ctx, param, value: _validate_path_exists(value),
              help='effector候选区BED(仅单文件模式)|Effector BED (single-sample)')
@click.option('--effector-gff', default=None,
              callback=lambda ctx, param, value: _validate_path_exists(value),
              help='effector候选区GFF3(仅单文件模式)|Effector GFF3 (single-sample)')
@click.option('--genome-name', default=None,
              help='输出前缀(仅单文件模式)|Output prefix (single-sample only)')
@click.option('--skip-completed/--no-skip-completed', default=True,
              show_default=True,
              help='断点续传(跳过已完成步骤)|Resume (skip completed steps)')
@click.option('--log-level', default='INFO', show_default=True,
              type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'ERROR']),
              help='日志级别|Log level')
def pathorepeat(input, output_dir, threads, masking_mode, ltr_struct,
                tesorter_db, db_hmm, famdb_dir, effector_bed, effector_gff,
                genome_name, skip_completed, log_level):
    """病原菌重复序列注释:RepeatModeler2+RepeatMasker(-xsmall)+TEsorter
    |Pathogen repeat annotation: RepeatModeler2 + RepeatMasker (-xsmall) + TEsorter

    示例|Examples: biopytools pathorepeat -i genome.fa -o out_dir/
    """
    pathorepeat_main = _lazy_import_main()

    # 构造参数列表(默认值显式透传)|Build argv (defaults always forwarded)
    args = ['pathorepeat', '-i', input, '-o', output_dir,
            '-t', str(threads), '--masking-mode', masking_mode,
            '--tesorter-db', tesorter_db]
    if not ltr_struct:
        args.append('--no-ltr-struct')
    if db_hmm:
        args.extend(['--db-hmm', db_hmm])
    if famdb_dir:
        args.extend(['--famdb-dir', famdb_dir])
    if effector_bed:
        args.extend(['--effector-bed', effector_bed])
    if effector_gff:
        args.extend(['--effector-gff', effector_gff])
    if genome_name:
        args.extend(['--genome-name', genome_name])
    if not skip_completed:
        args.append('--no-skip-completed')
    args.extend(['--log-level', log_level])

    original_argv = sys.argv
    sys.argv = args
    try:
        pathorepeat_main()
    except SystemExit as e:
        sys.exit(e.code)
    finally:
        sys.argv = original_argv
