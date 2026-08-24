"""
NLR-Annotator分析命令|NLR-Annotator Analysis Command
"""

import click
import sys
import os


def _lazy_import_nlr_annotator_main():
    """延迟加载nlr_annotator主函数|Lazy load nlr_annotator main function"""
    try:
        from ...nlr_annotator.main import main as nlr_annotator_main
        return nlr_annotator_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_path_exists(path):
    """验证路径存在(仅在非帮助模式)|Validate path exists (only in non-help mode)"""
    if not _is_help_request() and path and not os.path.exists(os.path.expanduser(path)):
        raise click.BadParameter(f"路径不存在|Path does not exist: {path}")
    return path


@click.command(
    short_help='NLR基因预测工具(支持DNA/CDS)|NLR gene prediction tool (DNA/CDS)',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('-i', '--input',
              required=True,
              callback=lambda ctx, param, value: _validate_path_exists(value) if value else None,
              help='输入DNA/CDS FASTA文件或目录|Input DNA/CDS FASTA file or directory')
@click.option('-o', '--output-dir',
              default='./output',
              help='输出目录|Output directory')
@click.option('-t', '--threads',
              type=int, default=12, show_default=True,
              help='线程数|Number of threads')
@click.option('--sample-suffix',
              default='*.fa', show_default=True,
              help='目录模式下文件匹配后缀|File match suffix for directory mode')
@click.option('--merge-only',
              is_flag=True,
              help='只合并已有结果TSV(*.nlr_annotator.tsv),不运行NLR-Annotator'
                   '|Merge existing result TSVs only, skip NLR-Annotator')
@click.option('--no-filter-contained',
              is_flag=True,
              help='关闭被包含冗余调用过滤(默认开启:剔除被完整基因完全包含的短片段调用,留档*.removed.tsv)'
                   '|Disable contained-call filtering (default ON)')
@click.option('--output-gff',
              is_flag=True,
              help='输出GFF文件|Output GFF file')
@click.option('--output-bed',
              is_flag=True,
              help='输出BED文件|Output BED file')
@click.option('--output-motifs',
              is_flag=True,
              help='输出motifs BED文件|Output motifs BED file')
@click.option('--output-alignment',
              is_flag=True,
              help='输出motif比对FASTA|Output motif alignment FASTA')
@click.option('--jar-path', default='',
              help='NLR-Annotator JAR文件路径|NLR-Annotator JAR file path')
@click.option('--mot-file', default='',
              help='mot.txt配置文件路径|mot.txt config file path')
@click.option('--store-file', default='',
              help='store.txt配置文件路径|store.txt config file path')
@click.option('--java-path', default='java', show_default=True,
              help='Java解释器路径(conda env用~/miniforge3/envs/xxx/bin/java)|Java interpreter path')
@click.option('--num-seqs-per-thread', type=int, default=1000, show_default=True,
              help='每线程处理序列数|Sequences per thread')
@click.option('--distance-within-motif-combination', type=int, default=500, show_default=True,
              help='motif组合内距离|Distance within motif combination')
@click.option('--distance-for-elongating', type=int, default=2500, show_default=True,
              help='延伸距离|Distance for elongating')
@click.option('--distance-between-motif-combinations', type=int, default=50000, show_default=True,
              help='motif组合间距离|Distance between motif combinations')
def nlr_annotator(input, output_dir, threads, sample_suffix, merge_only, no_filter_contained,
                   output_gff, output_bed, output_motifs, output_alignment, jar_path, mot_file,
                   store_file, java_path, num_seqs_per_thread,
                   distance_within_motif_combination, distance_for_elongating,
                   distance_between_motif_combinations):
    """从DNA/CDS序列预测NLR基因|Predict NLR genes from DNA/CDS sequences

    示例|Example: biopytools nlr-annotator -i genome.cds.fa -o output_dir/
    """
    nlr_main = _lazy_import_nlr_annotator_main()

    argv = ['nlr_annotator.py']
    argv.extend(['-i', input])
    argv.extend(['-o', output_dir])
    argv.extend(['-t', str(threads)])
    argv.extend(['--sample-suffix', sample_suffix])
    if merge_only:
        argv.append('--merge-only')
    if no_filter_contained:
        argv.append('--no-filter-contained')
    if output_gff:
        argv.append('--output-gff')
    if output_bed:
        argv.append('--output-bed')
    if output_motifs:
        argv.append('--output-motifs')
    if output_alignment:
        argv.append('--output-alignment')
    # 工具路径与高级参数透传(默认值不透传,main 用自身 default)|Forward tool paths & advanced params
    if jar_path:
        argv.extend(['--jar-path', jar_path])
    if mot_file:
        argv.extend(['--mot-file', mot_file])
    if store_file:
        argv.extend(['--store-file', store_file])
    if java_path != 'java':
        argv.extend(['--java-path', java_path])
    if num_seqs_per_thread != 1000:
        argv.extend(['--num-seqs-per-thread', str(num_seqs_per_thread)])
    if distance_within_motif_combination != 500:
        argv.extend(['--distance-within-motif-combination', str(distance_within_motif_combination)])
    if distance_for_elongating != 2500:
        argv.extend(['--distance-for-elongating', str(distance_for_elongating)])
    if distance_between_motif_combinations != 50000:
        argv.extend(['--distance-between-motif-combinations', str(distance_between_motif_combinations)])

    original_argv = sys.argv
    sys.argv = argv

    try:
        nlr_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
