"""
Minigraph泛基因组图构建和分析命令|Minigraph Pangenome Graph Construction and Analysis Command
"""

import click
import sys
import os


def _lazy_import_main():
    """延迟加载minigraph主函数|Lazy load minigraph main function"""
    try:
        from ....minigraph.main import main as minigraph_main
        return minigraph_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _validate_file_exists(file_path):
    """验证文件存在(仅在非帮助模式)|Validate file existence (only in non-help mode)"""
    help_flags = {'-h', '--help'}
    is_help = any(arg in help_flags for arg in sys.argv)

    if not is_help and file_path and not os.path.exists(file_path):
        raise click.BadParameter(f"文件不存在|File does not exist: {file_path}")
    return file_path


@click.group(
    short_help='Minigraph泛基因组图工具|Minigraph pangenome graph tool',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
def minigraph():
    """
    Minigraph泛基因组图构建和分析工具|Minigraph Pangenome Graph Construction and Analysis Tool

    构建泛基因组图、序列映射、SV调用和bubble提取|Build pangenome graphs, map sequences, call SVs and extract bubbles
    """
    pass


@minigraph.command(short_help='构建泛基因组图|Build pangenome graph')
@click.option('--ref',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='参考基因组FASTA文件|Reference genome FASTA file')
@click.option('--samples',
              required=True,
              multiple=True,
              callback=lambda ctx, param, value: [_validate_file_exists(f) for f in value] if value else [],
              help='样本基因组FASTA文件|Sample genome FASTA files')
@click.option('-o', '--output-gfa',
              default='./pangenome.gfa',
              show_default=True,
              help='输出GFA文件路径|Output GFA file path')
@click.option('--preset',
              default='ggs',
              type=click.Choice(['g', 'gs', 'ggs']),
              show_default=True,
              help='图构建预设|Graph building preset')
@click.option('--min-identity',
              type=float,
              default=0.9,
              show_default=True,
              help='最小序列相似度|Minimum sequence identity')
@click.option('--min-aln-len',
              type=int,
              default=100000,
              show_default=True,
              help='最小比对长度|Minimum alignment length')
@click.option('--max-gap',
              type=int,
              default=1000000,
              show_default=True,
              help='最大gap大小|Maximum gap size')
@click.option('-t', '--threads',
              type=int,
              default=16,
              show_default=True,
              help='线程数|Number of threads')
@click.option('--batch-size',
              type=int,
              help='批处理大小(MB)|Batch size (MB)')
@click.option('--minigraph-path',
              default='minigraph',
              show_default=True,
              help='minigraph工具路径|minigraph tool path')
@click.option('--gfatools-path',
              default='gfatools',
              show_default=True,
              help='gfatools工具路径|gfatools tool path')
@click.option('--keep-intermediate',
              is_flag=True,
              help='保留中间文件|Keep intermediate files')
@click.option('--append-mode',
              is_flag=True,
              help='追加模式|Append mode')
def build(ref, samples, output_gfa, preset, min_identity, min_aln_len, max_gap,
          threads, batch_size, minigraph_path, gfatools_path, keep_intermediate, append_mode):
    """
    构建泛基因组图|Build pangenome graph from multiple genomes

    示例|Example: biopytools minigraph build --ref ref.fa --samples s1.fa s2.fa -o graph.gfa
    """
    minigraph_main = _lazy_import_main()

    # 构建参数列表|Build argument list
    args = ['minigraph.py', 'build']

    args.extend(['--ref', ref])
    args.extend(['--samples'] + list(samples))
    args.extend(['--output-gfa', output_gfa])
    args.extend(['--preset', preset])

    if min_identity != 0.9:
        args.extend(['--min-identity', str(min_identity)])

    if min_aln_len != 100000:
        args.extend(['--min-aln-len', str(min_aln_len)])

    if max_gap != 1000000:
        args.extend(['--max-gap', str(max_gap)])

    if threads != 16:
        args.extend(['--threads', str(threads)])

    if batch_size:
        args.extend(['--batch-size', str(batch_size)])

    if minigraph_path != 'minigraph':
        args.extend(['--minigraph-path', minigraph_path])

    if gfatools_path != 'gfatools':
        args.extend(['--gfatools-path', gfatools_path])

    if keep_intermediate:
        args.append('--keep-intermediate')

    if append_mode:
        args.append('--append-mode')

    # 执行主程序|Execute main program
    original_argv = sys.argv
    sys.argv = args

    try:
        minigraph_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv


@minigraph.command(short_help='调用结构变异|Call structural variants')
@click.option('--graph-gfa',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='泛基因组图GFA文件|Pangenome graph GFA file')
@click.option('--samples',
              required=True,
              multiple=True,
              callback=lambda ctx, param, value: [_validate_file_exists(f) for f in value] if value else [],
              help='样本基因组FASTA文件|Sample genome FASTA files')
@click.option('-o', '--output-dir',
              default='./minigraph_call',
              show_default=True,
              help='输出目录|Output directory')
@click.option('--preset',
              default='asm',
              type=click.Choice(['asm']),
              show_default=True,
              help='SV调用预设|SV calling preset')
@click.option('-t', '--threads',
              type=int,
              default=16,
              show_default=True,
              help='线程数|Number of threads')
@click.option('--minigraph-path',
              default='minigraph',
              show_default=True,
              help='minigraph工具路径|minigraph tool path')
def call(graph_gfa, samples, output_dir, preset, threads, minigraph_path):
    """
    为样本生成BED文件|Generate BED files for samples

    示例|Example: biopytools minigraph call --graph-gfa graph.gfa --samples s1.fa s2.fa
    """
    minigraph_main = _lazy_import_main()

    # 构建参数列表|Build argument list
    args = ['minigraph.py', 'call']

    args.extend(['--graph-gfa', graph_gfa])
    args.extend(['--samples'] + list(samples))
    args.extend(['--output-dir', output_dir])
    args.extend(['--preset', preset])

    if threads != 16:
        args.extend(['--threads', str(threads)])

    if minigraph_path != 'minigraph':
        args.extend(['--minigraph-path', minigraph_path])

    # 执行主程序|Execute main program
    original_argv = sys.argv
    sys.argv = args

    try:
        minigraph_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv


@minigraph.command(short_help='提取SV bubbles|Extract SV bubbles')
@click.option('--graph-gfa',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='泛基因组图GFA文件|Pangenome graph GFA file')
@click.option('-o', '--output-bed',
              default='./sv_bubbles.bed',
              show_default=True,
              help='输出BED文件路径|Output BED file path')
@click.option('--gfatools-path',
              default='gfatools',
              show_default=True,
              help='gfatools工具路径|gfatools tool path')
def bubble(graph_gfa, output_bed, gfatools_path):
    """
    提取结构变异bubble|Extract structural variant bubbles from graph

    示例|Example: biopytools minigraph bubble --graph-gfa graph.gfa -o bubbles.bed
    """
    minigraph_main = _lazy_import_main()

    # 构建参数列表|Build argument list
    args = ['minigraph.py', 'bubble']

    args.extend(['--graph-gfa', graph_gfa])
    args.extend(['--output-bed', output_bed])

    if gfatools_path != 'gfatools':
        args.extend(['--gfatools-path', gfatools_path])

    # 执行主程序|Execute main program
    original_argv = sys.argv
    sys.argv = args

    try:
        minigraph_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv


@minigraph.command(short_help='序列到图映射|Map sequences to graph')
@click.option('--graph-gfa',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='泛基因组图GFA文件|Pangenome graph GFA file')
@click.option('--queries',
              required=True,
              multiple=True,
              callback=lambda ctx, param, value: [_validate_file_exists(f) for f in value] if value else [],
              help='查询序列FASTA文件|Query sequence FASTA files')
@click.option('-o', '--output-gaf',
              default='./mapping.gaf',
              show_default=True,
              help='输出GAF文件路径|Output GAF file path')
@click.option('--preset',
              default='lr',
              type=click.Choice(['sr', 'lr', 'map-pb', 'map-ont', 'asm']),
              show_default=True,
              help='映射预设|Mapping preset')
@click.option('--max-intron-len',
              type=int,
              help='最大内含子长度|Maximum intron length')
@click.option('-t', '--threads',
              type=int,
              default=16,
              show_default=True,
              help='线程数|Number of threads')
@click.option('--batch-size',
              type=int,
              help='批处理大小(MB)|Batch size (MB)')
@click.option('--minigraph-path',
              default='minigraph',
              show_default=True,
              help='minigraph工具路径|minigraph tool path')
def map(graph_gfa, queries, output_gaf, preset, max_intron_len, threads, batch_size, minigraph_path):
    """
    映射序列到泛基因组图|Map sequences to pangenome graph

    示例|Example: biopytools minigraph map --graph-gfa graph.gfa --queries reads.fa -o mapping.gaf
    """
    minigraph_main = _lazy_import_main()

    # 构建参数列表|Build argument list
    args = ['minigraph.py', 'map']

    args.extend(['--graph-gfa', graph_gfa])
    args.extend(['--queries'] + list(queries))
    args.extend(['--output-gaf', output_gaf])
    args.extend(['--preset', preset])

    if max_intron_len:
        args.extend(['--max-intron-len', str(max_intron_len)])

    if threads != 16:
        args.extend(['--threads', str(threads)])

    if batch_size:
        args.extend(['--batch-size', str(batch_size)])

    if minigraph_path != 'minigraph':
        args.extend(['--minigraph-path', minigraph_path])

    # 执行主程序|Execute main program
    original_argv = sys.argv
    sys.argv = args

    try:
        minigraph_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
