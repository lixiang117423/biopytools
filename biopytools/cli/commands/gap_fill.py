"""
TGS-GapCloser Gap填充命令|TGS-GapCloser Gap Filling Command
"""

import click
import sys
import os
from pathlib import Path


# 从config导入默认值|Import default values from config
def _get_default_tgsgapcloser_path():
    """获取默认TGS-GapCloser路径|Get default TGS-GapCloser path"""
    # 优先使用环境变量|Prefer environment variable
    env_path = os.environ.get('TGSGAPCLOSER_PATH')
    if env_path and os.path.exists(env_path):
        return env_path

    # 使用相对路径|Use relative path
    default_path = Path.home() / 'software/TGS-GapCloser2/TGS-GapCloser2-master/tgsgapcloser2'
    if default_path.exists():
        return str(default_path)

    # Fallback to old path
    old_path = Path.home() / 'software/TGS-GapCloser/tgsgapcloser'
    if old_path.exists():
        return str(old_path)

    # Return default anyway
    return str(default_path)


DEFAULT_TGSGAPCLOSER_PATH = _get_default_tgsgapcloser_path()


def _lazy_import_tgsgapcloser_main():
    """延迟加载tgsgapcloser主函数|Lazy load tgsgapcloser main function"""
    try:
        from ...tgsgapcloser.main import main as tgsgapcloser_main
        return tgsgapcloser_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_file_exists(file_path):
    """验证文件存在(仅在非帮助模式)|Validate file existence (only in non-help mode)"""
    if not _is_help_request() and file_path and not os.path.exists(file_path):
        raise click.BadParameter(f"文件不存在|File does not exist: {file_path}")
    return file_path


@click.command(
    short_help='TGS-GapCloser Gap填充工具|TGS-GapCloser gap filling tool',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('-s', '--scaff-file',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='输入scaffold文件|Input scaffold file')
@click.option('-t', '--tgstype',
              required=True,
              type=click.Choice(['ont', 'pb', 'hifi']),
              help='TGS类型|TGS type (ont/pb/hifi)')
@click.option('-ir', '--reads-file',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='输入TGS reads文件|Input TGS reads file')
@click.option('-o', '--output-prefix',
              required=True,
              help='输出前缀|Output prefix')
@click.option('-m', '--mode',
              type=click.Choice(['none', 'racon', 'pilon']),
              default='none',
              help='纠错模式|Error correction mode (default: none)')
@click.option('--tgsgapcloser-path',
              default=DEFAULT_TGSGAPCLOSER_PATH,
              show_default=False,
              help='TGS-GapCloser路径|TGS-GapCloser path (default: auto-detect)')
@click.option('-idy', '--min-idy',
              type=float,
              help='最小同一性|Min identity (auto-set)')
@click.option('-l', '--min-match',
              type=int,
              help='最小匹配长度|Min match length (auto-set)')
@click.option('-threads', '--threads',
              type=int,
              default=12,
              help='线程数|Number of threads (default: 12)')
@click.option('-chunk',
              type=int,
              default=3,
              help='分块数量|Chunk count (default: 3)')
@click.option('-g-check',
              is_flag=True,
              help='启用Gap大小差异检查|Enable gap size difference check')
@click.option('-min-nread',
              type=int,
              default=1,
              help='最小reads数量|Min read count (default: 1)')
@click.option('-max-nread',
              type=int,
              default=-1,
              help='最大reads数量|Max read count (default: -1)')
@click.option('-max-candidate',
              type=int,
              default=200,
              help='最大候选数|Max candidates (default: 200)')
@click.option('-racon', '--racon-path',
              help='Racon路径|Racon path')
@click.option('-racon-round',
              type=int,
              default=3,
              help='Racon轮数|Racon rounds (default: 3)')
@click.option('-pilon', '--pilon-path',
              help='Pilon路径|Pilon path')
@click.option('-ngs', '--ngs-file',
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='NGS reads文件|NGS reads file')
@click.option('-java', '--java-path',
              help='Java路径|Java path')
@click.option('-samtools', '--samtools-path',
              help='Samtools路径|Samtools path')
@click.option('-pilon-mem',
              default='300G',
              help='Pilon内存|Pilon memory (default: 300G)')
@click.option('-pilon-round',
              type=int,
              default=3,
              help='Pilon轮数|Pilon rounds (default: 3)')
@click.option('-minmap-arg',
              help='自定义minimap2参数|Custom minimap2 arguments')
@click.option('-ug', '--unitig-file',
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='hifiasm unitig/contig文件（第2轮填充）|hifiasm unitig/contig file (round 2)')
@click.option('-fl', '--flanking-len',
              type=int,
              default=5000,
              help='Flanking序列长度（bp）|Flanking sequence length (bp) (default: 5000)')
@click.option('-al', '--min-align-len',
              type=int,
              default=1000,
              help='最小比对长度（bp）|Min alignment length (bp) (default: 1000)')
@click.option('-ai', '--min-identity',
              type=int,
              default=40,
              help='最小比对同一性（%）|Min alignment identity (%) (default: 40)')
@click.option('-mfl', '--max-filling-len',
              type=int,
              default=1000000,
              help='最大填充长度（bp）|Max filling length (bp) (default: 1000000)')
def gap_fill(scaff_file, tgstype, reads_file, output_prefix, mode,
                 tgsgapcloser_path, min_idy, min_match, threads, chunk,
                 g_check, min_nread, max_nread, max_candidate,
                 racon_path, racon_round, pilon_path, ngs_file,
                 java_path, samtools_path, pilon_mem, pilon_round, minmap_arg,
                 unitig_file, flanking_len, min_align_len, min_identity, max_filling_len):
    """
    TGS-GapCloser Gap填充工具|TGS-GapCloser Gap Filling Tool

    使用三代测序数据填充基因组组装中的Gap|Fill gaps in genome assembly using TGS long reads

    示例|Examples: biopytools gap-fill -s s.fa -t ont -ir r.fa -o out
    """

    # 延迟加载|Lazy loading
    tgsgapcloser_main = _lazy_import_tgsgapcloser_main()

    # 构建参数列表|Build argument list
    args = ['tgsgapcloser.py']

    # 必需参数|Required parameters
    args.extend(['-s', scaff_file])
    args.extend(['-t', tgstype])
    args.extend(['-ir', reads_file])
    args.extend(['-o', output_prefix])

    # 可选参数|Optional parameters
    if mode != 'none':
        args.extend(['-m', mode])

    if tgsgapcloser_path != DEFAULT_TGSGAPCLOSER_PATH:
        args.extend(['--tgsgapcloser_path', tgsgapcloser_path])

    if min_idy is not None:
        args.extend(['-idy', str(min_idy)])

    if min_match is not None:
        args.extend(['-l', str(min_match)])

    if threads != 12:
        args.extend(['-threads', str(threads)])

    if chunk != 3:
        args.extend(['-chunk', str(chunk)])

    if g_check:
        args.append('-g_check')

    if min_nread != 1:
        args.extend(['-min_nread', str(min_nread)])

    if max_nread != -1:
        args.extend(['-max_nread', str(max_nread)])

    if max_candidate != 200:
        args.extend(['-max_candidate', str(max_candidate)])

    if racon_path is not None:
        args.extend(['-racon', racon_path])

    if racon_round != 3:
        args.extend(['-racon_round', str(racon_round)])

    if pilon_path is not None:
        args.extend(['-pilon', pilon_path])

    if ngs_file is not None:
        args.extend(['-ngs', ngs_file])

    if java_path is not None:
        args.extend(['-java', java_path])

    if samtools_path is not None:
        args.extend(['-samtools', samtools_path])

    if pilon_mem != '300G':
        args.extend(['-pilon_mem', pilon_mem])

    if pilon_round != 3:
        args.extend(['-pilon_round', str(pilon_round)])

    if minmap_arg is not None:
        args.extend(['-minmap_arg', minmap_arg])

    if unitig_file is not None:
        args.extend(['-unitig_file', unitig_file])

    if flanking_len != 5000:
        args.extend(['-flanking_len', str(flanking_len)])

    if min_align_len != 1000:
        args.extend(['-min_align_len', str(min_align_len)])

    if min_identity != 40:
        args.extend(['-min_identity', str(min_identity)])

    if max_filling_len != 1000000:
        args.extend(['-max_filling_len', str(max_filling_len)])

    # 执行主程序|Execute main program
    original_argv = sys.argv
    sys.argv = args

    try:
        tgsgapcloser_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
