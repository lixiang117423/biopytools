"""
Minibwa短读长比对分析命令|Minibwa Short-read Alignment Command
"""

import click
import sys
import os


def _lazy_import_minibwa_main():
    """延迟加载minibwa主函数|Lazy load minibwa main function"""
    try:
        from ...minibwa.main import main as minibwa_main
        return minibwa_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_path_exists(path):
    """验证路径存在(仅在非帮助模式)|Validate path existence (only in non-help mode)"""
    if not _is_help_request() and not os.path.exists(path):
        raise click.BadParameter(f"路径不存在|Path does not exist: {path}")
    return path


@click.command(
    short_help='Minibwa短读长比对|Minibwa short-read alignment',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('-g', '--genome',
              required=True,
              callback=lambda ctx, param, value: _validate_path_exists(value) if value else None,
              help='参考基因组FASTA|Reference genome FASTA')
@click.option('-i', '--input',
              required=True,
              callback=lambda ctx, param, value: _validate_path_exists(value) if value else None,
              help='FASTQ输入目录|FASTQ input directory')
@click.option('-p', '--pattern',
              default='_1.fq.gz',
              show_default=True,
              help='R1匹配模式|R1 pattern')
@click.option('-o', '--output-dir',
              default='./minibwa_output',
              show_default=True,
              help='输出目录|Output directory')
@click.option('--mode',
              type=click.Choice(['standard', 'hic', 'meth', 'long']),
              default='standard',
              show_default=True,
              help='比对模式|Alignment mode')
@click.option('-t', '--threads',
              type=int, default=12, show_default=True,
              help='线程数|Number of threads')
# Minibwa map 参数|Minibwa map parameters
@click.option('--preset',
              default='adap', show_default=True,
              help='-x 预设|preset (sr/lr/adap)')
@click.option('-k', '--min-seed',
              type=int, default=19, show_default=True,
              help='-k 最小种子长度|min seed length')
@click.option('-c', '--max-occ',
              type=int, default=250, show_default=True,
              help='-c 最大种子出现次数|max seed occurrences')
@click.option('--max-gap',
              type=int, default=100, show_default=True,
              help='minibwa -g 最大gap|max gap size')
@click.option('-w', '--bandwidth',
              type=int, default=100, show_default=True,
              help='-w 带宽|bandwidth')
@click.option('-W', '--bandwidth-long',
              type=int, default=20000, show_default=True,
              help='-W 长读带宽|long bandwidth')
@click.option('-m', '--min-chain-score',
              type=int, default=25, show_default=True,
              help='minibwa -m 最小链接分数|min chaining score')
@click.option('--sec-ratio',
              type=float, default=0.5, show_default=True,
              help='-p 次要/主要得分比|secondary-to-primary ratio')
@click.option('-N', '--max-sec',
              type=int, default=50, show_default=True,
              help='-N 保留次要比对数|retain N secondary alignments')
@click.option('-s', '--min-dp-score',
              type=int, default=30, show_default=True,
              help='-s 最小DP得分|min DP score')
# 打分参数|Scoring
@click.option('-A', '--match-score',
              type=int, default=2, show_default=True,
              help='-A 匹配得分|matching score')
@click.option('-B', '--mismatch-penalty',
              type=int, default=8, show_default=True,
              help='-B 错配罚分|mismatch penalty')
@click.option('-O', '--gap-open',
              default='12,23', show_default=True,
              help='-O gap开放罚分|gap open penalty')
@click.option('-E', '--gap-ext',
              default='2,1', show_default=True,
              help='-E gap延伸罚分|gap extension penalty')
# IO参数|I/O
@click.option('-R', '--read-group',
              default=None,
              help='-R SAM read group|read group line')
@click.option('-u', '--no-unmap',
              is_flag=True,
              help='-u 不输出未比对read|do not output unmapped')
@click.option('--outn',
              type=int, default=0, show_default=True,
              help='--outn 输出次要比对上限|max secondary to output')
@click.option('-y', '--copy-comment',
              is_flag=True,
              help='-y 复制FASTA/Q注释|copy comments')
@click.option('-Y', '--soft-clip-supp',
              is_flag=True,
              help='-Y 软剪切补充比对|soft clip supplementary')
@click.option('-K', '--batch-size',
              default='100m,1g', show_default=True,
              help='-K 批处理大小|batch size')
# 后处理|Post-processing
@click.option('--markdup',
              is_flag=True,
              help='标记重复|Mark duplicates')
@click.option('--remove-dup',
              is_flag=True,
              help='移除重复（需--markdup）|Remove duplicates (requires --markdup)')
# 覆盖度|Coverage
@click.option('--skip-coverage',
              is_flag=True,
              help='跳过覆盖度分析|Skip coverage analysis')
@click.option('--min-base-quality',
              type=int, default=0, show_default=True,
              help='最小碱基质量|Min base quality')
@click.option('--min-mapping-quality',
              type=int, default=0, show_default=True,
              help='最小比对质量|Min mapping quality')
@click.option('--max-depth',
              type=int, default=0, show_default=True,
              help='最大深度限制(0=无限)|Max depth (0=unlimited)')
@click.option('--window-size',
              type=int, default=1000000, show_default=True,
              help='窗口大小|Window size bp')
@click.option('--step-size',
              type=int, default=100000, show_default=True,
              help='步长|Step size bp')
# 工具路径|Tool paths
@click.option('--minibwa-path',
              default='~/software/minibwa/minibwa', show_default=True,
              help='minibwa二进制路径|minibwa binary path')
@click.option('--samtools-path',
              default='~/.local/bin/samtools', show_default=True,
              help='samtools二进制路径|samtools binary path')
# 运行控制|Run control
@click.option('--resume',
              is_flag=True,
              help='断点续传|Resume (skip completed samples)')
def minibwa(genome, input, pattern, output_dir, mode, threads,
            preset, min_seed, max_occ, max_gap, bandwidth, bandwidth_long,
            min_chain_score, sec_ratio, max_sec, min_dp_score,
            match_score, mismatch_penalty, gap_open, gap_ext,
            read_group, no_unmap, outn, copy_comment, soft_clip_supp, batch_size,
            markdup, remove_dup,
            skip_coverage, min_base_quality, min_mapping_quality, max_depth,
            window_size, step_size,
            minibwa_path, samtools_path, resume):
    """
    Minibwa短读长比对分析工具|Minibwa Short-read Alignment Tool

    支持标准短读长、Hi-C、BS-seq、长读四种模式，自动构建索引、批量比对、统计、覆盖度分析
    |Supports standard short-read, Hi-C, BS-seq, and long-read modes;
    auto-builds index, batch-aligns, generates stats and coverage

    示例|Example: biopytools minibwa -g genome.fa -i fastq_dir -o results -t 16
    """

    minibwa_main = _lazy_import_minibwa_main()

    # 构建参数列表|Build argument list
    args = ['minibwa.py']
    args.extend(['-g', genome])
    args.extend(['-i', input])

    # 仅非默认值传递，保持命令精简|Pass non-defaults only
    if pattern != '_1.fq.gz':
        args.extend(['-p', pattern])
    if output_dir != './minibwa_output':
        args.extend(['-o', output_dir])
    if mode != 'standard':
        args.extend(['--mode', mode])
    if threads != 12:
        args.extend(['-t', str(threads)])

    if preset != 'adap':
        args.extend(['--preset', preset])
    if min_seed != 19:
        args.extend(['-k', str(min_seed)])
    if max_occ != 250:
        args.extend(['-c', str(max_occ)])
    if max_gap != 100:
        args.extend(['-g', str(max_gap)])
    if bandwidth != 100:
        args.extend(['-w', str(bandwidth)])
    if bandwidth_long != 20000:
        args.extend(['-W', str(bandwidth_long)])
    if min_chain_score != 25:
        args.extend(['-m', str(min_chain_score)])
    if sec_ratio != 0.5:
        args.extend(['--sec-ratio', str(sec_ratio)])
    if max_sec != 50:
        args.extend(['-N', str(max_sec)])
    if min_dp_score != 30:
        args.extend(['-s', str(min_dp_score)])

    if match_score != 2:
        args.extend(['-A', str(match_score)])
    if mismatch_penalty != 8:
        args.extend(['-B', str(mismatch_penalty)])
    if gap_open != '12,23':
        args.extend(['-O', gap_open])
    if gap_ext != '2,1':
        args.extend(['-E', gap_ext])

    if read_group:
        args.extend(['-R', read_group])
    if no_unmap:
        args.append('-u')
    if outn != 0:
        args.extend(['--outn', str(outn)])
    if copy_comment:
        args.append('-y')
    if soft_clip_supp:
        args.append('-Y')
    if batch_size != '100m,1g':
        args.extend(['-K', batch_size])

    if markdup:
        args.append('--markdup')
    if remove_dup:
        args.append('--remove-dup')

    if skip_coverage:
        args.append('--skip-coverage')
    if min_base_quality != 0:
        args.extend(['--min-base-quality', str(min_base_quality)])
    if min_mapping_quality != 0:
        args.extend(['--min-mapping-quality', str(min_mapping_quality)])
    if max_depth != 0:
        args.extend(['--max-depth', str(max_depth)])
    if window_size != 1000000:
        args.extend(['--window-size', str(window_size)])
    if step_size != 100000:
        args.extend(['--step-size', str(step_size)])

    if minibwa_path != '~/software/minibwa/minibwa':
        args.extend(['--minibwa-path', minibwa_path])
    if samtools_path != '~/.local/bin/samtools':
        args.extend(['--samtools-path', samtools_path])

    if resume:
        args.append('--resume')

    original_argv = sys.argv
    sys.argv = args

    try:
        minibwa_main()
    except SystemExit as e:
        sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n用户中断|User interrupted", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
