"""
BWA比对分析命令|BWA Alignment Analysis Command
"""

import click
import sys
import os


def _lazy_import_bwa_main():
    """延迟加载BWA主函数|Lazy load BWA main function"""
    try:
        from ...bwa.main import main as bwa_main
        return bwa_main
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
    short_help='BWA比对分析|BWA alignment analysis',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--genome', '-g',
              required=True,
              callback=lambda ctx, param, value: _validate_path_exists(value) if value else None,
              help='参考基因组文件|Reference genome file')
@click.option('--input', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_path_exists(value) if value else None,
              help='输入FASTQ目录|Input FASTQ directory')
@click.option('--pattern', '-p',
              required=True,
              type=str,
              help='FASTQ文件匹配模式|FASTQ file pattern')
@click.option('--output-dir', '-o',
              default='./bwa_output',
              show_default=True,
              type=click.Path(),
              help='输出目录|Output directory')
@click.option('--threads', '-t',
              type=int,
              default=12,
              show_default=True,
              help='线程数|Number of threads')
@click.option('--bwa-k',
              type=int,
              default=19,
              show_default=True,
              help='最小种子长度|Minimum seed length')
@click.option('--bwa-w',
              type=int,
              default=100,
              show_default=True,
              help='带宽|Band width')
@click.option('--bwa-d',
              type=int,
              default=100,
              show_default=True,
              help='X-dropoff|Off-diagonal X-dropoff')
@click.option('--bwa-r',
              type=float,
              default=1.5,
              show_default=True,
              help='内部种子因子|Internal seed factor')
@click.option('--bwa-c',
              type=int,
              default=500,
              show_default=True,
              help='种子出现次数阈值|Seed occurrence threshold')
@click.option('--bwa-D', 'bwa_drop_ratio',
              type=float,
              default=0.50,
              show_default=True,
              help='短链丢弃比例|Short chain drop fraction')
@click.option('--bwa-W', 'bwa_min_chain',
              type=int,
              default=0,
              show_default=True,
              help='最小链长|Minimum chain length')
@click.option('--bwa-m',
              type=int,
              default=50,
              show_default=True,
              help='配对拯救轮数|Mate rescue rounds')
@click.option('--bwa-S', 'bwa_s',
              is_flag=True,
              help='跳过配对拯救|Skip mate rescue')
@click.option('--bwa-P', 'bwa_p',
              is_flag=True,
              help='跳过配对|Skip pairing')
@click.option('--bwa-A', 'bwa_score_match',
              type=int,
              default=1,
              show_default=True,
              help='匹配得分|Match score')
@click.option('--bwa-B', 'bwa_score_mismatch',
              type=int,
              default=4,
              show_default=True,
              help='错配罚分|Mismatch penalty')
@click.option('--bwa-O', 'bwa_gap_open',
              default="6,6",
              show_default=True,
              help='Gap开放罚分|Gap open penalty')
@click.option('--bwa-E', 'bwa_gap_ext',
              default="1,1",
              show_default=True,
              help='Gap延伸罚分|Gap extension penalty')
@click.option('--bwa-L', 'bwa_clip',
              default="5,5",
              show_default=True,
              help='末端剪切罚分|Clipping penalty')
@click.option('--bwa-U', 'bwa_unpaired',
              type=int,
              default=17,
              show_default=True,
              help='未配对罚分|Unpaired penalty')
@click.option('--bwa-M', 'bwa_mark_secondary',
              is_flag=True,
              help='标记次要比对|Mark shorter split hits as secondary')
@click.option('--bwa-T', 'bwa_min_score',
              type=int,
              default=30,
              show_default=True,
              help='最小输出得分|Minimum score to output')
@click.option('--bwa-a', 'bwa_all_align',
              is_flag=True,
              help='输出所有比对|Output all alignments')
@click.option('--bwa-C', 'bwa_append_comment',
              is_flag=True,
              help='附加FASTQ注释|Append FASTA/FASTQ comment')
@click.option('--bwa-V', 'bwa_ref_header',
              is_flag=True,
              help='输出参考序列头|Output reference FASTA header')
@click.option('--bwa-Y', 'bwa_soft_clip',
              is_flag=True,
              help='软剪切补充比对|Soft clipping for supplementary alignments')
@click.option('--markdup',
              is_flag=True,
              help='标记重复序列|Mark duplicate reads')
@click.option('--remove-dup',
              is_flag=True,
              help='移除重复序列|Remove duplicate reads')
@click.option('--min-base-quality',
              type=int,
              default=0,
              show_default=True,
              help='最小碱基质量|Minimum base quality')
@click.option('--min-mapping-quality',
              type=int,
              default=0,
              show_default=True,
              help='最小比对质量|Minimum mapping quality')
@click.option('--max-depth',
              type=int,
              default=0,
              show_default=True,
              help='最大深度限制|Max depth limit')
@click.option('--window-size',
              type=int,
              default=1000000,
              show_default=True,
              help='窗口大小|Window size in bp')
@click.option('--step-size',
              type=int,
              default=100000,
              show_default=True,
              help='步长|Step size in bp')
@click.option('--resume',
              is_flag=True,
              help='断点续传|Resume skip completed samples')
@click.option('--keep-sam',
              is_flag=True,
              help='保留SAM文件|Keep SAM files')
def bwa(genome, input, pattern, output_dir, threads,
        bwa_k, bwa_w, bwa_d, bwa_r, bwa_c,
        bwa_drop_ratio, bwa_min_chain, bwa_m, bwa_s, bwa_p,
        bwa_score_match, bwa_score_mismatch, bwa_gap_open, bwa_gap_ext, bwa_clip, bwa_unpaired,
        bwa_mark_secondary, bwa_min_score, bwa_all_align, bwa_append_comment, bwa_ref_header, bwa_soft_clip,
        markdup, remove_dup, min_base_quality, min_mapping_quality, max_depth,
        window_size, step_size, resume, keep_sam):
    """
    BWA比对分析工具|BWA Alignment Analysis Tool

    BWA-MEM全基因组比对分析|BWA-MEM whole genome alignment analysis

    示例|Examples: biopytools bwa -g genome.fa -i fastq_dir -p _1.clean.fq.gz
    """

    # 延迟加载|Lazy loading
    bwa_main = _lazy_import_bwa_main()

    # 构建参数列表|Build argument list
    args = ['bwa_align.py']
    args.extend(['-g', genome])
    args.extend(['-i', input])
    args.extend(['-p', pattern])

    if output_dir != './bwa_output':
        args.extend(['-o', output_dir])
    if threads != 88:
        args.extend(['-t', str(threads)])

    # BWA算法参数|BWA algorithm parameters
    if bwa_k != 19:
        args.extend(['--bwa-k', str(bwa_k)])
    if bwa_w != 100:
        args.extend(['--bwa-w', str(bwa_w)])
    if bwa_d != 100:
        args.extend(['--bwa-d', str(bwa_d)])
    if bwa_r != 1.5:
        args.extend(['--bwa-r', str(bwa_r)])
    if bwa_c != 500:
        args.extend(['--bwa-c', str(bwa_c)])
    if bwa_drop_ratio != 0.50:
        args.extend(['--bwa-D', str(bwa_drop_ratio)])
    if bwa_min_chain != 0:
        args.extend(['--bwa-W', str(bwa_min_chain)])
    if bwa_m != 50:
        args.extend(['--bwa-m', str(bwa_m)])
    if bwa_s:
        args.append('--bwa-S')
    if bwa_p:
        args.append('--bwa-P')

    # BWA打分参数|BWA scoring parameters
    if bwa_score_match != 1:
        args.extend(['--bwa-A', str(bwa_score_match)])
    if bwa_score_mismatch != 4:
        args.extend(['--bwa-B', str(bwa_score_mismatch)])
    if bwa_gap_open != "6,6":
        args.extend(['--bwa-O', bwa_gap_open])
    if bwa_gap_ext != "1,1":
        args.extend(['--bwa-E', bwa_gap_ext])
    if bwa_clip != "5,5":
        args.extend(['--bwa-L', bwa_clip])
    if bwa_unpaired != 17:
        args.extend(['--bwa-U', str(bwa_unpaired)])

    # BWA输出参数|BWA output parameters
    if bwa_mark_secondary:
        args.append('--bwa-M')
    if bwa_min_score != 30:
        args.extend(['--bwa-T', str(bwa_min_score)])
    if bwa_all_align:
        args.append('--bwa-a')
    if bwa_append_comment:
        args.append('--bwa-C')
    if bwa_ref_header:
        args.append('--bwa-V')
    if bwa_soft_clip:
        args.append('--bwa-Y')

    # 后处理参数|Post-processing parameters
    if markdup:
        args.append('--markdup')
    if remove_dup:
        args.append('--remove-dup')

    # 覆盖度参数|Coverage parameters
    if min_base_quality != 0:
        args.extend(['--min-base-quality', str(min_base_quality)])
    if min_mapping_quality != 0:
        args.extend(['--min-mapping-quality', str(min_mapping_quality)])
    if max_depth != 0:
        args.extend(['--max-depth', str(max_depth)])

    # 滑窗参数|Window parameters
    if window_size != 1000000:
        args.extend(['--window-size', str(window_size)])
    if step_size != 100000:
        args.extend(['--step-size', str(step_size)])

    # 其他参数|Other parameters
    if resume:
        args.append('--resume')
    if keep_sam:
        args.append('--keep-sam')

    # 保存和恢复sys.argv|Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        # 调用原始main函数|Call original main function
        bwa_main()
    except SystemExit as e:
        if e.code != 0:
            sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n用户中断|User interrupted", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
