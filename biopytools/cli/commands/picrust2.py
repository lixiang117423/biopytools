"""
PICRUSt2微生物群落功能丰度预测|PICRUSt2 Microbial Community Functional Abundance Prediction Command
"""

import click
import sys
import os


def _lazy_import_picrust2_main():
    """延迟加载picrust2主函数|Lazy load picrust2 main function"""
    try:
        from ...picrust2.main import main as picrust2_main
        return picrust2_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_file_exists(file_path):
    """验证文件存在|Validate file existence"""
    if not _is_help_request() and file_path and not os.path.exists(file_path):
        raise click.BadParameter(f"文件不存在|File does not exist: {file_path}")
    return file_path


@click.command(
    short_help='PICRUSt2微生物群落功能丰度预测|PICRUSt2 functional abundance prediction',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('-s', '--study-fasta',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='代表序列FASTA文件|FASTA of unaligned study sequences')
@click.option('-i', '--input',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='特征表文件(自动识别BIOM/TSV/Excel/Mothur)|Input table of sequence abundances')
@click.option('-o', '--output-dir',
              default='./picrust2_output',
              show_default=True,
              help='输出目录|Output directory')
@click.option('-t', '--threads',
              type=int,
              default=12,
              show_default=True,
              help='线程数|Number of threads')
@click.option('--max-nsti',
              type=float,
              default=2.0,
              show_default=True,
              help='最大NSTI阈值|Maximum NSTI value')
@click.option('--stratified',
              is_flag=True,
              help='生成分层输出表|Generate stratified output tables')
@click.option('--in-traits',
              default='EC,KO',
              show_default=True,
              help='功能数据库(EC,KO,GO,PFAM,BIGG,CAZY)|Gene families to predict')
@click.option('--placement-tool',
              type=click.Choice(['epa-ng', 'sepp']),
              default='epa-ng',
              show_default=True,
              help='序列放置工具|Placement tool')
@click.option('--hsp-method',
              type=click.Choice(['mp', 'emp_prob', 'pic', 'scp', 'subtree_average']),
              default='mp',
              show_default=True,
              help='隐状态预测方法|HSP method')
@click.option('--edge-exponent',
              type=float,
              default=0.5,
              show_default=True,
              help='HSP edge exponent')
@click.option('--pipeline',
              type=click.Choice(['auto', 'split', 'single']),
              default='auto',
              show_default=True,
              help='流程类型: auto/split/single|Pipeline type')
@click.option('--min-align',
              type=float,
              default=0.8,
              show_default=True,
              help='最小比对比例|Minimum alignment proportion')
@click.option('--min-reads',
              type=int,
              default=1,
              show_default=True,
              help='每ASV最小reads数|Minimum reads per ASV')
@click.option('--min-samples',
              type=int,
              default=1,
              show_default=True,
              help='每ASV最小样本数|Minimum samples per ASV')
@click.option('--no-pathways',
              is_flag=True,
              help='跳过通路推断|Skip pathway inference')
@click.option('--coverage',
              is_flag=True,
              help='计算通路覆盖度|Calculate pathway coverages')
@click.option('--skip-minpath',
              is_flag=True,
              help='跳过MinPath|Skip MinPath')
@click.option('--no-gap-fill',
              is_flag=True,
              help='跳过gap filling|Skip gap filling')
@click.option('--per-sequence-contrib',
              is_flag=True,
              help='逐序列运行MinPath|Run MinPath per sequence')
@click.option('--skip-norm',
              is_flag=True,
              help='跳过归一化|Skip normalization')
@click.option('--remove-intermediate',
              is_flag=True,
              help='移除中间文件|Remove intermediate files')
@click.option('--verbose',
              is_flag=True,
              help='详细输出|Verbose output')
def picrust2(study_fasta, input, output_dir, threads, max_nsti, stratified,
             in_traits, placement_tool, hsp_method, edge_exponent, pipeline,
             min_align, min_reads, min_samples, no_pathways, coverage,
             skip_minpath, no_gap_fill, per_sequence_contrib, skip_norm,
             remove_intermediate, verbose):
    """
    PICRUSt2微生物群落功能丰度预测|PICRUSt2 Microbial Community Functional Abundance Prediction

    基于16S rRNA标记基因序列预测微生物群落功能丰度和代谢通路|Predict functional abundances and metabolic pathways from marker gene sequences

    示例|Example: biopytools picrust2 -s study_seqs.fna -i seqabun.biom -o picrust2_out
    """

    picrust2_main = _lazy_import_picrust2_main()

    args = ['picrust2.py']
    args.extend(['-s', study_fasta])
    args.extend(['-i', input])

    if output_dir != './picrust2_output':
        args.extend(['-o', output_dir])

    if threads != 12:
        args.extend(['-t', str(threads)])

    if max_nsti != 2.0:
        args.extend(['--max-nsti', str(max_nsti)])

    if stratified:
        args.append('--stratified')

    if in_traits != 'EC,KO':
        args.extend(['--in-traits', in_traits])

    if placement_tool != 'epa-ng':
        args.extend(['--placement-tool', placement_tool])

    if hsp_method != 'mp':
        args.extend(['--hsp-method', hsp_method])

    if edge_exponent != 0.5:
        args.extend(['--edge-exponent', str(edge_exponent)])

    if pipeline != 'auto':
        args.extend(['--pipeline', pipeline])

    if min_align != 0.8:
        args.extend(['--min-align', str(min_align)])

    if min_reads != 1:
        args.extend(['--min-reads', str(min_reads)])

    if min_samples != 1:
        args.extend(['--min-samples', str(min_samples)])

    if no_pathways:
        args.append('--no-pathways')

    if coverage:
        args.append('--coverage')

    if skip_minpath:
        args.append('--skip-minpath')

    if no_gap_fill:
        args.append('--no-gap-fill')

    if per_sequence_contrib:
        args.append('--per-sequence-contrib')

    if skip_norm:
        args.append('--skip-norm')

    if remove_intermediate:
        args.append('--remove-intermediate')

    if verbose:
        args.append('--verbose')

    original_argv = sys.argv
    sys.argv = args

    try:
        picrust2_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
