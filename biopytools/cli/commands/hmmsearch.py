"""
HMMsearch分析命令|HMMsearch Analysis Command
"""

import click
import sys
import os


def _lazy_import_hmmsearch_main():
    """延迟加载hmmsearch主函数|Lazy load hmmsearch main function"""
    try:
        from ...hmmsearch.main import main as hmmsearch_main
        return hmmsearch_main
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
    short_help='HMMsearch分析工具|HMMsearch analysis tool',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
# 主输入参数（自动判断模式）|Main input parameter (auto-detect mode)
@click.option('-i', '--input',
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='输入文件：domtblout文件（模式1）或HMM文件（模式2，需同时指定-p）|Input file: domtblout (mode 1) or HMM file (mode 2, requires -p)')
@click.option('-p', '--protein-fasta',
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='蛋白序列FASTA文件或目录（模式2必需，模式1提取序列时需要）|Protein FASTA file or directory (required for mode 2, needed for mode 1 if extracting sequences)')
# 输出配置|Output configuration
@click.option('-o', '--output-dir',
              default='./hmmsearch_output',
              show_default=True,
              help='输出目录|Output directory')
@click.option('--output-prefix',
              default='hmmsearch_results',
              show_default=True,
              help='输出文件前缀|Output file prefix')
# HMMsearch软件配置|HMMsearch software configuration
@click.option('--hmmsearch-path',
              default='~/miniforge3/envs/resistify_v.1.3.0/bin/hmmsearch',
              show_default=True,
              help='hmmsearch程序路径|hmmsearch program path')
@click.option('-t', '--threads',
              default=12,
              show_default=True,
              type=int,
              help='线程数|Number of threads')
# HMMsearch运行参数|HMMsearch run parameters
@click.option('--evalue-cutoff',
              type=float,
              help='hmmsearch报告E-value阈值|hmmsearch reporting E-value threshold')
@click.option('--score-cutoff',
              type=float,
              help='hmmsearch报告分数阈值|hmmsearch reporting score threshold')
@click.option('--cut-tc',
              is_flag=True,
              help='使用模型的TC trusted cutoff|Use model TC trusted cutoff')
@click.option('--cut-ga',
              is_flag=True,
              help='使用模型的GA gathering cutoff|Use model GA gathering cutoff')
@click.option('--cut-nc',
              is_flag=True,
              help='使用模型的NC noise cutoff|Use model NC noise cutoff')
# 结果过滤参数|Result filtering parameters
@click.option('-e', '--evalue-threshold',
              type=float,
              help='Domain E-value阈值(保留小于等于该值的)|Domain E-value threshold (keep <= this value)')
@click.option('-s', '--score-threshold',
              type=float,
              help='Domain分数阈值(保留大于等于该值的)|Domain score threshold (keep >= this value)')
# 序列提取选项|Sequence extraction options
@click.option('--no-extract-proteins',
              is_flag=True,
              help='不提取蛋白序列|Do not extract protein sequences')
@click.option('--no-extract-domains',
              is_flag=True,
              help='不提取domain序列|Do not extract domain sequences')
# 输出格式选项|Output format options
@click.option('--no-csv',
              is_flag=True,
              help='不输出CSV文件|Do not output CSV file')
@click.option('--no-excel',
              is_flag=True,
              help='不输出Excel文件|Do not output Excel file')
def hmmsearch(input, protein_fasta, output_dir, output_prefix,
              hmmsearch_path, threads, evalue_cutoff, score_cutoff,
              cut_tc, cut_ga, cut_nc, evalue_threshold, score_threshold,
              no_extract_proteins, no_extract_domains, no_csv, no_excel):
    """
    HMMsearch分析工具：运行hmmsearch并处理结果|HMMsearch Analysis Tool: Run hmmsearch and process results

    模式1: 处理已有domtblout文件（不指定-p）|Mode 1: Process existing domtblout file (no -p)
    模式2: 运行hmmsearch并处理结果（指定-p）|Mode 2: Run hmmsearch and process results (with -p)

    示例|Examples: biopytools hmmsearch -i NB-ARC.hmm -p proteins.fa
    """

    # 延迟加载|Lazy loading
    hmmsearch_main = _lazy_import_hmmsearch_main()

    # 构建参数列表|Build argument list
    args = ['hmmsearch.py']

    # 主输入参数|Main input parameter
    if input:
        args.extend(['-i', input])

    # 蛋白序列参数|Protein fasta parameter
    if protein_fasta:
        args.extend(['-p', protein_fasta])

    # 输出配置|Output configuration
    if output_dir != './hmmsearch_output':
        args.extend(['-o', output_dir])

    if output_prefix != 'hmmsearch_results':
        args.extend(['--output-prefix', output_prefix])

    # HMMsearch软件配置|HMMsearch software configuration
    if hmmsearch_path != '~/miniforge3/envs/resistify_v.1.3.0/bin/hmmsearch':
        args.extend(['--hmmsearch-path', hmmsearch_path])

    if threads != 12:
        args.extend(['-t', str(threads)])

    # HMMsearch运行参数|HMMsearch run parameters
    if evalue_cutoff is not None:
        args.extend(['--evalue-cutoff', str(evalue_cutoff)])

    if score_cutoff is not None:
        args.extend(['--score-cutoff', str(score_cutoff)])

    if cut_tc:
        args.append('--cut-tc')

    if cut_ga:
        args.append('--cut-ga')

    if cut_nc:
        args.append('--cut-nc')

    # 结果过滤参数|Result filtering parameters
    if evalue_threshold is not None:
        args.extend(['-e', str(evalue_threshold)])

    if score_threshold is not None:
        args.extend(['-s', str(score_threshold)])

    # 布尔选项|Boolean options
    if no_extract_proteins:
        args.append('--no-extract-proteins')

    if no_extract_domains:
        args.append('--no-extract-domains')

    if no_csv:
        args.append('--no-csv')

    if no_excel:
        args.append('--no-excel')

    # 执行主程序|Execute main program
    original_argv = sys.argv
    sys.argv = args

    try:
        hmmsearch_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
