"""Phytophthora效应子鉴定命令|Phytophthora Effector Identification Command"""

import click
import sys
import os
import logging


def _is_help_request():
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_path(path):
    if not _is_help_request() and path and not os.path.exists(path):
        raise click.BadParameter(f"路径不存在|Path does not exist: {path}")
    return path


@click.group(
    short_help='Phytophthora效应子鉴定(RxLR+CRN+NLP+...)|Phytophthora effector identification',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120),
    invoke_without_command=True
)
@click.pass_context
def phyto_effector(ctx):
    """
    Phytophthora效应子鉴定(RxLR+CRN+NLP+...)|Phytophthora Effector Identification

    鉴定类型: RxLR, CRN, NLP, Protease, SCP, elicitin, YxSL
    流程: SignalP -> 各类型HMM搜索 + 基序注释
    """
    if ctx.invoked_subcommand is None:
        ctx.invoke(run)


@phyto_effector.command()
@click.option('-i', '--input', required=True,
              callback=lambda ctx, param, value: _validate_path(value) if value else None,
              help='输入FASTA文件或目录|Input FASTA file or directory')
@click.option('-o', '--output-dir', default='./phyto_effector_output',
              show_default=True, help='输出目录|Output directory')
@click.option('--skip-signalp', is_flag=True,
              help='跳过SignalP预测|Skip SignalP prediction')
@click.option('--signalp-path', default='~/miniforge3/envs/signalp6/bin/signalp6',
              show_default=True, help='SignalP程序路径|SignalP program path')
@click.option('--organism', default='eukarya', show_default=True,
              type=click.Choice(['eukarya', 'other']),
              help='生物类型|Organism type')
@click.option('--signalp-mode', default='slow-sequential', show_default=True,
              type=click.Choice(['fast', 'slow', 'slow-sequential']),
              help='SignalP运行模式|SignalP run mode')
@click.option('--signalp-version', default='both', show_default=True,
              type=click.Choice(['3', '6', 'both']),
              help='SignalP版本|SignalP version (3/6/both)')
@click.option('--signalp3-path',
              default='~/miniforge3/envs/signalp_v.3.0b/bin/signalp',
              show_default=True, help='SignalP 3.0程序路径|SignalP 3.0 program path')
@click.option('--signalp3-sprob-threshold', default=0.9, show_default=True, type=float,
              help='SignalP 3.0 HMM Sprob阈值|SignalP 3.0 HMM Sprob threshold')
@click.option('--hmmsearch-path', default='~/miniforge3/envs/resistify_v.1.3.0/bin/hmmsearch',
              help='hmmsearch程序路径|hmmsearch program path')
@click.option('--rxlr-hmm', default=None,
              help='RxLR HMM文件(默认内置)|RxLR HMM file (default: bundled)')
@click.option('--use-wy-domain', is_flag=True,
              help='同时搜索WY结构域|Also search WY domain')
@click.option('--rxlr-wy-hmm', default=None,
              help='WY HMM文件(默认内置)|WY HMM file (default: bundled)')
@click.option('--crn-hmm', default=None,
              help='CRN HMM文件(默认内置)|CRN HMM file (default: bundled)')
@click.option('--nlp-hmm', default=None,
              help='NLP HMM文件(默认内置)|NLP HMM file (default: bundled)')
@click.option('--protease-hmm', default=None,
              help='Protease HMM文件(默认内置)|Protease HMM file (default: bundled)')
@click.option('--scp-hmm', default=None,
              help='SCP HMM文件(默认内置)|SCP HMM file (default: bundled)')
@click.option('--elicitin-hmm', default=None,
              help='elicitin HMM文件(默认内置)|elicitin HMM file (default: bundled)')
@click.option('--yxsl-hmm', default=None,
              help='YxSL HMM文件(默认内置)|YxSL HMM file (default: bundled)')
@click.option('-e', '--evalue', default=1e-5, type=float,
              help='E-value阈值(已弃用)|E-value threshold (deprecated)')
@click.option('--score-threshold', default=0.0, type=float,
              help='HMM score阈值|HMM score threshold')
@click.option('-t', '--threads', default=12, type=int,
              help='线程数|Number of threads')
def run(input, output_dir, skip_signalp, signalp_path, organism, signalp_mode,
        signalp_version, signalp3_path, signalp3_sprob_threshold,
        hmmsearch_path, rxlr_hmm, use_wy_domain, rxlr_wy_hmm,
        crn_hmm, nlp_hmm, protease_hmm, scp_hmm, elicitin_hmm, yxsl_hmm,
        evalue, score_threshold, threads):
    """
    运行效应子鉴定流程|Run effector identification pipeline

    示例|Examples: biopytools phyto-effector -i proteins.fa -o effector_out
    """
    try:
        from ...phyto_effector.main import main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)

    args = ['phyto_effector.py']
    args.extend(['-i', input])
    if output_dir != './phyto_effector_output':
        args.extend(['-o', output_dir])
    if skip_signalp:
        args.append('--skip-signalp')
    if signalp_path != '~/miniforge3/envs/signalp6/bin/signalp6':
        args.extend(['--signalp-path', signalp_path])
    if organism != 'eukarya':
        args.extend(['--organism', organism])
    if signalp_mode != 'fast':
        args.extend(['--signalp-mode', signalp_mode])
    if signalp_version != 'both':
        args.extend(['--signalp-version', signalp_version])
    if signalp3_path != '~/miniforge3/envs/signalp_v.3.0b/bin/signalp':
        args.extend(['--signalp3-path', signalp3_path])
    if signalp3_sprob_threshold != 0.9:
        args.extend(['--signalp3-sprob-threshold', str(signalp3_sprob_threshold)])
    if hmmsearch_path != '~/miniforge3/envs/resistify_v.1.3.0/bin/hmmsearch':
        args.extend(['--hmmsearch-path', hmmsearch_path])
    if rxlr_hmm:
        args.extend(['--rxlr-hmm', rxlr_hmm])
    if use_wy_domain:
        args.append('--use-wy-domain')
    if rxlr_wy_hmm:
        args.extend(['--rxlr-wy-hmm', rxlr_wy_hmm])
    if crn_hmm:
        args.extend(['--crn-hmm', crn_hmm])
    if nlp_hmm:
        args.extend(['--nlp-hmm', nlp_hmm])
    if protease_hmm:
        args.extend(['--protease-hmm', protease_hmm])
    if scp_hmm:
        args.extend(['--scp-hmm', scp_hmm])
    if elicitin_hmm:
        args.extend(['--elicitin-hmm', elicitin_hmm])
    if yxsl_hmm:
        args.extend(['--yxsl-hmm', yxsl_hmm])
    if evalue != 1e-5:
        args.extend(['-e', str(evalue)])
    if score_threshold != 0.0:
        args.extend(['--score-threshold', str(score_threshold)])
    if threads != 12:
        args.extend(['-t', str(threads)])

    original_argv = sys.argv
    sys.argv = args
    try:
        main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv


@phyto_effector.command()
@click.option('-i', '--input-dir', required=True,
              callback=lambda ctx, param, value: _validate_path(value) if value else None,
              help='多个样品输出目录的父目录|Parent dir containing sample output dirs')
@click.option('-o', '--output-dir', default='./phyto_effector_merged',
              show_default=True, help='合并结果输出目录|Merged output directory')
def merge(input_dir, output_dir):
    """
    合并多个样品的效应子鉴定结果|Merge effector results from multiple samples

    将多个单独运行的样品结果合并为汇总文件。

    示例|Examples: biopytools phyto-effector merge -i ./sample_outputs -o ./merged
    """
    logger = logging.getLogger('phyto_effector_merge')
    logger.addHandler(logging.StreamHandler())
    logger.setLevel(logging.INFO)

    try:
        from ...phyto_effector.utils import merge_candidate_files
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)

    effector_types = ['rxlr', 'crn', 'nlp', 'protease', 'scp', 'elicitin', 'yxsl']

    sample_dirs = sorted([
        os.path.join(input_dir, d)
        for d in os.listdir(input_dir)
        if os.path.isdir(os.path.join(input_dir, d))
        and not d.startswith('.')
        and d not in ('99_logs', '00_pipeline_info')
    ])

    if not sample_dirs:
        click.echo(f"未找到样品目录|No sample directories found in: {input_dir}", err=True)
        sys.exit(1)

    sample_names = [os.path.basename(d) for d in sample_dirs]
    click.echo(f"找到|Found {len(sample_dirs)}个样品|samples: {', '.join(sample_names)}")
    click.echo(f"输出目录|Output directory: {output_dir}")
    click.echo("")

    os.makedirs(output_dir, exist_ok=True)

    for etype in effector_types:
        result = merge_candidate_files(sample_dirs, etype, output_dir, logger)
        if result:
            import pandas as pd
            df = pd.read_csv(result, sep='\t')
            samples = df['Sample'].nunique() if 'Sample' in df.columns else 0
            click.echo(f"  {etype:>10s}: {len(df):>6d}条记录|records ({samples}个样品|samples)")
        else:
            click.echo(f"  {etype:>10s}: 无数据|no data")

    click.echo(f"\n合并完成|Merge completed: {output_dir}")
