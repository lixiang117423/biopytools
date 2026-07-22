"""phylo-trim 整合命令|phylo-trim Integration Command"""

import click
import sys
import os


def _lazy_import_runner():
    """延迟加载 PhyloTrimRunner|Lazy load PhyloTrimRunner"""
    try:
        from ...phylo_trim.main import PhyloTrimRunner
        return PhyloTrimRunner
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    return any(arg in {'-h', '--help'} for arg in sys.argv)


def _validate_file_exists(file_path):
    if not _is_help_request() and file_path and not os.path.exists(file_path):
        raise click.BadParameter(f"文件不存在|File does not exist: {file_path}")
    return file_path


@click.command(
    short_help='整合 mafft-fasttree+trimal,出前后两棵树|Integrate mafft-fasttree+trimal, before/after trees',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--input', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              type=click.Path(exists=True, file_okay=True, dir_okay=False, resolve_path=True),
              help='输入原始序列 FASTA|Input raw sequence FASTA')
@click.option('--output-dir', '-o',
              default='./phylo_trim_output',
              show_default=True,
              type=click.Path(),
              help='输出目录|Output directory')
@click.option('--skip-trimal',
              is_flag=True,
              help='跳过 trimal,只出 trimal 前树|Skip trimal, before-tree only')
@click.option('--trimal-method',
              default='automated1',
              show_default=True,
              type=click.Choice(['automated1', 'gappyout', 'strict', 'strictplus', 'gt', 'cons']),
              help='trimal 方法|trimal method')
@click.option('--gt-threshold',
              default=0.9, show_default=True, type=float,
              help='trimal gap 阈值[0,1]|trimal gap threshold [0,1]')
@click.option('--cons-threshold',
              default=80, show_default=True, type=int,
              help='trimal 保守度[0,100]|trimal conservation [0,100]')
@click.option('--trimal-format',
              default='keep', show_default=True,
              type=click.Choice(['keep', 'fasta', 'phylip', 'phylip_paml', 'clustal', 'nexus', 'nbrf', 'mega']),
              help='trimal 输出格式|trimal output format (keep=沿用输入|input format)')
@click.option('--seq-type',
              type=click.Choice(['protein', 'nucleotide']),
              help='序列类型(默认自动检测)|Sequence type (auto-detect by default)')
@click.option('--threads', '-t',
              default=88, show_default=True, type=int,
              help='MAFFT 线程数|MAFFT threads')
@click.option('--mafft-params',
              default='--auto', show_default=True,
              help='MAFFT 额外参数|Additional MAFFT parameters')
@click.option('--fasttree-params',
              default='', show_default=True,
              help='FastTree 额外参数|Additional FastTree parameters')
@click.option('--mafft-path',
              default='mafft', show_default=True,
              help='MAFFT 路径|MAFFT path')
@click.option('--fasttree-path',
              default='fasttree', show_default=True,
              help='FastTree 路径|FastTree path')
@click.option('--sample-name',
              help='输出文件名前缀(默认输入 basename)|Output filename prefix (default: input basename)')
@click.option('--log-file',
              type=click.Path(),
              help='日志文件路径|Log file path')
@click.option('--verbose', '-v',
              is_flag=True,
              help='详细输出|Verbose output')
def phylo_trim(input, output_dir, skip_trimal, trimal_method, gt_threshold,
               cons_threshold, trimal_format, seq_type, threads, mafft_params,
               fasttree_params, mafft_path, fasttree_path, sample_name, log_file, verbose):
    """整合 mafft-fasttree + trimal|Integrate mafft-fasttree + trimal

    自动输出 trimal 前后两棵 FastTree 系统发育树|Output before/after-trimal FastTree trees

    示例|Examples: biopytools phylo-trim -i seqs.fa -o out/
    """

    try:
        PhyloTrimRunner = _lazy_import_runner()

        runner = PhyloTrimRunner(
            input_file=str(input),
            output_dir=str(output_dir),
            skip_trimal=skip_trimal,
            trimal_method=trimal_method,
            gt_threshold=gt_threshold,
            cons_threshold=cons_threshold,
            trimal_format=trimal_format,
            seq_type=seq_type,
            threads=threads,
            mafft_params=mafft_params,
            fasttree_params=fasttree_params,
            mafft_path=mafft_path,
            fasttree_path=fasttree_path,
            sample_name=sample_name,
            log_file=str(log_file) if log_file else None,
            verbose=verbose,
        )

        success = runner.run()

        if success:
            click.echo("phylo-trim 完成|phylo-trim completed successfully!")
        else:
            click.echo("phylo-trim 失败|phylo-trim failed!", err=True)
            sys.exit(1)

    except KeyboardInterrupt:
        click.echo("\n用户中断|Interrupted by user", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"发生错误|Error occurred: {e}", err=True)
        sys.exit(1)
