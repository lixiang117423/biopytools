"""
ALLHiC Pipeline CLI Wrapper|ALLHiC流程CLI包装器
"""

import click
import sys
import os
from pathlib import Path


def _lazy_import_allhic_main():
    """延迟加载ALLHiC主函数|Lazy load ALLHiC main function"""
    try:
        from ...allhic.main import main as allhic_main
        return allhic_main
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
    short_help='ALLHiC基因组scaffolding流程|ALLHiC genome scaffolding pipeline',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--reference', '-r',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              type=click.Path(exists=True),
              help='参考基因组文件|Reference genome file')
@click.option('--read1', '-1',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              type=click.Path(exists=True),
              help='Hi-C读段1文件|Hi-C read 1 file')
@click.option('--read2', '-2',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              type=click.Path(exists=True),
              help='Hi-C读段2文件|Hi-C read 2 file')
@click.option('--chr-num', '-k',
              required=True,
              type=int,
              help='染色体数目|Number of chromosomes')
@click.option('--enzyme', '-e',
              default='GATC',
              show_default=True,
              help='酶切位点|Enzyme motif')
@click.option('--threads', '-t',
              type=int,
              default=12,
              show_default=True,
              help='线程数|Number of threads')
@click.option('--workdir', '-w',
              default='./allhic_output',
              show_default=True,
              help='工作目录|Working directory')
@click.option('--mapq-step1',
              type=int,
              default=1,
              show_default=True,
              help='步骤1比对质量阈值|Step 1 MapQ threshold')
@click.option('--bin-size',
              default='500k',
              show_default=True,
              help='Bin大小|Bin size')
@click.option('--min-bin-size',
              default='50k',
              show_default=True,
              help='最小Bin大小|Minimum bin size')
@click.option('--skip-mapping',
              is_flag=True,
              help='跳过步骤1:比对|Skip Step 1: Mapping')
@click.option('--skip-allele',
              is_flag=True,
              help='跳过步骤1.5:等位基因检测|Skip Step 1.5: Allele Detection')
@click.option('--skip-prune',
              is_flag=True,
              help='跳过步骤2:修剪|Skip Step 2: Pruning')
@click.option('--skip-partition',
              is_flag=True,
              help='跳过步骤3:分区|Skip Step 3: Partition')
@click.option('--skip-extract',
              is_flag=True,
              help='跳过步骤3.5:提取矩阵|Skip Step 3.5: Extract Matrix')
@click.option('--skip-rescue',
              is_flag=True,
              help='跳过步骤4:救援|Skip Step 4: Rescue')
@click.option('--skip-optimize',
              is_flag=True,
              help='跳过步骤5:优化|Skip Step 5: Optimization')
@click.option('--skip-build',
              is_flag=True,
              help='跳过步骤6:构建FASTA|Skip Step 6: Build FASTA')
@click.option('--skip-plot',
              is_flag=True,
              help='跳过步骤7:绘制热图|Skip Step 7: Plot Heatmap')
@click.option('--skip-asmkit',
              is_flag=True,
              help='跳过步骤8:JBAT生成|Skip Step 8: JBAT Generation')
@click.option('--diagnose',
              is_flag=True,
              help='诊断模式|Diagnostic mode')
@click.option('--verbose', '-v',
              is_flag=True,
              help='详细输出|Verbose output')
def allhic(reference, read1, read2, chr_num, enzyme, threads, workdir,
           mapq_step1, bin_size, min_bin_size, skip_mapping, skip_allele,
           skip_prune, skip_partition, skip_extract, skip_rescue,
           skip_optimize, skip_build, skip_plot, skip_asmkit,
           diagnose, verbose):
    """
    ALLHiC基因组scaffolding流程|ALLHiC Genome Scaffolding Pipeline

    Hi-C辅助基因组装配流程|Hi-C assisted genome assembly pipeline

    示例|Examples: biopytools allhic -r genome.fa -1 hic_R1.fq.gz -2 hic_R2.fq.gz -k 12
    """

    # 延迟加载|Lazy loading
    allhic_main = _lazy_import_allhic_main()

    # 构建参数列表|Build argument list
    args = ['allhic.py']

    # 必需参数|Required parameters
    args.extend(['-r', reference])
    args.extend(['-1', read1])
    args.extend(['-2', read2])
    args.extend(['-k', str(chr_num)])

    # 可选参数|Optional parameters
    if enzyme != 'GATC':
        args.extend(['-e', enzyme])

    if threads != 88:
        args.extend(['-t', str(threads)])

    if workdir != './allhic_output':
        args.extend(['-w', workdir])

    if mapq_step1 != 1:
        args.extend(['--mapq-step1', str(mapq_step1)])

    if bin_size != '500k':
        args.extend(['--bin-size', bin_size])

    if min_bin_size != '50k':
        args.extend(['--min-bin-size', min_bin_size])

    # 布尔选项|Boolean options
    if skip_mapping:
        args.append('--skip-mapping')

    if skip_allele:
        args.append('--skip-allele')

    if skip_prune:
        args.append('--skip-prune')

    if skip_partition:
        args.append('--skip-partition')

    if skip_extract:
        args.append('--skip-extract')

    if skip_rescue:
        args.append('--skip-rescue')

    if skip_optimize:
        args.append('--skip-optimize')

    if skip_build:
        args.append('--skip-build')

    if skip_plot:
        args.append('--skip-plot')

    if skip_asmkit:
        args.append('--skip-asmkit')

    if diagnose:
        args.append('--diagnose')

    if verbose:
        args.append('--verbose')

    # 执行主程序|Execute main program
    original_argv = sys.argv
    sys.argv = args

    try:
        allhic_main()
    except SystemExit as e:
        if e.code != 0:
            sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n用户中断操作|Operation interrupted by user", err=True)
        sys.exit(130)
    except Exception as e:
        click.echo(f"执行失败|Execution failed: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
