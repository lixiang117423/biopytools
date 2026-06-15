"""
Parabricks WGS分析命令|Parabricks WGS Analysis Command
"""

import click
import sys
import os


def _lazy_import_parabricks_main():
    """延迟加载parabricks主函数|Lazy load parabricks main function"""
    try:
        from ...parabricks.main import main as parabricks_main
        return parabricks_main
    except ImportError as e:
        click.echo(f"导入错误|Import error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_path_exists(path):
    """验证路径存在(仅在非帮助模式)|Validate path existence (only in non-help mode)"""
    if _is_help_request():
        return path
    if not os.path.exists(path):
        raise click.BadParameter(f"路径不存在|Path does not exist: {path}")
    return path


@click.command(
    short_help='Parabricks GPU加速WGS分析|Parabricks GPU-accelerated WGS analysis',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--input-dir', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_path_exists(value) if value else None,
              help='输入目录(FASTQ文件)|Input directory containing FASTQ files')
@click.option('--output-dir', '-o',
              required=True,
              type=click.Path(),
              help='输出目录|Output directory')
@click.option('--reference', '-r',
              required=True,
              callback=lambda ctx, param, value: _validate_path_exists(value) if value else None,
              help='参考基因组文件|Reference genome file')
@click.option('--workflow', '-w',
              type=click.Choice(['fq2bam', 'haplotypecaller', 'genotypegvcf', 'all']),
              default='all',
              show_default=True,
              help='工作流程|Workflow: fq2bam, haplotypecaller, genotypegvcf, or all')
@click.option('--threads', '-t',
              default=12,
              type=int,
              show_default=True,
              help='线程数|Number of threads')
@click.option('--parabricks-path',
              default='~/software/containers/parabricks.sif',
              type=str,
              show_default=True,
              help='Parabricks程序路径|Parabricks program path')
@click.option('--tmp-dir',
              type=click.Path(),
              help='临时目录|Temporary directory')
@click.option('--gvcf/--no-gvcf',
              default=True,
              show_default=True,
              help='输出GVCF格式|Output GVCF format')
@click.option('--joint-calling/--no-joint-calling',
              default=True,
              show_default=True,
              help='启用Joint Calling|Enable Joint Calling')
@click.option('--combined-output',
              default='combined.g.vcf',
              show_default=True,
              help='Joint Calling输出文件名|Joint Calling output filename')
@click.option('--min-confidence',
              default=30,
              type=int,
              show_default=True,
              help='最小置信度阈值|Minimum confidence threshold')
@click.option('--min-base-quality',
              default=20,
              type=int,
              show_default=True,
              help='最小碱基质量阈值|Minimum base quality threshold')
@click.option('--ploidy',
              default=2,
              type=int,
              show_default=True,
              help='倍性|Ploidy')
@click.option('--pcr-indel-model',
              default='CONSERVATIVE',
              type=str,
              show_default=True,
              help='PCR indel模型|PCR indel model')
@click.option('--read1-pattern',
              default='*_1.clean.fq.gz',
              type=str,
              show_default=True,
              help='R1文件模式|R1 file pattern')
@click.option('--read2-pattern',
              default='*_2.clean.fq.gz',
              type=str,
              show_default=True,
              help='R2文件模式|R2 file pattern')
def parabricks(input_dir, output_dir, reference, workflow, threads, parabricks_path, tmp_dir,
               gvcf, joint_calling, combined_output,
               min_confidence, min_base_quality, ploidy, pcr_indel_model,
               read1_pattern, read2_pattern):
    """
    Parabricks WGS批量分析工具|Parabricks WGS Batch Analysis Tool

    基于NVIDIA Parabricks进行GPU加速的全基因组测序分析
    Perform GPU-accelerated whole genome sequencing analysis using NVIDIA Parabricks

    示例|Examples: biopytools parabricks -i /data/fastq -o /results -r /ref/genome.fa
    """

    # 延迟加载|Lazy load
    parabricks_main = _lazy_import_parabricks_main()

    # 构建主函数参数列表|Build argument list for main function
    args = ['parabricks.py']

    # 必需参数|Required parameters
    args.extend(['-i', input_dir])
    args.extend(['-o', output_dir])
    args.extend(['-r', reference])

    # 可选参数(仅在非默认时添加)|Optional parameters (add only when non-default)
    if workflow != 'all':
        args.extend(['-w', workflow])

    if threads != 88:
        args.extend(['-t', str(threads)])

    if parabricks_path != '~/software/containers/parabricks.sif':
        args.extend(['--parabricks-path', parabricks_path])

    if tmp_dir:
        args.extend(['--tmp-dir', tmp_dir])

    # GVCF和Joint Calling参数|GVCF and Joint Calling parameters
    if not gvcf:
        args.append('--no-gvcf')
    if not joint_calling:
        args.append('--no-joint-calling')
    if combined_output != 'combined.g.vcf':
        args.extend(['--combined-output', combined_output])

    # 质控参数|Quality control parameters
    if min_confidence != 30:
        args.extend(['--min-confidence', str(min_confidence)])
    if min_base_quality != 20:
        args.extend(['--min-base-quality', str(min_base_quality)])
    if ploidy != 2:
        args.extend(['--ploidy', str(ploidy)])
    if pcr_indel_model != 'CONSERVATIVE':
        args.extend(['--pcr-indel-model', pcr_indel_model])

    # 文件模式参数|File pattern parameters
    if read1_pattern != '*_1.clean.fq.gz':
        args.extend(['--read1-pattern', read1_pattern])
    if read2_pattern != '*_2.clean.fq.gz':
        args.extend(['--read2-pattern', read2_pattern])

    # 保存并恢复sys.argv|Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        # 调用主函数|Call main function
        parabricks_main()
    except SystemExit as e:
        # 处理系统退出|Handle system exit
        if e.code != 0:
            sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n用户中断|User interrupted", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"运行时错误|Runtime error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
