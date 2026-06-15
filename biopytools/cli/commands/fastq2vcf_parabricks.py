"""
Fastq到VCF (Parabricks)命令|Fastq to VCF (Parabricks) Command
"""

import click
import sys
import os


def _lazy_import_fastq2vcf_parabricks_main():
    """延迟加载fastq2vcf_parabricks主函数|Lazy load fastq2vcf_parabricks main function"""
    try:
        from ...fastq2vcf_parabricks.main import main as fastq2vcf_parabricks_main
        return fastq2vcf_parabricks_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_file_exists(file_path):
    """验证文件存在性(仅在非帮助模式下)|Validate file existence (only in non-help mode)"""
    if not _is_help_request() and file_path and not os.path.exists(file_path):
        raise click.BadParameter(f"文件不存在|File does not exist: {file_path}")
    return file_path


def _validate_directory_exists(dir_path):
    """验证目录存在性(仅在非帮助模式下)|Validate directory existence (only in non-help mode)"""
    if not _is_help_request() and dir_path and not os.path.exists(dir_path):
        raise click.BadParameter(f"目录不存在|Directory does not exist: {dir_path}")
    return dir_path


@click.command(short_help="Fastq到VCF (Parabricks)全流程|Fastq to VCF (Parabricks) Pipeline",
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--version', '-V',
              is_flag=True,
              help='显示版本信息|Show version information')
@click.option('--input', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_directory_exists(value) if value else None,
              type=click.Path(exists=True, file_okay=False, dir_okay=True, resolve_path=True),
              help='原始FASTQ文件目录路径|Raw FASTQ files directory path')
@click.option('--genome', '-g',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              type=click.Path(exists=True, file_okay=True, dir_okay=False, resolve_path=True),
              help='参考基因组文件路径|Reference genome file path')
@click.option('--project-base', '-p',
              required=True,
              callback=lambda ctx, param, value: _validate_directory_exists(value) if value else None,
              type=click.Path(exists=True, file_okay=False, dir_okay=True, resolve_path=True),
              help='项目根目录路径|Project base directory path')
@click.option('--clean-fastq-dir',
              type=click.Path(file_okay=False, dir_okay=True, resolve_path=True),
              help='清洁FASTQ文件目录路径|Clean FASTQ files directory path')
@click.option('--mapping-dir',
              type=click.Path(file_okay=False, dir_okay=True, resolve_path=True),
              help='比对结果目录路径|Mapping results directory path')
@click.option('--joint-dir',
              type=click.Path(file_okay=False, dir_okay=True, resolve_path=True),
              help='联合检测结果目录路径|Joint calling results directory path')
@click.option('--filter-dir',
              type=click.Path(file_okay=False, dir_okay=True, resolve_path=True),
              help='过滤结果目录路径|Filtering results directory path')
@click.option('--output-dir', '-o',
              type=click.Path(file_okay=False, dir_okay=True, resolve_path=True),
              help='输出目录路径|Output directory path')
@click.option('--threads-mapping',
              type=int,
              default=12,
              show_default=True,
              help='比对线程数|Number of mapping threads')
@click.option('--threads-gtx',
              type=int,
              default=12,
              show_default=True,
              help='GTX线程数|Number of GTX threads')
@click.option('--threads-filter',
              type=int,
              default=12,
              show_default=True,
              help='过滤线程数|Number of filtering threads')
@click.option('--snp-min-dp',
              type=int,
              default=5,
              show_default=True,
              help='SNP最小深度|SNP minimum depth')
@click.option('--snp-min-qual',
              type=int,
              default=30,
              show_default=True,
              help='SNP最小质量|SNP minimum quality')
@click.option('--indel-min-dp',
              type=int,
              default=5,
              show_default=True,
              help='InDel最小深度|InDel minimum depth')
@click.option('--indel-min-qual',
              type=int,
              default=30,
              show_default=True,
              help='InDel最小质量|InDel minimum quality')
@click.option('--gatk-threshold',
              type=int,
              default=4,
              show_default=True,
              help='GATK模式样本数阈值|GATK sample count threshold')
@click.option('--gtx-single-threshold',
              type=int,
              default=200,
              show_default=True,
              help='GTX单机模式样本数阈值|GTX single machine sample count threshold')
@click.option('--gtx-window-size',
              type=int,
              default=20000000,
              show_default=True,
              help='GTX分块窗口大小|GTX chunk window size in bp')
@click.option('--gtx-bin',
              default='~/software/gtx/bin/gtx',
              show_default=True,
              type=click.Path(exists=True, file_okay=True, dir_okay=False, resolve_path=True),
              help='GTX可执行文件路径|GTX executable path')
@click.option('--step', '-s',
              type=click.Choice(['1', '2', '3', '4', '5']),
              help='只运行指定步骤|Run only specified step')
@click.option('--no-checkpoint',
              is_flag=True,
              help='禁用断点续传|Disable checkpoint resume')
@click.option('--dry-run',
              is_flag=True,
              help='测试模式|Test mode')
@click.option('--verbose', '-v',
              count=True,
              help='详细输出模式|Verbose output mode')
@click.option('--quiet',
              is_flag=True,
              help='静默模式|Quiet mode')
@click.option('--log-level',
              type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL']),
              default='INFO',
              show_default=True,
              help='日志级别|Log level')
@click.option('--log-file',
              type=click.Path(),
              help='日志文件路径|Log file path')
@click.option('--force', '-f',
              is_flag=True,
              help='强制覆盖已存在的结果|Force overwrite existing results')
@click.option('--skip-qc',
              is_flag=True,
              help='跳过质控步骤|Skip QC step')
@click.option('--skip-mapping',
              is_flag=True,
              help='跳过比对步骤|Skip mapping step')
@click.option('--read1-pattern-fastp',
              default="_1.fq.gz",
              show_default=True,
              help='质控R1文件匹配模式|QC R1 file pattern')
@click.option('--read2-pattern-fastp',
              default="_2.fq.gz",
              show_default=True,
              help='质控R2文件匹配模式|QC R2 file pattern')
def fastq2vcf_parabricks(version, input, genome, project_base, clean_fastq_dir,
                         mapping_dir, joint_dir, filter_dir, output_dir,
                         threads_mapping, threads_gtx, threads_filter,
                         snp_min_dp, snp_min_qual, indel_min_dp, indel_min_qual,
                         gatk_threshold, gtx_single_threshold, gtx_window_size,
                         gtx_bin, step, no_checkpoint, dry_run, verbose, quiet,
                         log_level, log_file, force,
                         skip_qc, skip_mapping, read1_pattern_fastp, read2_pattern_fastp):
    """
    基于Parabricks的重测序全基因组变异检测全流程|Whole genome resequencing variant detection pipeline based on Parabricks

    示例|Examples: biopytools fastq2vcf-parabricks -i /path/to/raw_fastq -g /path/to/genome.fa -p /path/to/project
    """

    # 处理版本信息|Handle version information
    if version:
        from ...fastq2vcf_parabricks.main import VERSION
        click.echo(f"Fastq to VCF (Parabricks) Pipeline Version: {VERSION}")
        sys.exit(0)

    # 延迟加载|Lazy loading: import only when actually called
    fastq2vcf_parabricks_main = _lazy_import_fastq2vcf_parabricks_main()

    # 构建参数列表传递给原始main函数|Build argument list for original main function
    args = ['fastq2vcf_parabricks.py']

    # 必需参数|Required parameters
    args.extend(['--input', input])
    args.extend(['-g', genome])
    args.extend(['--project-base', project_base])

    # 可选参数（只在非默认值时添加）|Optional parameters (add only when non-default)
    if clean_fastq_dir:
        args.extend(['--clean-fastq-dir', clean_fastq_dir])

    if mapping_dir:
        args.extend(['--mapping-dir', mapping_dir])

    if joint_dir:
        args.extend(['--joint-dir', joint_dir])

    if filter_dir:
        args.extend(['--filter-dir', filter_dir])

    if output_dir:
        args.extend(['--output-dir', output_dir])

    if threads_mapping != 88:
        args.extend(['--threads-mapping', str(threads_mapping)])

    if threads_gtx != 88:
        args.extend(['--threads-gtx', str(threads_gtx)])

    if threads_filter != 88:
        args.extend(['--threads-filter', str(threads_filter)])

    if snp_min_dp != 5:
        args.extend(['--snp-min-dp', str(snp_min_dp)])

    if snp_min_qual != 30:
        args.extend(['--snp-min-qual', str(snp_min_qual)])

    if indel_min_dp != 5:
        args.extend(['--indel-min-dp', str(indel_min_dp)])

    if indel_min_qual != 30:
        args.extend(['--indel-min-qual', str(indel_min_qual)])

    if gatk_threshold != 4:
        args.extend(['--gatk-threshold', str(gatk_threshold)])

    if gtx_single_threshold != 200:
        args.extend(['--gtx-single-threshold', str(gtx_single_threshold)])

    if gtx_window_size != 20000000:
        args.extend(['--gtx-window-size', str(gtx_window_size)])

    if gtx_bin != '~/software/gtx/bin/gtx':
        args.extend(['--gtx-bin', gtx_bin])

    # 步骤控制|Step control
    if step:
        args.extend(['--step', step])

    # 处理选项（布尔标志）|Processing options (boolean flags)
    if no_checkpoint:
        args.append('--no-checkpoint')

    if dry_run:
        args.append('--dry-run')

    # 处理verbose (count)|Handle verbose (count)
    if verbose > 0:
        args.extend(['--verbose'] * verbose)

    if quiet:
        args.append('--quiet')

    if log_level != 'INFO':
        args.extend(['--log-level', log_level])

    if log_file:
        args.extend(['--log-file', log_file])

    if force:
        args.append('--force')

    if skip_qc:
        args.append('--skip-qc')

    if skip_mapping:
        args.append('--skip-mapping')

    if read1_pattern_fastp != "_1.fq.gz":
        args.extend(['--read1-pattern-fastp', read1_pattern_fastp])

    if read2_pattern_fastp != "_2.fq.gz":
        args.extend(['--read2-pattern-fastp', read2_pattern_fastp])

    # 保存并恢复sys.argv|Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        # 调用原始的main函数|Call original main function
        fastq2vcf_parabricks_main()
    except SystemExit as e:
        # 处理程序正常退出|Handle normal program exit
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
