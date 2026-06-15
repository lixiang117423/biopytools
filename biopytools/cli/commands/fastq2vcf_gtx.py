"""
Fastq to VCF (GTX) Command|基于GTX的Fastq转VCF命令
"""

import click
import sys
import os


def _lazy_import_fastq2vcf_gtx_main():
    """延迟加载fastq2vcf_gtx主函数|Lazy load fastq2vcf_gtx main function"""
    try:
        from ...fastq2vcf_gtx.main import main as fastq2vcf_gtx_main
        return fastq2vcf_gtx_main
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


@click.command(short_help="FastqVCF (GTX)|Fastq to VCF (GTX) Command",
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--input', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_directory_exists(value) if value else None,
              type=click.Path(exists=True, file_okay=False, dir_okay=True, resolve_path=True),
              help='输入目录|Raw FASTQ files directory path')
@click.option('--genome', '-g',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              type=click.Path(exists=True, file_okay=True, dir_okay=False, resolve_path=True),
              help='参考基因组|Reference genome file path')
@click.option('--output-dir', '-o',
              default='.',
              show_default=True,
              type=click.Path(file_okay=False, dir_okay=True, resolve_path=True),
              help='输出目录|Output directory path')
@click.option('--clean-fastq-dir',
              type=click.Path(file_okay=False, dir_okay=True, resolve_path=True),
              help='清洗后FASTQ目录|Clean FASTQ files directory path')
@click.option('--mapping-dir',
              type=click.Path(file_okay=False, dir_okay=True, resolve_path=True),
              help='比对结果目录|Mapping results directory path')
@click.option('--gvcf-dir',
              type=click.Path(file_okay=False, dir_okay=True, resolve_path=True),
              help='gVCF文件目录|gVCF files directory path')
@click.option('--bam-dir',
              type=click.Path(file_okay=False, dir_okay=True, resolve_path=True),
              help='BAM文件目录|BAM files directory path')
@click.option('--joint-dir',
              type=click.Path(file_okay=False, dir_okay=True, resolve_path=True),
              help='联合检测目录|Joint calling results directory path')
@click.option('--filter-dir',
              type=click.Path(file_okay=False, dir_okay=True, resolve_path=True),
              help='过滤结果目录|Filtering results directory path')
@click.option('--threads', '-t',
              type=int,
              default=12,
              show_default=True,
              help='线程数|Number of threads')
@click.option('--min-depth', '--snp-min-dp',
              type=int,
              default=5,
              show_default=True,
              help='最小测序深度|Minimum depth for SNP/InDel')
@click.option('--quality', '-q', '--min-qual', '--snp-min-qual',
              type=int,
              default=30,
              show_default=True,
              help='最小质量值|Minimum quality for SNP/InDel')
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
@click.option('--gtx-single-threshold',
              type=int,
              default=200,
              show_default=True,
              help='GTX单机样本阈值|GTX single machine sample count threshold')
@click.option('--gtx-window-size',
              type=int,
              default=20000000,
              show_default=True,
              help='GTX窗口大小|GTX chunk window size in bp')
@click.option('--gtx-bin',
              default='~/software/gtx/bin/gtx',
              show_default=True,
              type=click.Path(file_okay=True, dir_okay=False),
              help='GTX可执行文件路径|GTX executable path')
@click.option('--use-gtx-wgs',
              is_flag=True,
              default=True,
              show_default=True,
              help='使用GTX WGS模式|Use GTX WGS')
@click.option('--gtx-pcr-indel-model',
              default='CONSERVATIVE',
              show_default=True,
              help='GTX PCR InDel模型|GTX PCR InDel model')
@click.option('--gtx-min-confidence',
              type=int,
              default=30,
              show_default=True,
              help='GTX最小置信度|GTX minimum confidence')
@click.option('--gtx-min-base-qual',
              type=int,
              default=20,
              show_default=True,
              help='GTX最小碱基质量|GTX minimum base quality')
@click.option('--gtx-ploidy',
              type=int,
              default=2,
              show_default=True,
              help='GTX倍性|GTX ploidy')
@click.option('--gtx-cmd-gen-script',
              default="${HOME}/software/scripts/51.生成GTX按染色体合并gVCF的脚本.sh",
              show_default=True,
              help='GTX命令生成脚本|GTX command generation script path')
@click.option('--step', '-s',
              type=click.Choice(['1', '2', '3', '4', '5']),
              help='运行指定步骤|Run only specified step')
@click.option('--no-checkpoint',
              is_flag=True,
              help='禁用检查点恢复|Disable checkpoint resume')
@click.option('--dry-run',
              is_flag=True,
              help='试运行|Dry run')
@click.option('--force', '-f',
              is_flag=True,
              help='强制覆盖|Force overwrite')
@click.option('--keep-intermediate',
              is_flag=True,
              help='保留中间文件|Keep intermediate files')
@click.option('--verbose', '-v',
              count=True,
              help='详细输出|Verbose output')
@click.option('--quiet',
              is_flag=True,
              help='仅输出错误|Quiet mode')
@click.option('--log-file',
              type=click.Path(),
              help='日志文件|Log file path')
@click.option('--log-level',
              type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL']),
              default='INFO',
              show_default=True,
              help='日志级别|Log level')
@click.option('--skip-qc',
              is_flag=True,
              help='跳过质控|Skip QC')
@click.option('--skip-mapping',
              is_flag=True,
              help='跳过比对|Skip mapping')
@click.option('--read1-pattern-fastp',
              default="_1.fq.gz",
              show_default=True,
              help='R1文件模式|QC R1 file pattern')
@click.option('--read2-pattern-fastp',
              default="_2.fq.gz",
              show_default=True,
              help='R2文件模式|QC R2 file pattern')
def fastq2vcf_gtx(input, genome, output_dir,
                         clean_fastq_dir, mapping_dir, gvcf_dir, bam_dir, joint_dir, filter_dir,
                         threads,
                         min_depth, quality, indel_min_dp, indel_min_qual,
                         gtx_single_threshold, gtx_window_size,
                         gtx_bin, use_gtx_wgs, gtx_pcr_indel_model,
                         gtx_min_confidence, gtx_min_base_qual, gtx_ploidy, gtx_cmd_gen_script,
                         step, no_checkpoint, dry_run, force, keep_intermediate,
                         verbose, quiet, log_file, log_level,
                         skip_qc, skip_mapping, read1_pattern_fastp, read2_pattern_fastp):
    """
    基于GTX的全基因组重测序变异检测流程|Fastq to VCF (GTX) Pipeline

    示例|Examples: biopytools fastq2vcf-gtx -i /path/to/raw_fastq -g /path/to/genome.fa -o /path/to/output
    """

    # 延迟加载:仅在实际调用时导入|Lazy loading: import only when actually called
    fastq2vcf_gtx_main = _lazy_import_fastq2vcf_gtx_main()

    # 为原始main函数构建参数列表|Build argument list for original main function
    args = ['fastq2vcf_gtx.py']

    # 必需参数|Required parameters
    args.extend(['--input', input])
    args.extend(['-g', genome])
    args.extend(['--output-dir', output_dir])

    # 可选路径参数|Optional path parameters
    if clean_fastq_dir:
        args.extend(['--clean-fastq-dir', clean_fastq_dir])

    if mapping_dir:
        args.extend(['--mapping-dir', mapping_dir])

    if gvcf_dir:
        args.extend(['--gvcf-dir', gvcf_dir])

    if bam_dir:
        args.extend(['--bam-dir', bam_dir])

    if joint_dir:
        args.extend(['--joint-dir', joint_dir])

    if filter_dir:
        args.extend(['--filter-dir', filter_dir])

    # 线程参数|Thread parameters
    if threads != 12:
        args.extend(['--threads', str(threads)])

    # 过滤参数|Filtering parameters
    if min_depth != 5:
        args.extend(['--min-depth', str(min_depth)])

    if quality != 30:
        args.extend(['--quality', str(quality)])

    if indel_min_dp != 5:
        args.extend(['--indel-min-dp', str(indel_min_dp)])

    if indel_min_qual != 30:
        args.extend(['--indel-min-qual', str(indel_min_qual)])

    # 样本阈值参数|Sample thresholds
    if gtx_single_threshold != 200:
        args.extend(['--gtx-single-threshold', str(gtx_single_threshold)])

    if gtx_window_size != 20000000:
        args.extend(['--gtx-window-size', str(gtx_window_size)])

    # 工具路径参数|Tool paths
    if gtx_bin != '~/software/gtx/bin/gtx':
        args.extend(['--gtx-bin', gtx_bin])

    if gtx_cmd_gen_script != "${HOME}/software/scripts/51.生成GTX按染色体合并gVCF的脚本.sh":
        args.extend(['--gtx-cmd-gen-script', gtx_cmd_gen_script])

    # GTX WGS参数处理|GTX WGS parameters handling
    if not use_gtx_wgs:
        args.append('--no-gtx-wgs')

    if gtx_pcr_indel_model != 'CONSERVATIVE':
        args.extend(['--gtx-pcr-indel-model', gtx_pcr_indel_model])

    if gtx_min_confidence != 30:
        args.extend(['--gtx-min-confidence', str(gtx_min_confidence)])

    if gtx_min_base_qual != 20:
        args.extend(['--gtx-min-base-qual', str(gtx_min_base_qual)])

    if gtx_ploidy != 2:
        args.extend(['--gtx-ploidy', str(gtx_ploidy)])

    # 步骤控制|Step control
    if step:
        args.extend(['--step', step])

    # 执行控制|Execution control
    if no_checkpoint:
        args.append('--no-checkpoint')

    if dry_run:
        args.append('--dry-run')

    if force:
        args.append('--force')

    if keep_intermediate:
        args.append('--keep-intermediate')

    # 日志选项|Logging options
    if verbose > 0:
        for _ in range(verbose):
            args.append('--verbose')

    if quiet:
        args.append('--quiet')

    if log_file:
        args.extend(['--log-file', log_file])

    if log_level != 'INFO':
        args.extend(['--log-level', log_level])

    # 其他选项|Other options
    if skip_qc:
        args.append('--skip-qc')

    if skip_mapping:
        args.append('--skip-mapping')

    if read1_pattern_fastp != "_1.fq.gz":
        args.extend(['--read1-pattern-fastp', read1_pattern_fastp])

    if read2_pattern_fastp != "_2.fq.gz":
        args.extend(['--read2-pattern-fastp', read2_pattern_fastp])

    # 保存和恢复sys.argv|Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        # 调用原始main函数|Call original main function
        fastq2vcf_gtx_main()
    except SystemExit as e:
        # 处理正常程序退出|Handle normal program exit
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
