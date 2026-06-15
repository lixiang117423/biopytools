"""
BWA-GATK变异检测流程命令|BWA-GATK Pipeline Analysis Command
"""

import click
import sys
import os


def _lazy_import_bwa_gatk_main():
    """延迟加载BWA-GATK主函数|Lazy load BWA-GATK main function"""
    try:
        from ...bwa_gatk.main import main as bwa_gatk_main
        return bwa_gatk_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_path_exists(path, path_type="file"):
    """验证路径存在(仅在非帮助模式)|Validate path existence (only in non-help mode)"""
    if not _is_help_request() and not os.path.exists(path):
        raise click.BadParameter(
            f"{path_type}不存在|{path_type.capitalize()} does not exist: {path}"
        )
    return path


@click.command(
    short_help='BWA-GATK变异检测|BWA-GATK variant calling pipeline',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--input', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_path_exists(value, "input") if value else None,
              help='输入FASTQ文件目录（包含原始或清洁的FASTQ文件）|Input FASTQ directory (containing raw or clean FASTQ files)')
@click.option('--genome', '-g',
              required=True,
              callback=lambda ctx, param, value: _validate_path_exists(value, "genome") if value else None,
              help='参考基因组FASTA|Reference genome FASTA file')
@click.option('--output-dir', '-o',
              default='./bwa_gatk_output',
              show_default=True,
              type=click.Path(),
              help='输出目录|Output directory')
@click.option('--sample-name',
              help='样本名|Sample name')
@click.option('--threads', '-t',
              type=int,
              default=12,
              show_default=True,
              help='线程数|Number of threads')
@click.option('--memory',
              default='300G',
              show_default=True,
              help='最大内存使用|Maximum memory usage')
@click.option('--bwa-algorithm',
              type=click.Choice(['mem', 'aln']),
              default='mem',
              show_default=True,
              help='BWA比对算法|BWA alignment algorithm')
@click.option('--min-seed-length',
              type=int,
              default=19,
              show_default=True,
              help='最小种子长度|Minimum seed length')
@click.option('--band-width',
              type=int,
              default=100,
              show_default=True,
              help='带宽|Band width')
@click.option('--min-base-quality',
              type=int,
              default=20,
              show_default=True,
              help='最小碱基质量|Minimum base quality')
@click.option('--min-mapping-quality',
              type=int,
              default=20,
              show_default=True,
              help='最小比对质量|Minimum mapping quality')
@click.option('--stand-call-conf',
              type=float,
              default=30.0,
              show_default=True,
              help='调用置信度阈值|Stand call confidence threshold')
@click.option('--filter-expression',
              help='变异过滤表达式|Variant filter expression')
@click.option('--min-depth',
              type=int,
              default=10,
              show_default=True,
              help='最小测序深度|Minimum sequencing depth')
@click.option('--max-depth',
              type=int,
              help='最大测序深度|Maximum sequencing depth')
@click.option('--skip-bwa',
              is_flag=True,
              help='跳过BWA比对|Skip BWA alignment step')
@click.option('--skip-gatk',
              is_flag=True,
              help='跳过GATK变异检测|Skip GATK variant calling step')
@click.option('--only-align',
              is_flag=True,
              help='仅执行比对|Only perform alignment step')
@click.option('--only-call',
              is_flag=True,
              help='仅执行变异检测|Only perform variant calling step')
@click.option('--resume',
              is_flag=True,
              default=True,
              show_default=True,
              help='启用断点恢复|Enable resume from checkpoint')
@click.option('--force',
              is_flag=True,
              help='强制重新运行所有步骤|Force rerun all steps')
@click.option('--known-sites',
              type=click.Path(exists=True),
              help='已知变异位点VCF|Known variant sites VCF file')
@click.option('--bqsr',
              is_flag=True,
              help='启用碱基质量分数重校正|Enable Base Quality Score Recalibration')
@click.option('--remove-duplicates',
              is_flag=True,
              default=True,
              show_default=True,
              help='移除PCR重复|Remove PCR duplicates')
@click.option('--emit-ref-confidence',
              type=click.Choice(['NONE', 'BP_RESOLUTION', 'GVCF']),
              default='NONE',
              show_default=True,
              help='参考置信度模式|Reference confidence mode')
@click.option('--bwa-path',
              help='BWA工具路径|BWA tool path')
@click.option('--gatk-path',
              help='GATK工具路径|GATK tool path')
@click.option('--samtools-path',
              help='Samtools工具路径|Samtools tool path')
@click.option('--verbose', '-v',
              is_flag=True,
              help='详细输出模式|Verbose output mode')
@click.option('--quiet', '-q',
              is_flag=True,
              help='静默模式|Quiet mode')
@click.option('--keep-intermediate',
              is_flag=True,
              help='保留中间文件|Keep intermediate files')
@click.option('--skip-qc',
              is_flag=True,
              help='跳过质控步骤|Skip quality control step')
@click.option('--fastp-path',
              default='fastp',
              show_default=True,
              help='fastp可执行文件路径|fastp executable path')
@click.option('--qc-threads',
              type=int,
              default=12,
              show_default=True,
              help='质控线程数|QC threads')
@click.option('--qc-quality-threshold',
              type=int,
              default=20,
              show_default=True,
              help='质控质量阈值|QC quality threshold')
@click.option('--qc-min-length',
              type=int,
              default=50,
              show_default=True,
              help='质控最小长度|QC minimum length')
def bwa_gatk(input, genome, output_dir, sample_name, threads, memory,
             bwa_algorithm, min_seed_length, band_width,
             min_base_quality, min_mapping_quality, stand_call_conf,
             filter_expression, min_depth, max_depth,
             skip_bwa, skip_gatk, only_align, only_call, resume, force,
             known_sites, bqsr, remove_duplicates, emit_ref_confidence,
             bwa_path, gatk_path, samtools_path,
             verbose, quiet, keep_intermediate,
             skip_qc, fastp_path, qc_threads, qc_quality_threshold, qc_min_length):
    """
    BWA-GATK变异检测流程工具|BWA-GATK Variant Calling Pipeline Tool

    基于BWA比对和GATK变异检测的完整流程（包含质控）|Complete pipeline based on BWA alignment and GATK variant calling (with QC)

    示例|Examples: biopytools bwa-gatk -i fastq_dir/ -g ref.fasta -o output/
    """

    # 延迟加载|Lazy loading
    bwa_gatk_main = _lazy_import_bwa_gatk_main()

    # 构建参数列表|Build argument list
    args = ['bwa_gatk.py']
    args.extend(['-i', input])
    args.extend(['-g', genome])
    args.extend(['-o', output_dir])

    # 可选参数|Optional parameters
    if sample_name:
        args.extend(['--sample-name', sample_name])

    if threads != 88:
        args.extend(['-t', str(threads)])
    if memory != '300G':
        args.extend(['--memory', memory])

    # BWA参数|BWA parameters
    if bwa_algorithm != 'mem':
        args.extend(['--bwa-algorithm', bwa_algorithm])
    if min_seed_length != 19:
        args.extend(['--min-seed-length', str(min_seed_length)])
    if band_width != 100:
        args.extend(['--band-width', str(band_width)])

    # GATK参数|GATK parameters
    if min_base_quality != 20:
        args.extend(['--min-base-quality', str(min_base_quality)])
    if min_mapping_quality != 20:
        args.extend(['--min-mapping-quality', str(min_mapping_quality)])
    if stand_call_conf != 30.0:
        args.extend(['--stand-call-conf', str(stand_call_conf)])

    # 过滤参数|Filtering parameters
    if filter_expression:
        args.extend(['--filter-expression', filter_expression])
    if min_depth != 10:
        args.extend(['--min-depth', str(min_depth)])
    if max_depth:
        args.extend(['--max-depth', str(max_depth)])

    # 流程控制|Pipeline control
    if skip_bwa:
        args.append('--skip-bwa')
    if skip_gatk:
        args.append('--skip-gatk')
    if only_align:
        args.append('--only-align')
    if only_call:
        args.append('--only-call')
    if force:
        args.append('--force')
    if not resume:
        args.append('--no-resume')

    # 高级选项|Advanced options
    if known_sites:
        args.extend(['--known-sites', known_sites])
    if bqsr:
        args.append('--bqsr')
    if not remove_duplicates:
        args.append('--no-remove-duplicates')
    if emit_ref_confidence != 'NONE':
        args.extend(['--emit-ref-confidence', emit_ref_confidence])

    # 工具路径|Tool paths
    if bwa_path:
        args.extend(['--bwa-path', bwa_path])
    if gatk_path:
        args.extend(['--gatk-path', gatk_path])
    if samtools_path:
        args.extend(['--samtools-path', samtools_path])

    # 其他选项|Other options
    if verbose:
        args.append('-v')
    if quiet:
        args.append('-q')
    if keep_intermediate:
        args.append('--keep-intermediate')

    # 质控选项|QC options
    if skip_qc:
        args.append('--skip-qc')
    if fastp_path != 'fastp':
        args.extend(['--fastp-path', fastp_path])
    if qc_threads != 12:
        args.extend(['--qc-threads', str(qc_threads)])
    if qc_quality_threshold != 20:
        args.extend(['--qc-quality-threshold', str(qc_quality_threshold)])
    if qc_min_length != 50:
        args.extend(['--qc-min-length', str(qc_min_length)])

    # 保存和恢复sys.argv|Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        # 调用原始main函数|Call original main function
        bwa_gatk_main()
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
