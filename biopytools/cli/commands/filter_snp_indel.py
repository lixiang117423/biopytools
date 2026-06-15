"""
VCF SNP/INDEL过滤命令|VCF SNP/INDEL Filtering Command
Optimized version: lazy loading for faster response
"""

import click
import sys
import os


def _lazy_import_vcf_filter_main():
    """延迟加载VCF过滤主函数|Lazy load VCF filter main function"""
    try:
        from ...filter_snp_indel.main import main as vcf_filter_main
        return vcf_filter_main
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
    short_help='VCF变异过滤工具|VCF variant filtering tool',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--input', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='输入VCF文件路径(支持.vcf/.vcf.gz)|Input VCF file path (supports .vcf/.vcf.gz)')
@click.option('--output-dir', '-o',
              default='./filtered_vcf',
              show_default=True,
              help='输出目录|Output directory')
@click.option('--threads', '-t',
              type=int,
              default=12,
              show_default=True,
              help='线程数|Number of threads')
@click.option('--snp-qual',
              type=float,
              default=30.0,
              show_default=True,
              help='SNP最小质量值|SNP minimum QUAL')
@click.option('--snp-dp',
              type=int,
              default=10,
              show_default=True,
              help='SNP最小测序深度|SNP minimum DP')
@click.option('--snp-mq',
              type=float,
              default=40.0,
              show_default=True,
              help='SNP最小比对质量|SNP minimum MQ')
@click.option('--snp-qd',
              type=float,
              default=2.0,
              show_default=True,
              help='SNP最小质量/深度比|SNP minimum QD')
@click.option('--snp-fs',
              type=float,
              default=60.0,
              show_default=True,
              help='SNP最大FisherStrand值|SNP maximum FS')
@click.option('--snp-sor',
              type=float,
              default=3.0,
              show_default=True,
              help='SNP最大StrandOddsRatio|SNP maximum SOR')
@click.option('--snp-mqrs',
              type=float,
              default=-12.5,
              show_default=True,
              help='SNP最小MappingQualityRankSum|SNP minimum MQRankSum')
@click.option('--snp-rprs',
              type=float,
              default=-8.0,
              show_default=True,
              help='SNP最小ReadPosRankSum|SNP minimum ReadPosRankSum')
@click.option('--snp-maf',
              type=float,
              default=0.05,
              show_default=True,
              help='SNP最小次等位基因频率|SNP minimum MAF')
@click.option('--snp-biallelic',
              is_flag=True,
              default=True,
              show_default=True,
              help='只保留双等位位点SNP|Keep only biallelic SNPs')
@click.option('--no-snp-biallelic',
              is_flag=True,
              flag_value=False,
              default=False,
              help='禁用双等位点过滤|Disable biallelic filtering')
@click.option('--indel-qual',
              type=float,
              default=30.0,
              show_default=True,
              help='INDEL最小质量值|INDEL minimum QUAL')
@click.option('--indel-dp',
              type=int,
              default=10,
              show_default=True,
              help='INDEL最小测序深度|INDEL minimum DP')
@click.option('--indel-mq',
              type=float,
              default=40.0,
              show_default=True,
              help='INDEL最小比对质量|INDEL minimum MQ')
@click.option('--indel-qd',
              type=float,
              default=2.0,
              show_default=True,
              help='INDEL最小质量/深度比|INDEL minimum QD')
@click.option('--indel-fs',
              type=float,
              default=200.0,
              show_default=True,
              help='INDEL最大FisherStrand值|INDEL maximum FS')
@click.option('--indel-sor',
              type=float,
              default=10.0,
              show_default=True,
              help='INDEL最大StrandOddsRatio|INDEL maximum SOR')
@click.option('--indel-rprs',
              type=float,
              default=-20.0,
              show_default=True,
              help='INDEL最小ReadPosRankSum|INDEL minimum ReadPosRankSum')
@click.option('--variant-type',
              type=click.Choice(['both', 'snp_only', 'indel_only']),
              default='both',
              show_default=True,
              help='输入VCF文件的变异类型|Variant type in input VCF')
@click.option('--bcftools-path',
              default='bcftools',
              show_default=True,
              help='BCFtools软件路径|BCFtools software path')
@click.option('--verbose', '-v',
              count=True,
              help='详细输出模式(-v: INFO, -vv: DEBUG)|Verbose mode (-v: INFO, -vv: DEBUG)')
@click.option('--quiet',
              is_flag=True,
              help='静默模式(仅输出ERROR)|Quiet mode (ERROR only)')
@click.option('--log-level',
              help='日志级别(DEBUG/INFO/WARNING/ERROR/CRITICAL)|Log level')
@click.option('--log-file',
              help='日志文件路径|Log file path')
@click.option('--force', '-f',
              is_flag=True,
              help='强制覆盖已存在文件|Force overwrite existing files')
@click.option('--dry-run',
              is_flag=True,
              help='模拟运行(不实际执行)|Dry run without execution')
@click.option('--repair-vcf',
              is_flag=True,
              help='自动修复损坏的VCF文件（列数不匹配等问题）|Auto-repair corrupted VCF files (column mismatch, etc.)')
def filter_snp_indel(input, output_dir, threads,
                     snp_qual, snp_dp, snp_mq, snp_qd, snp_fs, snp_sor, snp_mqrs, snp_rprs, snp_maf,
                     snp_biallelic, no_snp_biallelic,
                     indel_qual, indel_dp, indel_mq, indel_qd, indel_fs, indel_sor, indel_rprs,
                     variant_type,
                     bcftools_path,
                     verbose, quiet, log_level, log_file,
                     force, dry_run, repair_vcf):
    """
    VCF SNP/INDEL过滤工具|VCF SNP/INDEL Filtering Tool

    使用GATK最佳实践过滤VCF变异，自动分离SNP和INDEL并应用不同质控标准|Use GATK best practices for VCF variant filtering, automatically separate SNPs and INDELs and apply different quality control standards

    示例|Examples: biopytools filter-snp-indel -i variants.vcf -o filtered_output/
    """

    # 延迟加载|Lazy load
    vcf_filter_main = _lazy_import_vcf_filter_main()

    # 构建参数列表|Build argument list
    args = ['filter_snp_indel.py']
    args.extend(['-i', input])

    # 输出控制|Output control
    if output_dir != './filtered_vcf':
        args.extend(['-o', output_dir])
    if threads != 64:
        args.extend(['-t', str(threads)])

    # SNP过滤参数|SNP filtering parameters
    if snp_qual != 30.0:
        args.extend(['--snp-qual', str(snp_qual)])
    if snp_dp != 10:
        args.extend(['--snp-dp', str(snp_dp)])
    if snp_mq != 40.0:
        args.extend(['--snp-mq', str(snp_mq)])
    if snp_qd != 2.0:
        args.extend(['--snp-qd', str(snp_qd)])
    if snp_fs != 60.0:
        args.extend(['--snp-fs', str(snp_fs)])
    if snp_sor != 3.0:
        args.extend(['--snp-sor', str(snp_sor)])
    if snp_mqrs != -12.5:
        args.extend(['--snp-mqrs', str(snp_mqrs)])
    if snp_rprs != -8.0:
        args.extend(['--snp-rprs', str(snp_rprs)])
    if snp_maf != 0.05:
        args.extend(['--snp-maf', str(snp_maf)])

    # 处理双等位点过滤参数|Handle biallelic filtering parameters
    if no_snp_biallelic:
        args.append('--no-snp-biallelic')

    # INDEL过滤参数|INDEL filtering parameters
    if indel_qual != 30.0:
        args.extend(['--indel-qual', str(indel_qual)])
    if indel_dp != 10:
        args.extend(['--indel-dp', str(indel_dp)])
    if indel_mq != 40.0:
        args.extend(['--indel-mq', str(indel_mq)])
    if indel_qd != 2.0:
        args.extend(['--indel-qd', str(indel_qd)])
    if indel_fs != 200.0:
        args.extend(['--indel-fs', str(indel_fs)])
    if indel_sor != 10.0:
        args.extend(['--indel-sor', str(indel_sor)])
    if indel_rprs != -20.0:
        args.extend(['--indel-rprs', str(indel_rprs)])

    # 变异类型参数|Variant type parameter
    if variant_type != 'both':
        args.extend(['--variant-type', variant_type])

    # 工具路径|Tool path
    if bcftools_path != 'bcftools':
        args.extend(['--bcftools-path', bcftools_path])

    # 日志和执行控制参数|Logging and execution control parameters
    if verbose:
        args.extend(['-' + 'v' * verbose])

    if quiet:
        args.extend(['--quiet'])

    if log_level:
        args.extend(['--log-level', log_level])

    if log_file:
        args.extend(['--log-file', log_file])

    if force:
        args.extend(['--force'])

    if dry_run:
        args.extend(['--dry-run'])

    if repair_vcf:
        args.extend(['--repair-vcf'])

    # 执行主程序|Execute main program
    original_argv = sys.argv
    sys.argv = args

    try:
        vcf_filter_main()
    except SystemExit as e:
        if e.code != 0:
            sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n用户中断操作|Filtering interrupted by user", err=True)
        sys.exit(130)
    except Exception as e:
        click.echo(f"执行失败|Execution failed: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
