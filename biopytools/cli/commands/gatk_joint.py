"""
GATK Joint Genotyping Command|GATK Joint Genotyping Command
"""

import click
import sys
import os


def _lazy_import_joint_main():
    """延迟加载Joint Genotyping主函数|Lazy load joint genotyping main function"""
    try:
        from ...gatk_joint.main import main as joint_main
        return joint_main
    except ImportError as e:
        click.echo(f"Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_path_exists(path):
    """验证路径存在性(仅在非帮助模式下)|Validate path existence (only in non-help mode)"""
    if not _is_help_request() and path and not os.path.exists(path):
        raise click.BadParameter(f"Path does not exist: {path}")
    return path


@click.command(
    short_help='GATK Joint Genotyping|GATK Joint Genotyping',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
# Required Parameters|必需参数
@click.option('--input', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_path_exists(value) if value else None,
              help='输入目录(包含VCF/GVCF文件)|Input directory (containing VCF/GVCF files)')
@click.option('--genome', '-g',
              required=True,
              callback=lambda ctx, param, value: _validate_path_exists(value) if value else None,
              help='参考基因组文件(.fasta/.fa)|Reference genome file (.fasta/.fa)')

# Output Control|输出控制
@click.option('--output-dir', '-o',
              default='./joint_genotyping_output',
              show_default=True,
              type=click.Path(),
              help='输出目录|Output directory')

# Computing Resources|计算资源
@click.option('--threads', '-t',
              type=int,
              default=12,
              show_default=True,
              help='线程数|Number of threads')
@click.option('--memory', '-m',
              default='100g',
              show_default=True,
              help='Java内存设置|Java memory setting')

# Interval Settings|区间设置
@click.option('--intervals', '-L',
              help='分析区间(染色体或区间文件)|Analysis intervals (chromosome or interval file)')

# SNP Filtering Parameters|SNP过滤参数
@click.option('--snp-qd',
              type=float,
              default=2.0,
              show_default=True,
              help='SNP QD阈值|SNP QD threshold')
@click.option('--snp-fs',
              type=float,
              default=60.0,
              show_default=True,
              help='SNP FS阈值|SNP FS threshold')
@click.option('--snp-mq',
              type=float,
              default=40.0,
              show_default=True,
              help='SNP MQ阈值|SNP MQ threshold')
@click.option('--snp-mqrs',
              type=float,
              default=-12.5,
              show_default=True,
              help='SNP MQRankSum阈值|SNP MQRankSum threshold')
@click.option('--snp-rprs',
              type=float,
              default=-8.0,
              show_default=True,
              help='SNP ReadPosRankSum阈值|SNP ReadPosRankSum threshold')
@click.option('--snp-sor',
              type=float,
              default=3.0,
              show_default=True,
              help='SNP SOR阈值|SNP SOR threshold')

# INDEL Filtering Parameters|INDEL过滤参数
@click.option('--indel-qd',
              type=float,
              default=2.0,
              show_default=True,
              help='INDEL QD阈值|INDEL QD threshold')
@click.option('--indel-fs',
              type=float,
              default=200.0,
              show_default=True,
              help='INDEL FS阈值|INDEL FS threshold')
@click.option('--indel-rprs',
              type=float,
              default=-20.0,
              show_default=True,
              help='INDEL ReadPosRankSum阈值|INDEL ReadPosRankSum threshold')
@click.option('--indel-sor',
              type=float,
              default=10.0,
              show_default=True,
              help='INDEL SOR阈值|INDEL SOR threshold')

# Tool Paths|工具路径
@click.option('--gatk-path',
              default='gatk',
              show_default=True,
              help='GATK软件路径|GATK software path')
@click.option('--bcftools-path',
              default='bcftools',
              show_default=True,
              help='BCFtools软件路径|BCFtools software path')
def gatk_joint(input, genome, output_dir, threads, memory, intervals,
                     snp_qd, snp_fs, snp_mq, snp_mqrs, snp_rprs, snp_sor,
                     indel_qd, indel_fs, indel_rprs, indel_sor,
                     gatk_path, bcftools_path):
    """
    GATK Joint Genotyping分析工具|GATK Joint Genotyping Tool

    示例|Examples: biopytools gatk-joint -i gvcf_folder/ -g ref.fasta -o results/
    """

    # 延迟加载|Lazy loading: import only when actually called
    joint_main = _lazy_import_joint_main()

    # 构建参数列表|Build argument list
    args = ['gatk_joint.py']
    args.extend(['-i', input])
    args.extend(['-g', genome])

    # Output parameters|输出参数
    if output_dir != './joint_genotyping_output':
        args.extend(['-o', output_dir])

    # Computing resources|计算资源
    if threads != 88:
        args.extend(['-t', str(threads)])
    if memory != '100g':
        args.extend(['-m', memory])

    # Interval settings|区间设置
    if intervals:
        args.extend(['-L', intervals])

    # SNP filtering parameters|SNP过滤参数
    if snp_qd != 2.0:
        args.extend(['--snp-qd', str(snp_qd)])
    if snp_fs != 60.0:
        args.extend(['--snp-fs', str(snp_fs)])
    if snp_mq != 40.0:
        args.extend(['--snp-mq', str(snp_mq)])
    if snp_mqrs != -12.5:
        args.extend(['--snp-mqrs', str(snp_mqrs)])
    if snp_rprs != -8.0:
        args.extend(['--snp-rprs', str(snp_rprs)])
    if snp_sor != 3.0:
        args.extend(['--snp-sor', str(snp_sor)])

    # INDEL filtering parameters|INDEL过滤参数
    if indel_qd != 2.0:
        args.extend(['--indel-qd', str(indel_qd)])
    if indel_fs != 200.0:
        args.extend(['--indel-fs', str(indel_fs)])
    if indel_rprs != -20.0:
        args.extend(['--indel-rprs', str(indel_rprs)])
    if indel_sor != 10.0:
        args.extend(['--indel-sor', str(indel_sor)])

    # Tool paths|工具路径
    if gatk_path != 'gatk':
        args.extend(['--gatk-path', gatk_path])
    if bcftools_path != 'bcftools':
        args.extend(['--bcftools-path', bcftools_path])

    # 保存和恢复sys.argv|Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        joint_main()
    except SystemExit as e:
        if e.code != 0:
            sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n用户中断|Interrupted by user", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"执行失败|Execution failed: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
