"""
HiFi+Hi-C基因组组装CLI包装器|HiFi+Hi-C Genome Assembly CLI Wrapper
"""

import click
import sys
import os
from pathlib import Path


def _lazy_import_assembler_main():
    """延迟加载assembler主函数|Lazy load assembler main function"""
    try:
        from ...hifi_hic.main import main as assembler_main
        return assembler_main
    except ImportError as e:
        def error_func():
            click.echo(f"导入错误|Import error: {e}", err=True)
            sys.exit(1)
        return error_func


@click.command(short_help='HiFi基因组组装(hifiasm)|HiFi Genome Assembly (hifiasm)',
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--hifi', '-i',
              required=True,
              type=click.Path(exists=True),
              help='HiFi数据文件路径|Path to HiFi data file')
@click.option('--hic-r1', '-1',
              type=click.Path(exists=True),
              help='Hi-C R1文件路径（可选）|Path to Hi-C R1 file (optional)')
@click.option('--hic-r2', '-2',
              type=click.Path(exists=True),
              help='Hi-C R2文件路径（可选）|Path to Hi-C R2 file (optional)')
@click.option('--prefix', '-p',
              default="genome_sample",
              type=str,
              show_default=True,
              help='样本前缀|Sample prefix')
@click.option('--threads', '-t',
              default=12,
              type=int,
              show_default=True,
              help='线程数|Number of threads')
@click.option('--genome-size', '-g',
              default="1.45g",
              type=str,
              show_default=True,
              help='预估基因组大小|Estimated genome size (e.g., 1.45g, 250m)')
@click.option('--n-hap',
              default=2,
              type=int,
              show_default=True,
              help='倍性|Ploidy (haploid count)')
@click.option('--purge-level', '-l',
              default=None,
              type=int,
              help='Purge level (0=no purging, 1=light, 2/3=aggressive) [default: 3 for unzip]')
@click.option('--hom-cov',
              default=None,
              type=int,
              help='Homozygous read coverage (--hom-cov) [default: auto]')
@click.option('--output', '-o',
              default="./assembly_output",
              type=click.Path(),
              show_default=True,
              help='输出目录|Output directory')
@click.option('--ngs',
              type=click.Path(exists=True, file_okay=False, dir_okay=True),
              help='NGS二代数据目录（可选）|NGS second-generation data directory (optional)')
@click.option('--ngs-pattern',
              default="_1.clean.fq.gz",
              type=str,
              show_default=True,
              help='NGS文件匹配模式|NGS file matching pattern (default: _1.clean.fq.gz)')
@click.option('--high-cov',
              default=95.0,
              type=float,
              show_default=True,
              help='高质量contig覆盖度阈值|High quality contig coverage threshold (default: 95.0)')
@click.option('--medium-cov-min',
              default=30.0,
              type=float,
              show_default=True,
              help='中等质量contig最小覆盖度|Medium quality contig minimum coverage (default: 30.0)')
@click.option('--no-purge-dups',
              is_flag=True,
              default=False,
              help='禁用Purge_Dups去冗余|Disable Purge_Dups deduplication (enabled by default)')
@click.option('--purge-dups-path',
              default='~/miniforge3/envs/purge_dups_v.1.2.6',
              type=str,
              show_default=True,
              help='Purge_Dups软件路径|Purge_Dups software path (default: ~/miniforge3/envs/purge_dups_v.1.2.6)')
@click.option('--purge-dups-threads',
              type=int,
              default=None,
              help='去冗余线程数|Deduplication threads (default: same as assembly threads)')
@click.option('--purge-dups-read-type',
              type=click.Choice(['pacbio', 'hifi', 'illumina']),
              default='hifi',
              show_default=True,
              help='去冗余reads类型|Deduplication reads type (default: hifi)')
@click.option('--no-resume',
              is_flag=True,
              default=False,
              help='禁用断点续传（强制重新运行所有步骤）|Disable resume mode (force rerun all steps)')
@click.option('--resume',
              is_flag=True,
              default=False,
              help='启用断点续传（默认已启用）|Enable resume mode (enabled by default)')
@click.option('--verbose', '-v',
              is_flag=True,
              help='详细输出模式|Verbose output mode')
def hifi_hic(hifi, hic_r1, hic_r2, prefix, threads, genome_size,
             n_hap, purge_level, hom_cov, output, ngs, ngs_pattern, high_cov, medium_cov_min,
             no_purge_dups, purge_dups_path, purge_dups_threads, purge_dups_read_type,
             no_resume, resume, verbose):
    """
    HiFi基因组组装流程|HiFi Genome Assembly Pipeline

    使用hifiasm进行基于HiFi数据的基因组组装（可选Hi-C数据）
    Genome assembly using hifiasm with HiFi data (optional Hi-C data)

    示例|Examples: biopytools hifi-hic -i hifi.fq -p sample1
    高级用法（禁用去冗余）|Advanced (disable dedup): biopytools hifi-hic -i hifi.fq -p sample1 --no-purge-dups
    """

    # 延迟加载|Lazy load
    assembler_main = _lazy_import_assembler_main()

    # 构建主函数参数列表|Build argument list for main function
    args = ['hifi_hic.py']

    # 必需参数|Required parameters
    args.extend(['--hifi', hifi])

    # 可选参数：Hi-C数据|Optional parameters: Hi-C data
    if hic_r1:
        args.extend(['--hic-r1', hic_r1])

    if hic_r2:
        args.extend(['--hic-r2', hic_r2])

    # 可选参数(仅在非默认时添加)|Optional parameters (add only when non-default)
    if prefix != "genome_sample":
        args.extend(['--prefix', prefix])

    if threads != 88:
        args.extend(['--threads', str(threads)])

    if genome_size != "1.45g":
        args.extend(['--genome-size', genome_size])

    if n_hap != 2:
        args.extend(['--n-hap', str(n_hap)])

    if purge_level is not None:
        args.extend(['--purge-level', str(purge_level)])

    if hom_cov is not None:
        args.extend(['--hom-cov', str(hom_cov)])

    if output != "./assembly_output":
        args.extend(['--output', output])

    # NGS polish参数|NGS polish parameters
    if ngs:
        args.extend(['--ngs', ngs])

    if ngs_pattern != "_1.clean.fq.gz":
        args.extend(['--ngs-pattern', ngs_pattern])

    if high_cov != 95.0:
        args.extend(['--high-cov', str(high_cov)])

    if medium_cov_min != 30.0:
        args.extend(['--medium-cov-min', str(medium_cov_min)])

    # Purge_Dups去冗余参数|Purge_Dups deduplication parameters
    if no_purge_dups:
        args.append('--no-purge-dups')

    if purge_dups_path != '~/miniforge3/envs/purge_dups_v.1.2.6':
        args.extend(['--purge-dups-path', purge_dups_path])

    if purge_dups_threads is not None:
        args.extend(['--purge-dups-threads', str(purge_dups_threads)])

    if purge_dups_read_type != 'hifi':
        args.extend(['--purge-dups-read-type', purge_dups_read_type])

    # 断点续传参数|Resume parameters
    if no_resume:
        args.append('--no-resume')

    if resume:
        args.append('--resume')

    if verbose:
        args.append('--verbose')

    # 保存并恢复sys.argv|Save and restore sys.argv
    original_argv = sys.argv

    try:
        # 调用主函数|Call main function
        sys.argv = args
        assembler_main()

    except SystemExit as e:
        # 处理系统退出|Handle system exit
        if e.code != 0:
            sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n用户中断|Assembly interrupted by user", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"运行时错误|Runtime error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
