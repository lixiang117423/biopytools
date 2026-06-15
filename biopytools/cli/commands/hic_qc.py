"""
Hi-C数据质量控制评估命令|Hi-C Data Quality Control Assessment Command

"""

import click
import sys
import os


def _lazy_import_hicpro_qc_main():
    """延迟加载hicpro_qc主函数|Lazy load hicpro_qc main function"""
    try:
        from ...hicpro_qc.main import main as hicpro_qc_main
        return hicpro_qc_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _lazy_import_pairtools_qc_main():
    """延迟加载pairtools_qc主函数|Lazy load pairtools_qc main function"""
    try:
        from ...pairtools_qc.main import main as pairtools_qc_main
        return pairtools_qc_main
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


def _validate_dir_exists(dir_path):
    """验证目录存在(仅在非帮助模式)|Validate directory exists (only in non-help mode)"""
    if not _is_help_request() and dir_path and not os.path.exists(dir_path):
        raise click.BadParameter(f"目录不存在|Directory does not exist: {dir_path}")
    return dir_path


@click.command(
    short_help='Hi-C数据质量评估工具|Hi-C data quality assessment tool',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--input', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_dir_exists(value) if value and ctx.params.get('input_type') == 'hicpro' else (_validate_file_exists(value) if value and ctx.params.get('input_type') != 'hicpro' else value),
              help='输入文件或目录|Input file or directory (hicpro: HiC-Pro输出目录; pairs/bam: pairs或BAM文件)')
@click.option('--input-type',
              type=click.Choice(['hicpro', 'pairs', 'bam'],
                                case_sensitive=False),
              default='hicpro',
              show_default=True,
              help='输入类型|Input type (hicpro: HiC-Pro输出目录; pairs: .pairs.gz文件; bam: .bam文件)')
@click.option('--output-dir', '-o',
              default='./hic_qc_output',
              show_default=True,
              help='输出目录|Output directory')
@click.option('--sample-name', '-s',
              help='样本名称（仅hicpro模式，可选，默认自动检测）|Sample name (hicpro mode only, optional, auto-detect by default)')

# HiC-Pro模式阈值|HiC-Pro mode thresholds
@click.option('--min-mapping-rate',
              type=float,
              default=70.0,
              show_default=True,
              help='最低比对率阈值%%(仅hicpro模式)|Minimum mapping rate threshold%% (hicpro mode only)')
@click.option('--min-unique-rate',
              type=float,
              default=60.0,
              show_default=True,
              help='最低唯一比对率阈值%%(仅hicpro模式)|Minimum unique mapping rate threshold%% (hicpro mode only)')
@click.option('--min-valid-pairs-rate',
              type=float,
              default=50.0,
              show_default=True,
              help='最低valid pairs比例阈值%%(仅hicpro模式)|Minimum valid pairs rate threshold%% (hicpro mode only)')
@click.option('--max-dangling-ends-rate',
              type=float,
              default=15.0,
              show_default=True,
              help='最高dangling ends比例阈值%%(仅hicpro模式)|Maximum dangling ends rate threshold%% (hicpro mode only)')
@click.option('--max-self-ligation-rate',
              type=float,
              default=5.0,
              show_default=True,
              help='最高self-ligation比例阈值%%(仅hicpro模式)|Maximum self-ligation rate threshold%% (hicpro mode only)')
@click.option('--max-religation-rate',
              type=float,
              default=10.0,
              show_default=True,
              help='最高religation比例阈值%%(仅hicpro模式)|Maximum religation rate threshold%% (hicpro mode only)')

# Pairtools模式阈值|Pairtools mode thresholds
@click.option('--pairtools-path', '-p',
              default='~/miniforge3/envs/pairtools_v.1.1.3/bin/pairtools',
              show_default=True,
              help='Pairtools可执行文件路径（仅pairs/bam模式）|Pairtools executable path (pairs/bam mode only)')
@click.option('--chroms-path', '-c',
              help='Chromosome sizes文件路径（BAM输入时必需）|Chromosome sizes file path (required for BAM input)')
@click.option('--max-unmapped-rate',
              type=float,
              default=20.0,
              show_default=True,
              help='未比对reads比例阈值%%(仅pairs/bam模式)|Threshold for unmapped reads rate%% (pairs/bam mode only)')
@click.option('--max-single-sided-rate',
              type=float,
              default=10.0,
              show_default=True,
              help='单端比对比例阈值%%(仅pairs/bam模式)|Threshold for single-sided mapping rate%% (pairs/bam mode only)')
@click.option('--min-mapped-rate',
              type=float,
              default=80.0,
              show_default=True,
              help='双端比对率阈值%%(仅pairs/bam模式)|Threshold for paired mapping rate%% (pairs/bam mode only)')
@click.option('--max-dup-rate',
              type=float,
              default=30.0,
              show_default=True,
              help='PCR重复率阈值%%(仅pairs/bam模式)|Threshold for PCR duplication rate%% (pairs/bam mode only)')

# 通用阈值|Common thresholds
@click.option('--min-cis-trans-ratio',
              type=float,
              default=5.0,
              show_default=True,
              help='最低cis/trans比例阈值|Minimum cis/trans ratio threshold')
@click.option('--max-duplication-rate',
              type=float,
              default=30.0,
              show_default=True,
              help='最高PCR重复率阈值%%|Maximum PCR duplication rate threshold%%')
def hic_qc(input, input_type, output_dir, sample_name,
            min_mapping_rate, min_unique_rate, min_valid_pairs_rate,
            max_dangling_ends_rate, max_self_ligation_rate, max_religation_rate,
            pairtools_path, chroms_path, max_unmapped_rate,
            max_single_sided_rate, min_mapped_rate, max_dup_rate,
            min_cis_trans_ratio, max_duplication_rate):
    """
    Hi-C数据质量评估工具|Hi-C Data Quality Assessment Tool

    支持两种模式|Supports two modes:

    1. HiC-Pro模式（默认）|HiC-Pro mode (default):
       评估HiC-Pro输出目录中的统计文件
       Assess statistics from HiC-Pro output directory

       示例|Example: biopytools hic-qc -i /path/to/hicpro_output

    2. Pairtools模式|Pairtools mode:
       使用pairtools评估pairs或BAM文件
       Assess pairs or BAM files using pairtools

       示例|Example: biopytools hic-qc -i sample.pairs.gz --input-type pairs
    """

    # 根据输入类型选择不同的QC工具|Select QC tool based on input type
    if input_type.lower() == 'hicpro':
        # 使用HiCPro QC|Use HiCPro QC
        hicpro_qc_main_func = _lazy_import_hicpro_qc_main()

        # 构建参数列表|Build argument list
        args = ['hicpro_qc.py']

        # 必需参数|Required parameters
        args.extend(['-i', input])

        # 可选参数|Optional parameters
        if output_dir != './hic_qc_output':
            args.extend(['-o', output_dir])

        if sample_name:
            args.extend(['-s', sample_name])

        if min_mapping_rate != 70.0:
            args.extend(['--min-mapping-rate', str(min_mapping_rate)])

        if min_unique_rate != 60.0:
            args.extend(['--min-unique-rate', str(min_unique_rate)])

        if min_valid_pairs_rate != 50.0:
            args.extend(['--min-valid-pairs-rate', str(min_valid_pairs_rate)])

        if max_dangling_ends_rate != 15.0:
            args.extend(['--max-dangling-ends-rate', str(max_dangling_ends_rate)])

        if max_self_ligation_rate != 5.0:
            args.extend(['--max-self-ligation-rate', str(max_self_ligation_rate)])

        if max_religation_rate != 10.0:
            args.extend(['--max-religation-rate', str(max_religation_rate)])

        if min_cis_trans_ratio != 5.0:
            args.extend(['--min-cis-trans-ratio', str(min_cis_trans_ratio)])

        if max_duplication_rate != 30.0:
            args.extend(['--max-duplication-rate', str(max_duplication_rate)])

        # 执行主程序|Execute main program
        original_argv = sys.argv
        sys.argv = args

        try:
            hicpro_qc_main_func()
        except SystemExit as e:
            sys.exit(e.code)
        except Exception as e:
            click.echo(f"错误|Error: {e}", err=True)
            sys.exit(1)
        finally:
            sys.argv = original_argv

    else:  # pairs or bam mode
        # 使用Pairtools QC|Use Pairtools QC
        pairtools_qc_main_func = _lazy_import_pairtools_qc_main()

        # 构建参数列表|Build argument list
        args = ['pairtools_qc.py']

        # 必需参数|Required parameters
        args.extend(['-i', input])

        # 可选参数|Optional parameters
        if output_dir != './hic_qc_output':
            args.extend(['-o', output_dir])

        if pairtools_path != '~/miniforge3/envs/pairtools_v.1.1.3/bin/pairtools':
            args.extend(['-p', pairtools_path])

        if chroms_path:
            args.extend(['-c', chroms_path])

        if max_unmapped_rate != 20.0:
            args.extend(['--max-unmapped-rate', str(max_unmapped_rate)])

        if max_single_sided_rate != 10.0:
            args.extend(['--max-single-sided-rate', str(max_single_sided_rate)])

        if min_mapped_rate != 80.0:
            args.extend(['--min-mapped-rate', str(min_mapped_rate)])

        if max_dup_rate != 30.0:
            args.extend(['--max-dup-rate', str(max_dup_rate)])

        if min_cis_trans_ratio != 5.0:  # Note: pairtools默认是4.0，但这里统一用5.0
            args.extend(['--min-cis-trans-ratio', str(min_cis_trans_ratio)])

        # 执行主程序|Execute main program
        original_argv = sys.argv
        sys.argv = args

        try:
            pairtools_qc_main_func()
        except SystemExit as e:
            sys.exit(e.code)
        except Exception as e:
            click.echo(f"错误|Error: {e}", err=True)
            sys.exit(1)
        finally:
            sys.argv = original_argv
