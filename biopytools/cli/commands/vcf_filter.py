"""
VCF筛选命令|VCF Filtering Command
"""

import click
import sys
import os


def _lazy_import_vcf_filter_main():
    """延迟加载vcf_filter主函数|Lazy load vcf_filter main function"""
    try:
        from ...vcf_filter.main import main as vcf_filter_main
        return vcf_filter_main
    except ImportError as e:
        click.echo(f"导入错误|Import error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_file_exists(file_path):
    """验证文件存在(仅在非帮助模式)|Validate file existence (only in non-help mode)"""
    if _is_help_request():
        return file_path
    if file_path and not os.path.exists(file_path):
        raise click.BadParameter(f"文件不存在|File does not exist: {file_path}")
    return file_path


@click.command(
    short_help='VCF筛选工具|VCF filtering tool',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--input', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='输入VCF文件路径|Input VCF file path')
@click.option('--output', '-o',
              type=click.Path(),
              help='输出VCF文件路径|Output VCF file path')
@click.option('--chr', '-c',
              required=True,
              type=str,
              help='染色体名称(支持逗号分隔多个)|Chromosome name(s) (comma-separated for multiple)')
@click.option('--start', '-s',
              type=int,
              help='起始位置|Start position')
@click.option('--end', '-e',
              type=int,
              help='结束位置|End position')
@click.option('--convert-format',
              is_flag=True,
              help='使用PLINK进行格式转换|Use PLINK for format conversion')
@click.option('--plink-path',
              default='plink',
              type=str,
              show_default=True,
              help='PLINK可执行文件路径|PLINK executable path')
@click.option('--allow-extra-chr',
              is_flag=True,
              default=True,
              show_default=True,
              help='允许额外染色体|Allow extra chromosomes')
@click.option('--maf',
              type=float,
              help='最小等位基因频率|Minimum allele frequency')
@click.option('--max-missing',
              type=float,
              help='最大缺失率|Maximum missing rate')
@click.option('--quality-threshold',
              type=float,
              help='质量阈值|Quality threshold')
@click.option('--min-depth',
              type=int,
              help='最小深度|Minimum depth')
@click.option('--max-depth',
              type=int,
              help='最大深度|Maximum depth')
@click.option('--keep-samples',
              type=str,
              help='保留样本名称(逗号分隔)|Sample names to keep (comma-separated)')
@click.option('--remove-samples',
              type=str,
              help='移除样本名称(逗号分隔)|Sample names to remove (comma-separated)')
@click.option('--keep-ids',
              type=str,
              help='保留变异位点ID(逗号分隔)|Variant IDs to keep (comma-separated)')
@click.option('--remove-ids',
              type=str,
              help='移除变异位点ID(逗号分隔)|Variant IDs to remove (comma-separated)')
@click.option('--biallelic-only',
              is_flag=True,
              help='只保留双等位基因位点|Keep only biallelic sites')
@click.option('--remove-indels',
              is_flag=True,
              help='移除插入缺失变异|Remove indel variants')
@click.option('--skip-validation',
              is_flag=True,
              default=True,
              show_default=True,
              help='跳过输入验证以提高速度|Skip input validation for speed')
@click.option('--force-validation',
              is_flag=True,
              help='强制执行输入验证|Force input validation')
@click.option('--verbose', '-v',
              is_flag=True,
              help='显示详细信息|Show verbose information')
def vcf_filter(input, output, chr, start, end, convert_format, plink_path, allow_extra_chr,
               maf, max_missing, quality_threshold, min_depth, max_depth,
               keep_samples, remove_samples, keep_ids, remove_ids,
               biallelic_only, remove_indels, skip_validation, force_validation, verbose):
    """
    VCF筛选工具|VCF Filtering Tool

    筛选VCF文件，支持区域提取、质量过滤和格式转换|Filter VCF files with support for region extraction, quality filtering, and format conversion

    示例|Examples: biopytools vcf-filter -i input.vcf -c chr1 -s 1000 -e 2000
    """

    # 延迟加载|Lazy load
    vcf_filter_main = _lazy_import_vcf_filter_main()

    # 构建主函数参数列表|Build argument list for main function
    args = ['vcf_filter.py']

    # 必需参数|Required parameters
    args.extend(['-i', input])
    args.extend(['-c', chr])

    # 可选参数(仅在非默认时添加)|Optional parameters (add only when non-default)
    if output:
        args.extend(['-o', output])

    if start is not None:
        args.extend(['-s', str(start)])

    if end is not None:
        args.extend(['-e', str(end)])

    if convert_format:
        args.append('--convert-format')

    if plink_path != 'plink':
        args.extend(['--plink-path', plink_path])

    if not allow_extra_chr:
        args.append('--allow-extra-chr')

    # 质量过滤参数|Quality filtering parameters
    if maf is not None:
        args.extend(['--maf', str(maf)])

    if max_missing is not None:
        args.extend(['--max-missing', str(max_missing)])

    if quality_threshold is not None:
        args.extend(['--quality-threshold', str(quality_threshold)])

    if min_depth is not None:
        args.extend(['--min-depth', str(min_depth)])

    if max_depth is not None:
        args.extend(['--max-depth', str(max_depth)])

    # 样本筛选参数|Sample filtering parameters
    if keep_samples:
        args.extend(['--keep-samples', keep_samples])

    if remove_samples:
        args.extend(['--remove-samples', remove_samples])

    # 变异位点筛选参数|Variant filtering parameters
    if keep_ids:
        args.extend(['--keep-ids', keep_ids])

    if remove_ids:
        args.extend(['--remove-ids', remove_ids])

    if biallelic_only:
        args.append('--biallelic-only')

    if remove_indels:
        args.append('--remove-indels')

    # 验证参数|Validation parameters
    if not skip_validation or force_validation:
        if force_validation:
            args.append('--force-validation')

    if verbose:
        args.append('--verbose')

    # 保存并恢复sys.argv|Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        # 调用主函数|Call main function
        vcf_filter_main()
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
