"""
SNP index计算和分析命令|SNP Index Calculation and Analysis Command
"""

import click
import sys
import os


def _lazy_import_snp_index_main():
    """延迟加载snp_index主函数|Lazy load snp_index main function"""
    try:
        from ...snp_index.main import main as snp_index_main
        return snp_index_main
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


def _validate_directory_exists(dir_path):
    """验证目录存在，不存在则创建|Validate directory exists, create if not"""
    if not _is_help_request() and dir_path:
        os.makedirs(dir_path, exist_ok=True)
    return dir_path


@click.command(
    short_help='SNP index计算和分析工具|SNP index calculation and analysis tool',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--input', '-i',
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='输入VCF文件路径|Input VCF file path')
@click.option('--output', '-o',
              callback=lambda ctx, param, value: _validate_directory_exists(value) if value else None,
              default='./snp_index_output',
              show_default=True,
              help='输出目录路径|Output directory path')
@click.option('--result-file', '-r',
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='已有结果文件路径(用于分析模式)|Existing result file path (for analysis mode)')
@click.option('--prefix', '-p',
              default='snp_index',
              show_default=True,
              help='输出文件前缀|Output file prefix')
@click.option('--min-depth',
              type=int,
              default=10,
              show_default=True,
              help='最小测序深度|Minimum sequencing depth')
@click.option('--min-quality',
              type=int,
              default=20,
              show_default=True,
              help='最小质量值|Minimum quality value')
@click.option('--min-mapping-quality',
              type=int,
              default=20,
              show_default=True,
              help='最小mapping质量|Minimum mapping quality')
@click.option('--sample-names',
              nargs=2,
              type=str,
              help='指定样本名称(sample1 sample2)|Specify sample names (sample1 sample2)')
@click.option('--extreme-threshold',
              type=float,
              default=0.8,
              show_default=True,
              help='极端DeltaSNP index阈值|Extreme DeltaSNP index threshold')
@click.option('--region-threshold',
              type=float,
              default=0.5,
              show_default=True,
              help='区域检测阈值|Region detection threshold')
@click.option('--min-region-snps',
              type=int,
              default=5,
              show_default=True,
              help='区域最少SNP数量|Minimum SNPs for region')
@click.option('--max-region-gap',
              type=int,
              default=10000,
              show_default=True,
              help='区域最大gap(bp)|Maximum gap in region (bp)')
@click.option('--window-size',
              type=int,
              default=1000000,
              show_default=True,
              help='滑动窗口大小(bp)|Sliding window size in bp')
@click.option('--step-size',
              type=int,
              default=100000,
              show_default=True,
              help='滑动步长(bp)|Sliding step size in bp')
@click.option('--min-window-snps',
              type=int,
              default=5,
              show_default=True,
              help='窗口最少SNP数|Minimum SNPs per window')
@click.option('--confidence-level',
              type=float,
              default=0.95,
              show_default=True,
              help='置信水平|Confidence level')
@click.option('--disable-sliding-window-plot',
              is_flag=True,
              help='禁用滑动窗口折线图|Disable sliding window line plot')
@click.option('--create-multi-chrom-plot',
              is_flag=True,
              help='创建多染色体分离图|Create multi-chromosome separated plot')
@click.option('--calculate-only',
              is_flag=True,
              help='只计算SNP index，不分析|Calculate SNP index only, no analysis')
@click.option('--analyze-only',
              is_flag=True,
              help='只分析已有结果，不计算|Analyze existing results only, no calculation')
@click.option('--skip-visualization',
              is_flag=True,
              help='跳过可视化|Skip visualization')
@click.option('--threads', '-t',
              type=int,
              default=12,
              show_default=True,
              help='线程数|Number of threads')
@click.option('--force', '-f',
              is_flag=True,
              help='强制覆盖已存在文件|Force overwrite existing files')
@click.option('--log-file',
              type=str,
              help='日志文件路径|Log file path')
@click.option('--verbose', '-v',
              count=True,
              help='详细输出模式|Verbose mode (-v: INFO, -vv: DEBUG)')
@click.option('--quiet',
              is_flag=True,
              help='静默模式|Quiet mode (ERROR only)')
def snp_index(input, output, result_file, prefix, min_depth, min_quality, min_mapping_quality,
              sample_names, extreme_threshold, region_threshold, min_region_snps, max_region_gap,
              window_size, step_size, min_window_snps, confidence_level, disable_sliding_window_plot,
              create_multi_chrom_plot, calculate_only, analyze_only, skip_visualization,
              threads, force, log_file, verbose, quiet):
    """
    SNP index计算和分析工具|SNP Index Calculation and Analysis Tool

    从VCF文件计算SNP index和DeltaSNP index，用于BSA分析|Calculate SNP index and DeltaSNP index from VCF file for BSA analysis

    计算原理|Calculation Principle:
    - SNP index = Alt_depth / (Ref_depth + Alt_depth)
    - DeltaSNP index = SNP_index_pool1 - SNP_index_pool2

    示例|Examples: biopytools snp-index -i variation.filtered.snp.vcf.gz -o results/ -p my_analysis

    """

    # 延迟加载|Lazy load
    snp_index_main = _lazy_import_snp_index_main()

    # 构建主函数参数列表|Build argument list for main function
    args = ['snp_index.py']

    # 输入输出参数|Input/output parameters
    if input is not None:
        args.extend(['-i', input])

    if output != './snp_index_output':
        args.extend(['-o', output])

    if prefix != 'snp_index':
        args.extend(['-p', prefix])

    if result_file is not None:
        args.extend(['-r', result_file])

    # 过滤参数|Filtering parameters
    if min_depth != 10:
        args.extend(['--min-depth', str(min_depth)])

    if min_quality != 20:
        args.extend(['--min-quality', str(min_quality)])

    if min_mapping_quality != 20:
        args.extend(['--min-mapping-quality', str(min_mapping_quality)])

    if sample_names:
        args.extend(['--sample-names'] + list(sample_names))

    # 分析参数|Analysis parameters
    if extreme_threshold != 0.8:
        args.extend(['--extreme-threshold', str(extreme_threshold)])

    if region_threshold != 0.5:
        args.extend(['--region-threshold', str(region_threshold)])

    if min_region_snps != 5:
        args.extend(['--min-region-snps', str(min_region_snps)])

    if max_region_gap != 10000:
        args.extend(['--max-region-gap', str(max_region_gap)])

    # 滑动窗口参数|Sliding window parameters
    if window_size != 1000000:
        args.extend(['--window-size', str(window_size)])

    if step_size != 100000:
        args.extend(['--step-size', str(step_size)])

    if min_window_snps != 5:
        args.extend(['--min-window-snps', str(min_window_snps)])

    if confidence_level != 0.95:
        args.extend(['--confidence-level', str(confidence_level)])

    # 可视化参数|Visualization parameters
    if disable_sliding_window_plot:
        args.append('--disable-sliding-window-plot')

    if create_multi_chrom_plot:
        args.append('--create-multi-chrom-plot')

    # 模式选择|Mode selection
    if calculate_only:
        args.append('--calculate-only')

    if analyze_only:
        args.append('--analyze-only')

    if skip_visualization:
        args.append('--skip-visualization')

    # 执行控制|Execution control
    if threads != 1:
        args.extend(['-t', str(threads)])

    if force:
        args.append('--force')

    # 日志控制|Logging control
    if log_file is not None:
        args.extend(['--log-file', log_file])

    if verbose > 0:
        args.extend(['-v'] * verbose)

    if quiet:
        args.append('--quiet')

    # 保存并恢复sys.argv|Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        # 调用主函数|Call main function
        snp_index_main()
    except SystemExit as e:
        # 处理系统退出|Handle system exit
        sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n用户中断|User interrupted", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"运行时错误|Runtime error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
