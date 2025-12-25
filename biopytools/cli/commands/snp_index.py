"""
🧬 SNP index计算和分析命令 | SNP Index Calculation and Analysis Command
"""

import click
import sys
import os


def _lazy_import_snp_index_main():
    """懒加载snp_index main函数 | Lazy load snp_index main function"""
    try:
        from ...snp_index.main import main as snp_index_main
        return snp_index_main
    except ImportError as e:
        click.echo(f"❌ 导入错误 | Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否是帮助请求 | Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_file_exists(file_path):
    """验证文件是否存在（仅在非帮助模式下）| Validate file existence (only in non-help mode)"""
    if not _is_help_request() and file_path and not os.path.exists(file_path):
        raise click.BadParameter(f"文件不存在 | File does not exist: {file_path}")
    return file_path


def _validate_directory_exists(dir_path):
    """验证目录是否存在，不存在则创建 | Validate directory exists, create if not"""
    if not _is_help_request() and dir_path:
        os.makedirs(dir_path, exist_ok=True)
    return dir_path


@click.command(short_help="SNP index计算和分析工具",
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--input', '-i',
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='📄 输入VCF文件路径 | Input VCF file path')
@click.option('--output', '-o',
              callback=lambda ctx, param, value: _validate_directory_exists(value) if value else None,
              default='./snp_index_output',
              help='📁 输出目录路径 | Output directory path (default: ./snp_index_output)')
@click.option('--result-file', '-r',
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='📊 已有结果文件路径（用于分析模式）| Existing result file path (for analysis mode)')
@click.option('--prefix', '-p',
              default='snp_index',
              help='🏷️ 输出文件前缀 | Output file prefix (default: snp_index)')

# 过滤参数 | Filtering parameters
@click.option('--min-depth',
              type=int, default=10,
              help='🔍 最小测序深度 | Minimum sequencing depth (default: 10)')
@click.option('--min-quality',
              type=int, default=20,
              help='⭐ 最小质量值 | Minimum quality value (default: 20)')
@click.option('--min-mapping-quality',
              type=int, default=20,
              help='🎯 最小mapping质量 | Minimum mapping quality (default: 20)')
@click.option('--sample-names',
              nargs=2, type=str,
              help='👥 指定样本名称（样本1 样本2）| Specify sample names (sample1 sample2)')

# 分析参数 | Analysis parameters
@click.option('--extreme-threshold',
              type=float, default=0.8,
              help='📈 极端ΔSNP index阈值 | Extreme ΔSNP index threshold (default: 0.8)')
@click.option('--region-threshold',
              type=float, default=0.5,
              help='🗺️ 区域检测阈值 | Region detection threshold (default: 0.5)')
@click.option('--min-region-snps',
              type=int, default=5,
              help='🔢 区域最少SNP数量 | Minimum SNPs for region (default: 5)')
@click.option('--max-region-gap',
              type=int, default=10000,
              help='📏 区域最大gap(bp) | Maximum gap in region (default: 10000)')

# 滑动窗口参数 | Sliding window parameters
@click.option('--window-size',
              type=int, default=1000000,
              help='📊 滑动窗口大小(bp) | Sliding window size in bp (default: 1000000)')
@click.option('--step-size',
              type=int, default=100000,
              help='📈 滑动步长(bp) | Sliding step size in bp (default: 100000)')
@click.option('--min-window-snps',
              type=int, default=5,
              help='🔢 窗口最少SNP数 | Minimum SNPs per window (default: 5)')
@click.option('--confidence-level',
              type=float, default=0.95,
              help='📊 置信水平 | Confidence level (default: 0.95)')

# 可视化参数 | Visualization parameters
@click.option('--disable-sliding-window-plot',
              is_flag=True,
              help='📈 禁用滑动窗口折线图 | Disable sliding window line plot')
@click.option('--create-multi-chrom-plot',
              is_flag=True,
              help='🧬 创建多染色体分离图 | Create multi-chromosome separated plot')

# 模式选择 | Mode selection
@click.option('--calculate-only',
              is_flag=True,
              help='🧮 只计算SNP index，不进行分析 | Calculate SNP index only, no analysis')
@click.option('--analyze-only',
              is_flag=True,
              help='📊 只分析已有结果，不计算 | Analyze existing results only, no calculation')
@click.option('--skip-visualization',
              is_flag=True,
              help='📈 跳过可视化图表生成 | Skip visualization plot generation')

# 执行控制 | Execution control
@click.option('--threads', '-t',
              type=int, default=1,
              help='🔧 线程数 | Number of threads (default: 1)')
@click.option('--force', '-f',
              is_flag=True,
              help='💪 强制覆盖已存在文件 | Force overwrite existing files')
@click.option('--log-file',
              type=str,
              help='📝 日志文件路径 | Log file path')

# 日志控制 | Logging control
@click.option('--verbose', '-v',
              count=True,
              help='📢 详细输出模式 | Verbose mode (-v: INFO, -vv: DEBUG)')
@click.option('--quiet',
              is_flag=True,
              help='🤫 静默模式 | Quiet mode (ERROR only)')

def snp_index(input, output, result_file, prefix, min_depth, min_quality, min_mapping_quality,
              sample_names, extreme_threshold, region_threshold, min_region_snps, max_region_gap,
              window_size, step_size, min_window_snps, confidence_level, disable_sliding_window_plot,
              create_multi_chrom_plot, calculate_only, analyze_only, skip_visualization,
              threads, force, log_file, verbose, quiet):
    """
    SNP index计算和分析工具.

    从VCF文件计算SNP index和ΔSNP index，提供完整的结果分析和可视化功能。
    支持BSA（Bulked Segregant Analysis）分析。

    计算原理 | Calculation Principle:
    - SNP index = Alt_depth / (Ref_depth + Alt_depth)
    - ΔSNP index = SNP_index_pool1 - SNP_index_pool2

    功能特点 | Features:
    - 🧮 从VCF文件计算SNP index和ΔSNP index
    - 📊 详细的统计分析
    - 🗺️ 潜在目标区域检测
    - 📈 多种可视化图表
    - 🎯 **滑动窗口分析**：默认启用滑动窗口折线图分析
    - 📋 **滑动窗口结果**：自动保存滑动窗口计算结果
    - 🔧 灵活的过滤参数
    - 📋 支持多种输出格式

    示例 | Examples:

    \b
    # 🎯 完整分析流程（推荐使用）
    biopytools snp-index \\
        -i variation.filtered.snp.vcf.gz \\
        -o results/ \\
        -p my_analysis \\
        -v

    \b
    # 🔧 自定义过滤参数
    biopytools snp-index \\
        -i input.vcf.gz \\
        -o results/ \\
        --min-depth 20 \\
        --min-quality 30 \\
        --extreme-threshold 0.9

    \b
    # 📊 只分析已有结果文件
    biopytools snp-index \\
        --analyze-only \\
        -r existing_results.tsv \\
        -o analysis_output/

    \b
    # 🧮 只计算，不分析和可视化
    biopytools snp-index \\
        -i input.vcf.gz \\
        -o results/ \\
        --calculate-only \\
        --skip-visualization

    \b
    # 👥 指定特定样本进行分析
    biopytools snp-index \\
        -i input.vcf.gz \\
        -o results/ \\
        --sample-names pool_A pool_B

    输出文件 | Output Files:
    - 📄 {prefix}_results.tsv: SNP index计算结果
    - 📊 {prefix}_extreme_sites.tsv: 极端位点列表
    - 🗺️ {prefix}_potential_regions.tsv: 潜在目标区域
    - 📊 {prefix}_sliding_windows.tsv: **滑动窗口计算结果** | **Sliding window calculation results**
    - 📋 {prefix}_confidence_intervals.txt: 置信区间分析结果
    - 🗺️ {prefix}_candidate_regions.tsv: **候选区域识别结果** | **Candidate region identification results**
    - 📈 {prefix}_comprehensive.png: 综合分析图
    - 📊 {prefix}_manhattan.png: 曼哈顿图
    - 📉 {prefix}_delta_distribution.png: ΔSNP index分布
    - 📈 {prefix}_sliding_window.png: **滑动窗口折线图** | **Sliding window line plot**
    - 🧬 {prefix}_multi_chrom_sliding.png: 多染色体分离图（可选）
    """

    # 🚀 懒加载：只有在实际调用时才导入模块 | Lazy loading: import only when actually called
    snp_index_main = _lazy_import_snp_index_main()

    # 构建参数列表传递给原始main函数 🔄 | Build argument list for original main function
    args = ['snp_index.py']

    # 输入输出参数 | Input/output parameters
    if input:
        args.extend(['-i', input])

    args.extend(['-o', output])
    args.extend(['-p', prefix])

    if result_file:
        args.extend(['-r', result_file])

    # 过滤参数 | Filtering parameters
    if min_depth != 10:
        args.extend(['--min-depth', str(min_depth)])

    if min_quality != 20:
        args.extend(['--min-quality', str(min_quality)])

    if min_mapping_quality != 20:
        args.extend(['--min-mapping-quality', str(min_mapping_quality)])

    if sample_names:
        args.extend(['--sample-names'] + sample_names)

    # 分析参数 | Analysis parameters
    if extreme_threshold != 0.8:
        args.extend(['--extreme-threshold', str(extreme_threshold)])

    if region_threshold != 0.5:
        args.extend(['--region-threshold', str(region_threshold)])

    if min_region_snps != 5:
        args.extend(['--min-region-snps', str(min_region_snps)])

    if max_region_gap != 10000:
        args.extend(['--max-region-gap', str(max_region_gap)])

    # 滑动窗口参数 | Sliding window parameters
    if window_size != 1000000:
        args.extend(['--window-size', str(window_size)])

    if step_size != 100000:
        args.extend(['--step-size', str(step_size)])

    if min_window_snps != 5:
        args.extend(['--min-window-snps', str(min_window_snps)])

    if confidence_level != 0.95:
        args.extend(['--confidence-level', str(confidence_level)])

    # 可视化参数 | Visualization parameters
    if disable_sliding_window_plot:
        args.append('--disable-sliding-window-plot')

    if create_multi_chrom_plot:
        args.append('--create-multi-chrom-plot')

    # 模式选择 | Mode selection
    if calculate_only:
        args.append('--calculate-only')

    if analyze_only:
        args.append('--analyze-only')

    if skip_visualization:
        args.append('--skip-visualization')

    # 执行控制 | Execution control
    if threads != 1:
        args.extend(['-t', str(threads)])

    if force:
        args.append('--force')

    # 日志控制 | Logging control
    if log_file:
        args.extend(['--log-file', log_file])

    if verbose:
        args.extend(['-v'] * verbose)

    if quiet:
        args.append('--quiet')

    # 保存并恢复sys.argv 💾 | Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        # 调用原始的main函数 🚀 | Call original main function
        snp_index_main()
    except SystemExit as e:
        # 处理程序正常退出 ✅ | Handle normal program exit
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"❌ 错误 | Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv