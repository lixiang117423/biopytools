"""
🧬 K-mer分析CLI包装器 | K-mer Analysis CLI Wrapper
优化版本：使用懒加载解决响应速度问题
"""

import click
import sys
import os


def _lazy_import_kmer_main():
    """懒加载kmer main函数 | Lazy load kmer main function"""
    try:
        from ...metagraph_kmer.main import main as kmer_main
        return kmer_main
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


def _validate_path(path_value):
    """验证路径（文件或文件夹）| Validate path (file or folder)"""
    if not _is_help_request() and path_value:
        if not os.path.exists(path_value):
            raise click.BadParameter(f"路径不存在 | Path does not exist: {path_value}")
    return path_value


@click.command(short_help="K-mer库构建与查询分析工具",
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('-r', '--reference',
              required=True,
              callback=lambda ctx, param, value: _validate_path(value) if value else None,
              help='🧬 参考序列文件或文件夹 | Reference sequence file or folder (FASTA/FASTQ)')
@click.option('-q', '--query',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='🔍 查询序列文件 | Query sequence file (FASTA/FASTQ)')
@click.option('-o', '--output',
              default='./kmer_output',
              type=click.Path(),
              help='📁 输出目录 | Output directory (default: ./kmer_output)')
@click.option('-k', '--kmer',
              type=int,
              default=31,
              help='🔢 k-mer长度 | K-mer length (default: 31)')
@click.option('-t', '--threads',
              type=int,
              default=88,
              help='🧵 线程数 | Number of threads (default: 88)')
@click.option('-m', '--memory',
              type=int,
              default=64,
              help='💾 内存限制(GB) | Memory limit in GB (default: 64)')
@click.option('--no-canonical',
              is_flag=True,
              help='❌ 禁用canonical k-mer处理 | Disable canonical k-mer processing')
@click.option('--min-count',
              type=int,
              default=1,
              help='📊 最小k-mer计数阈值 | Minimum k-mer count threshold (default: 1)')
@click.option('--kmc-path',
              default='kmc',
              help='🔧 KMC工具路径 | KMC tool path (default: kmc)')
@click.option('--kmc-tools-path',
              default='kmc_tools',
              help='🔧 KMC-tools路径 | KMC-tools path (default: kmc_tools)')
def metagraph_kmer(reference, query, output, kmer, threads, memory, no_canonical, 
         min_count, kmc_path, kmc_tools_path):
    """
    K-mer库构建与查询分析工具
    
    使用KMC进行高效的k-mer计数和分析，支持大规模基因组数据处理。
    自动处理canonical k-mer，提供坐标信息和丰度统计。
    
    功能特点 | Features:
    - 高效的k-mer库构建
    - 查询序列k-mer丰度统计
    - Canonical k-mer自动处理
    - 坐标信息保存（FASTA输入）
    - 多线程并行处理
    - 内存优化
    
    分析流程 | Analysis Pipeline:
    1. 构建参考k-mer库（使用KMC）
    2. 统计查询文件k-mer丰度
    3. 合并结果并生成最终输出
    
    输出文件 | Output Files:
    - final_output.txt: 最终的k-mer分析结果
    - coords.txt: k-mer坐标信息（FASTA参考）
    - counts.txt: k-mer计数统计
    - *.kmc_*: KMC中间文件
    - analysis.log: 分析日志
    
    示例 | Examples:
    
    \b
    # 基本用法
    biopytools metagraph-kmer -r reference.fasta -q query.fastq -o kmer_results
    
    \b
    # 自定义k-mer长度和线程数
    biopytools metagraph-kmer -r ref.fa -q query.fq -k 21 -t 64 -o results
    
    \b
    # 使用文件夹作为参考输入
    biopytools metagraph-kmer -r ref_folder/ -q reads.fastq.gz -o output
    
    \b
    # 高性能配置
    biopytools metagraph-kmer -r genome.fa -q sample.fq \\
        -k 31 -t 88 -m 64 -o kmer_analysis
    
    \b
    # 禁用canonical模式并设置计数阈值
    biopytools metagraph-kmer -r ref.fa -q query.fq \\
        --no-canonical --min-count 5
    
    \b
    # 自定义KMC工具路径
    biopytools metagraph-kmer -r reference.fasta -q query.fastq \\
        --kmc-path /opt/kmc/kmc \\
        --kmc-tools-path /opt/kmc/kmc_tools
    
    输入文件格式 | Input File Formats:
    - FASTA: 未压缩或gzip压缩 (.fa, .fasta, .fa.gz, .fasta.gz)
    - FASTQ: 未压缩或gzip压缩 (.fq, .fastq, .fq.gz, .fastq.gz)
    - 支持多个参考文件（使用文件夹输入）
    
    参考输入说明 | Reference Input:
    - 文件: 单个FASTA/FASTQ文件
    - 文件夹: 包含多个FASTA/FASTQ文件的目录
    - FASTA输入会保存k-mer坐标信息
    - FASTQ输入只进行k-mer计数
    
    K-mer参数建议 | K-mer Parameter Guidelines:
    
    k-mer长度选择:
    - k=21: 适用于短reads，较快的速度
    - k=31: 默认值，平衡特异性和灵敏度
    - k=51: 高特异性，适用于大基因组
    
    Canonical k-mer:
    - 默认启用，考虑正反链
    - 减少k-mer数量，节省内存
    - 适用于大多数应用场景
    
    计数阈值:
    - min-count=1: 保留所有k-mer（默认）
    - min-count=2+: 过滤低频k-mer，减少噪音
    
    性能优化 | Performance Optimization:
    
    线程数设置:
    - 建议使用可用CPU核心数
    - 默认88线程适用于高性能服务器
    - 根据实际硬件调整
    
    内存管理:
    - 默认64GB，适用于大多数基因组
    - 可根据数据规模和可用内存调整
    - KMC会自动优化内存使用
    
    存储建议:
    - 使用SSD提升I/O性能
    - 确保足够的临时空间
    - 中间文件可能占用较大空间
    
    应用场景 | Use Cases:
    - 基因组比较分析
    - 序列相似性检测
    - 污染序列筛查
    - 重复序列识别
    - 测序深度评估
    - 基因组组装辅助
    
    质量控制 | Quality Control:
    - 确保输入文件格式正确
    - 检查文件完整性
    - 验证k-mer长度合理性
    - 监控内存和磁盘使用
    - 查看日志文件排查问题
    
    常见问题 | Troubleshooting:
    - 内存不足: 减少k-mer长度或增加内存限制
    - 速度慢: 增加线程数或使用SSD存储
    - 磁盘空间不足: 清理临时文件或增加存储
    - KMC错误: 检查工具路径和版本
    
    依赖软件 | Dependencies:
    - KMC (K-mer Counter): 主要计数工具
    - KMC-tools: 辅助工具集
    
    注意事项 | Notes:
    - 确保KMC工具已正确安装并在PATH中
    - 大数据集可能需要较长处理时间
    - 建议预先进行数据质量控制
    - 输出文件可用于下游分析
    """
    
    # 🚀 懒加载：只有在实际调用时才导入模块 | Lazy loading: import only when actually called
    kmer_main = _lazy_import_kmer_main()
    
    # 构建参数列表传递给原始main函数 🔄 | Build argument list for original main function
    args = ['metagraph_kmer.py']
    
    # 必需参数 📋 | Required parameters
    args.extend(['-r', reference])
    args.extend(['-q', query])
    
    # 可选参数（只在非默认值时添加）⚙️ | Optional parameters (add only when non-default)
    if output != './kmer_output':
        args.extend(['-o', output])
    
    if kmer != 31:
        args.extend(['-k', str(kmer)])
    
    if threads != 88:
        args.extend(['-t', str(threads)])
    
    if memory != 64:
        args.extend(['-m', str(memory)])
    
    if min_count != 1:
        args.extend(['--min-count', str(min_count)])
    
    if kmc_path != 'kmc':
        args.extend(['--kmc-path', kmc_path])
    
    if kmc_tools_path != 'kmc_tools':
        args.extend(['--kmc-tools-path', kmc_tools_path])
    
    # 处理选项（布尔标志）🚩 | Processing options (boolean flags)
    if no_canonical:
        args.append('--no-canonical')
    
    # 保存并恢复sys.argv 💾 | Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args
    
    try:
        # 调用原始的main函数 🚀 | Call original main function
        kmer_main()
    except SystemExit as e:
        # 处理程序正常退出 ✅ | Handle normal program exit
        sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\nK-mer分析被用户中断 | K-mer analysis interrupted by user", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"❌ K-mer分析失败 | K-mer analysis failed: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv