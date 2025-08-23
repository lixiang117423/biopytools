"""
📊 BAM覆盖度分析命令 | BAM Depth Analysis Command
"""

import click
import sys
from ...depth_analyzer.main import main as depth_main


@click.command(short_help = "BAM/SAM文件覆盖度分析工具")
@click.option('--input', '-i',
              multiple=True,
              required=True,
              type=click.Path(exists=True),
              help='🎯 输入BAM/SAM文件路径或文件夹 | Input BAM/SAM file paths or directories')
@click.option('--output', '-o',
              default='./depth_results.txt',
              type=click.Path(),
              help='📄 输出文件路径或目录 (默认: ./depth_results.txt) | Output file path or directory (default: ./depth_results.txt)')
@click.option('--chromosome', '-c',
              default='all',
              help='🧬 目标染色体名称 (默认: all) | Target chromosome name (default: all)')
@click.option('--region', '-r',
              default='all',
              help='📍 染色体区间，格式: start:end (默认: all) | Chromosome region, format: start:end (default: all)')
@click.option('--threads', '-t',
              type=int,
              default=88,
              help='⚡ 线程数 (默认: 88) | Number of threads (default: 88)')
@click.option('--quality', '-q',
              type=int,
              default=0,
              help='⭐ 最小碱基质量阈值 (默认: 0) | Minimum base quality threshold (default: 0)')
@click.option('--mapping-quality', '-Q',
              type=int,
              default=0,
              help='🎯 最小比对质量阈值 (默认: 0) | Minimum mapping quality threshold (default: 0)')
@click.option('--samtools-path',
              default='samtools',
              help='🔧 samtools软件路径 (默认: samtools) | samtools software path (default: samtools)')
@click.option('--output-format',
              type=click.Choice(['txt']),
              default='txt',
              help='📄 输出格式 (默认: txt) | Output format (default: txt)')
@click.option('--compress',
              is_flag=True,
              help='🗜️ 压缩输出文件 | Compress output file')
@click.option('--enable-windows',
              is_flag=True,
              help='📊 启用滑窗分析 | Enable sliding window analysis')
@click.option('--window-size',
              type=int,
              default=1000,
              help='📏 窗口大小(bp) (默认: 1000) | Window size (bp) (default: 1000)')
@click.option('--window-step',
              type=int,
              default=0,
              help='👣 窗口步长(bp) (默认: 0) | Window step size (bp) (default: 0)')
def coverage(input, output, chromosome, region, threads, quality, mapping_quality,
             samtools_path, output_format, compress, enable_windows, window_size, window_step):
    """
    BAM/SAM文件覆盖度分析工具.
    
    分析BAM/SAM文件的测序覆盖度，支持指定染色体和区间分析，
    可选滑窗分析模式，自动处理索引文件。
    
    示例 | Examples:
    
    \b
    # 🎯 分析单个BAM文件
    biopytools bam-depth -i sample.bam -o depth_results.txt
    
    \b
    # 🧬 指定染色体和区间
    biopytools bam-depth -i sample.bam -o results.txt \\
        -c chr12 -r 136491092:138554123
    
    \b
    # 📁 批量分析多个文件
    biopytools bam-depth -i sample1.bam sample2.bam \\
        -o results.txt --threads 16
    
    \b
    # 📊 启用滑窗分析
    biopytools bam-depth -i sample.bam -o results.txt \\
        --enable-windows --window-size 1000 --window-step 500
    
    \b
    # 🎯 高质量比对过滤
    biopytools bam-depth -i data.bam -o results.txt \\
        -q 20 -Q 30 --compress
    """
    
    # 构建参数列表传递给原始main函数 🔄 | Build argument list for original main function
    args = ['biopytools', 'bam-depth']
    
    # 必需参数 📋 | Required parameters
    args.extend(['-i'] + list(input))
    
    # 可选参数 ⚙️ | Optional parameters
    if output != './depth_results.txt':
        args.extend(['-o', output])
    
    if chromosome != 'all':
        args.extend(['-c', chromosome])
    
    if region != 'all':
        args.extend(['-r', region])
    
    if threads != 88:
        args.extend(['-t', str(threads)])
    
    if quality != 0:
        args.extend(['-q', str(quality)])
    
    if mapping_quality != 0:
        args.extend(['-Q', str(mapping_quality)])
    
    if samtools_path != 'samtools':
        args.extend(['--samtools-path', samtools_path])
    
    if output_format != 'txt':
        args.extend(['--output-format', output_format])
    
    # 布尔选项 🚩 | Boolean options
    if compress:
        args.append('--compress')
    
    if enable_windows:
        args.append('--enable-windows')
    
    if window_size != 1000:
        args.extend(['--window-size', str(window_size)])
    
    if window_step != 0:
        args.extend(['--window-step', str(window_step)])
    
    # 保存并恢复sys.argv 💾 | Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args
    
    try:
        # 调用原始的main函数 🚀 | Call original main function
        depth_main()
    except SystemExit as e:
        # 处理程序正常退出 ✅ | Handle normal program exit
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"❌ 错误 | Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv