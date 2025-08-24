"""
FASTQ质量控制命令 | FASTQ Quality Control Command
"""

import click
import sys
from ...fastp.main import main as fastp_main


@click.command(short_help = '批量运行fastp',
               context_settings=dict(help_option_names=['-h', '--help'],max_content_width=120))
@click.option('--input-dir', '-i',
              required=True,
              type=click.Path(exists=True, file_okay=False),
              help='输入原始FASTQ数据目录 | Input raw FASTQ data directory')
@click.option('--output-dir', '-o',
              required=True,
              type=click.Path(file_okay=False),
              help='输出清洁FASTQ数据目录 | Output clean FASTQ data directory')
@click.option('--fastp-path',
              default='fastp',
              help='fastp可执行文件路径 (默认: fastp) | fastp executable path (default: fastp)')
@click.option('--threads', '-t',
              type=int,
              default=12,
              help='线程数 (默认: 12) | Number of threads (default: 12)')
@click.option('--quality-threshold', '-q',
              type=int,
              default=30,
              help='质量阈值 (默认: 30) | Quality threshold (default: 30)')
@click.option('--min-length', '-l',
              type=int,
              default=50,
              help='最小长度 (默认: 50) | Minimum length (default: 50)')
@click.option('--unqualified-percent', '-u',
              type=int,
              default=40,
              help='不合格碱基百分比阈值 (默认: 40) | Unqualified base percentage threshold (default: 40)')
@click.option('--n-base-limit', '-n',
              type=int,
              default=10,
              help='N碱基数量限制 (默认: 10) | N base count limit (default: 10)')
@click.option('--read1-suffix',
              default='_1.fq.gz',
              help='Read1文件后缀 (默认: _1.fq.gz) | Read1 file suffix (default: _1.fq.gz)')
@click.option('--read2-suffix',
              default='_2.fq.gz',
              help='Read2文件后缀 (默认: _2.fq.gz) | Read2 file suffix (default: _2.fq.gz)')
def fastp(input_dir, output_dir, fastp_path, threads, quality_threshold, 
          min_length, unqualified_percent, n_base_limit, read1_suffix, read2_suffix):
    """
    FASTQ数据质量控制批处理工具
    
    使用fastp对FASTQ文件进行批量质量控制，支持双端测序数据的自动配对，
    生成质控后的清洁数据和详细的质控报告。
    
    示例 | Examples:
    
    \b
    # 基本质控处理
    biopytools fastp -i raw_data/ -o clean_data/
    
    \b
    # 自定义参数的质控处理
    biopytools fastp -i raw_data/ -o clean_data/ \\
        -t 24 -q 25 -l 100 -u 30 -n 5
    
    \b
    # 自定义文件后缀模式
    biopytools fastp -i raw_data/ -o clean_data/ \\
        --read1-suffix _R1.fastq.gz --read2-suffix _R2.fastq.gz
    
    \b
    # 指定fastp程序路径
    biopytools fastp -i raw_data/ -o clean_data/ \\
        --fastp-path /usr/local/bin/fastp -t 16
    """
    
    # 构建参数列表传递给原始main函数
    args = ['fastp.py']  # 模拟脚本名，避免传递biopytools子命令
    
    # 必需参数
    args.extend(['-i', input_dir])
    args.extend(['-o', output_dir])
    
    # 可选参数（只在非默认值时添加）
    if fastp_path != 'fastp':
        args.extend(['--fastp-path', fastp_path])
    
    if threads != 12:
        args.extend(['-t', str(threads)])
    
    if quality_threshold != 30:
        args.extend(['-q', str(quality_threshold)])
    
    if min_length != 50:
        args.extend(['-l', str(min_length)])
    
    if unqualified_percent != 40:
        args.extend(['-u', str(unqualified_percent)])
    
    if n_base_limit != 10:
        args.extend(['-n', str(n_base_limit)])
    
    if read1_suffix != '_1.fq.gz':
        args.extend(['--read1-suffix', read1_suffix])
    
    if read2_suffix != '_2.fq.gz':
        args.extend(['--read2-suffix', read2_suffix])
    
    # 保存并恢复sys.argv
    original_argv = sys.argv
    sys.argv = args
    
    try:
        # 调用原始的main函数
        fastp_main()
    except SystemExit as e:
        # 处理程序正常退出
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv