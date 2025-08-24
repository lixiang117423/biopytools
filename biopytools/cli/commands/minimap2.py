"""
Minimap2分析命令 | Minimap2 Analysis Command
"""

import click
import sys
from ...minimap2.main import main as minimap2_main


@click.command(short_help='Minimap2全基因组比对与未比对区域提取',
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--target', '-t',
              required=True,
              type=click.Path(exists=True),
              help='目标基因组文件路径 | Target genome file path')
@click.option('--query', '-q',
              required=True,
              type=click.Path(exists=True),
              help='查询基因组文件路径 | Query genome file path')
@click.option('--output-dir', '-o',
              default='./minimap2_output',
              type=click.Path(),
              help='输出目录 | Output directory (default: ./minimap2_output)')
@click.option('--preset', '-x',
              default='asm5',
              type=click.Choice(['asm5', 'asm10', 'asm20', 'map-ont', 'map-pb']),
              help='Minimap2预设参数 | Minimap2 preset parameters (default: asm5)')
@click.option('--threads', '-p',
              default=8,
              type=int,
              help='线程数 | Number of threads (default: 8)')
@click.option('--min-match', '-m',
              default=1000,
              type=int,
              help='最小匹配长度阈值 | Minimum match length threshold (default: 1000)')
@click.option('--min-unmapped', '-u',
              default=1000,
              type=int,
              help='最小未比对区间长度阈值 | Minimum unmapped region length threshold (default: 1000)')
@click.option('--tp-type',
              default='P',
              type=click.Choice(['S', 'P', 'SP']),
              help='保留的tp类型 | tp type to keep: S(secondary), P(primary), SP(both) (default: P)')
@click.option('--minimap2-path', '-M',
              default='minimap2',
              type=str,
              help='minimap2可执行文件路径 | minimap2 executable path (default: minimap2)')
@click.option('--seqkit-path', '-S',
              default='seqkit',
              type=str,
              help='seqkit可执行文件路径 | seqkit executable path (default: seqkit)')
def minimap2(target, query, output_dir, preset, threads, min_match, min_unmapped,
             tp_type, minimap2_path, seqkit_path):
    """
    Minimap2全基因组比对和未比对区间提取工具
    
    使用Minimap2进行全基因组序列比对，识别和提取未比对的基因组区域。
    该工具特别适用于比较基因组学分析，帮助发现目标基因组中缺失的序列片段。
    
    功能特点 | Features:
    - 高效的全基因组序列比对
    - 智能识别未比对区域
    - 多种预设参数适配不同数据类型
    - 序列提取和格式化输出
    - 详细的统计报告生成
    
    输出文件 | Output Files:
    - alignment.paf: Minimap2比对结果文件
    - unmapped_regions.bed: 未比对区域的BED格式文件
    - unmapped_sequences.fasta: 未比对区域的序列文件
    - summary_report.txt: 分析统计报告
    
    示例 | Examples:
    
    \b
    # 基本基因组比对
    biopytools minimap2 -t target_genome.fasta -q query_genome.fasta -o results/
    
    \b
    # 高精度组装比对
    biopytools minimap2 --target reference.fa --query assembly.fa \\
        --preset asm5 --threads 16 -o comparison/
    
    \b
    # ONT长读长数据比对
    biopytools minimap2 -t genome.fasta -q ont_reads.fasta \\
        --preset map-ont --threads 32 -o ont_analysis/
    
    \b
    # PacBio长读长数据比对
    biopytools minimap2 -t reference.fa -q pacbio_reads.fa \\
        --preset map-pb -p 24 -o pacbio_results/
    
    \b
    # 自定义过滤参数
    biopytools minimap2 -t genome1.fa -q genome2.fa \\
        --min-match 5000 --min-unmapped 2000 \\
        --tp-type SP --threads 64 -o detailed_analysis/
    
    \b
    # 指定软件路径
    biopytools minimap2 -t target.fa -q query.fa \\
        --minimap2-path /opt/minimap2/minimap2 \\
        --seqkit-path /usr/local/bin/seqkit -o custom_analysis/
    
    \b
    # 低精度快速比对
    biopytools minimap2 -t large_genome.fa -q contigs.fa \\
        --preset asm20 --min-match 500 --threads 128 -o quick_scan/
    
    预设参数说明 | Preset Parameters:
    - asm5: 高精度基因组组装比对 (~0.1% 序列差异)
    - asm10: 中等精度基因组组装比对 (~1% 序列差异)  
    - asm20: 低精度基因组组装比对 (~5% 序列差异)
    - map-ont: Oxford Nanopore长读长数据比对
    - map-pb: PacBio长读长数据比对
    
    tp类型说明 | tp Type Description:
    - P (Primary): 主要比对，每个查询序列的最佳比对
    - S (Secondary): 次要比对，每个查询序列的次优比对
    - SP (Both): 同时保留主要和次要比对
    
    处理流程 | Processing Pipeline:
    1. 使用Minimap2进行全基因组比对生成PAF文件
    2. 解析PAF文件并根据参数筛选比对结果
    3. 识别查询基因组中的未比对区域
    4. 生成未比对区域的BED格式坐标文件
    5. 提取未比对区域的FASTA格式序列
    6. 生成详细的统计和总结报告
    
    应用场景 | Use Cases:
    - 基因组组装质量评估
    - 物种特异性序列识别
    - 基因组缺失片段检测
    - 比较基因组学分析
    - 长读长测序数据比对
    
    性能建议 | Performance Tips:
    - 大基因组建议使用更多线程 (--threads 32+)
    - 根据序列相似性选择合适的预设参数
    - 调整最小匹配长度以平衡敏感性和特异性
    - 使用SSD存储可提高I/O性能
    """
    
    # 构建参数列表传递给原始main函数
    args = ['minimap2.py']
    
    # 必需参数
    args.extend(['-t', target])
    args.extend(['-q', query])
    
    # 可选参数（只在非默认值时添加）
    if output_dir != './minimap2_output':
        args.extend(['-o', output_dir])
    
    if preset != 'asm5':
        args.extend(['-x', preset])
    
    if threads != 8:
        args.extend(['-p', str(threads)])
    
    if min_match != 1000:
        args.extend(['-m', str(min_match)])
    
    if min_unmapped != 1000:
        args.extend(['-u', str(min_unmapped)])
    
    if tp_type != 'P':
        args.extend(['--tp-type', tp_type])
    
    if minimap2_path != 'minimap2':
        args.extend(['-M', minimap2_path])
    
    if seqkit_path != 'seqkit':
        args.extend(['-S', seqkit_path])
    
    # 保存并恢复sys.argv
    original_argv = sys.argv
    sys.argv = args
    
    try:
        # 调用原始的main函数
        minimap2_main()
    except SystemExit as e:
        # 处理程序正常退出
        if e.code != 0:
            sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n分析流程被用户中断 | Analysis pipeline interrupted by user", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"运行错误 | Runtime error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv