"""
基因组组装CLI包装器 | Genome Assembly CLI Wrapper
"""

import click
import sys
import os
from pathlib import Path

# 懒加载导入 | Lazy import
def get_assembler_main():
    """懒加载获取assembler主函数 | Lazy load assembler main function"""
    try:
        from ...genome_assembler.main import main as assembler_main
        return assembler_main
    except ImportError:
        try:
            from genome_assembler.main import main as assembler_main
            return assembler_main
        except ImportError as e:
            def error_func():
                click.echo(f"❌ 导入错误 | Import Error: {e}", err=True)
                click.echo("请确保genome_assembler模块已正确安装", err=True)
                sys.exit(1)
            return error_func


@click.command(short_help='HiFi + Hi-C基因组组装工具',
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--hifi', '-i',
              required=True,
              type=click.Path(exists=True),
              help='🧬 HiFi数据文件路径 | Path to HiFi data file')
@click.option('--hic-r1', '-1',
              required=True,
              type=click.Path(exists=True),
              help='🧬 Hi-C R1文件路径 | Path to Hi-C R1 file')
@click.option('--hic-r2', '-2',
              required=True,
              type=click.Path(exists=True),
              help='🧬 Hi-C R2文件路径 | Path to Hi-C R2 file')
@click.option('--prefix', '-p',
              default="genome_sample",
              type=str,
              help='📝 样本前缀 | Sample prefix (default: genome_sample)')
@click.option('--threads', '-t',
              default=88,
              type=int,
              help='🧵 线程数 | Number of threads (default: 88)')
@click.option('--genome-size', '-g',
              default="1.45g",
              type=str,
              help='🧬 预估基因组大小 | Estimated genome size (default: 1.45g)')
@click.option('--n-hap',
              default=2,
              type=int,
              help='🔢 倍性 | Ploidy (default: 2)')
@click.option('--output', '-o',
              default="./assembly_output",
              type=click.Path(),
              help='📁 输出目录 | Output directory (default: ./assembly_output)')
@click.option('--verbose', '-v',
              is_flag=True,
              help='📝 详细输出 | Verbose output')
def hifi_hic(hifi, hic_r1, hic_r2, prefix, threads, genome_size, 
                    n_hap, output, verbose):
    """
    🧬 HiFi + Hi-C 基因组组装工具
    
    使用 hifiasm 进行高效的HiFi结合Hi-C基因组组装，支持二倍体和多倍体基因组，
    生成高质量的单倍型组装和主要组装结果。
    
    功能特点 | Features:
    - 基于hifiasm的高效组装
    - HiFi长读长数据支持
    - Hi-C scaffolding支持
    - 单倍型分离（二倍体）
    - GFA到FASTA格式转换
    - 详细的组装统计报告
    
    输出文件 | Output Files:
    - hap1.primary.fa: 单倍型1主要序列
    - hap2.primary.fa: 单倍型2主要序列
    - primary.fa: 主要组装序列
    - alternate.fa: 备选组装序列
    - assembly_statistics.txt: 组装统计报告
    
    目录结构 | Directory Structure:
    - 01.raw_output/: 原始GFA等文件
    - 02.fasta/: FASTA格式文件
    - 03.logs/: 日志文件
    - 04.statistics/: 统计信息
    
    示例 | Examples:
    
    \b
    # 基本组装
    biopytools hifi-hic -i hifi_data.fq \\
        -1 hic_R1.fastq.gz -2 hic_R2.fastq.gz \\
        -p sample1
    
    \b
    # 自定义参数组装
    biopytools hifi-hic --hifi long_reads.fq \\
        --hic-r1 hic_R1.fq.gz --hic-r2 hic_R2.fq.gz \\
        --prefix sample2 --threads 64 --genome-size 2.5g
    
    \b
    # 指定输出目录
    biopytools hifi-hic -i hifi.fq \\
        -1 hic_R1.fq.gz -2 hic_R2.fq.gz \\
        -p sample3 -o /data/assembly_results
    
    \b
    # 四倍体基因组组装
    biopytools hifi-hic --hifi tetraploid_hifi.fq \\
        --hic-r1 hic_R1.fq.gz --hic-r2 hic_R2.fq.gz \\
        --prefix tetra_genome --n-hap 4
    
    \b
    # 高性能组装（使用高速存储）
    biopytools hifi-hic -i hifi_data.fq \\
        -1 hic_R1.fq.gz -2 hic_R2.fq.gz \\
        -p large_genome -t 128 -g 10g \\
        -o /fast_storage/assembly/
    
    输入文件要求 | Input File Requirements:
    
    HiFi数据:
    - 支持FASTQ/FASTA格式
    - 建议N50 > 20kb
    - 数据量建议 > 30X覆盖度
    
    Hi-C数据:
    - 双端FASTQ格式
    - 建议150bp读长
    - 数据量建议 > 100X有效reads
    
    性能建议 | Performance Tips:
    - 大基因组推荐使用高线程数 (-t 88或更多)
    - 使用高速存储作为输出目录可提升性能
    - 确保有足够的RAM（建议 > 100GB）
    - 对于非常大的基因组，考虑增加--genome-size值
    
    输出解读 | Output Interpretation:
    
    单倍型组装:
    - hap1.primary.fa: 来自父亲单倍型
    - hap2.primary.fa: 来自母亲单倍型
    
    主要组装:
    - primary.fa: 最佳质量组装
    - alternate.fa: 备选组装版本
    
    统计信息:
    - N50: 组装连续性指标
    - 总长度: 组装基因组大小
    - 序列数: contig/scaffold数量
    
    故障排除 | Troubleshooting:
    - 组装失败：检查输入文件格式和路径
    - 内存不足：减少线程数或增加RAM
    - 磁盘空间不足：清理临时文件或使用更大存储
    """
    
    # 懒加载主函数
    assembler_main = get_assembler_main()
    
    # 构建参数列表传递给原始main函数
    args = ['hifi_hic.py']
    
    # 必需参数
    args.extend(['--hifi', hifi])
    args.extend(['--hic-r1', hic_r1])
    args.extend(['--hic-r2', hic_r2])
    
    # 可选参数（只在非默认值时添加）
    if prefix != "genome_sample":
        args.extend(['--prefix', prefix])
    
    if threads != 88:
        args.extend(['--threads', str(threads)])
    
    if genome_size != "1.45g":
        args.extend(['--genome-size', genome_size])
    
    if n_hap != 2:
        args.extend(['--n-hap', str(n_hap)])
    
    if output != "./assembly_output":
        args.extend(['--output', output])
    
    if verbose:
        args.append('--verbose')
    
    # 保存并恢复sys.argv
    original_argv = sys.argv
    
    try:
        # 设置新的sys.argv
        sys.argv = args
        
        # 调用原始的main函数
        assembler_main()
        
    except SystemExit as e:
        # 处理程序正常退出
        if e.code != 0:
            sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n组装被用户中断 | Assembly interrupted by user", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"运行错误 | Runtime error: {e}", err=True)
        sys.exit(1)
    finally:
        # 恢复原始sys.argv
        sys.argv = original_argv


@click.group(invoke_without_command=True)
@click.pass_context
def cli(ctx):
    """基因组组装CLI工具组 | Genome Assembly CLI Tools"""
    if ctx.invoked_subcommand is None:
        click.echo(ctx.get_help())


@cli.command('hifi-hic')
@click.option('--hifi', '-i',
              required=True,
              type=click.Path(exists=True),
              help='🧬 HiFi数据文件路径 | Path to HiFi data file')
@click.option('--hic-r1', '-1',
              required=True,
              type=click.Path(exists=True),
              help='🧬 Hi-C R1文件路径 | Path to Hi-C R1 file')
@click.option('--hic-r2', '-2',
              required=True,
              type=click.Path(exists=True),
              help='🧬 Hi-C R2文件路径 | Path to Hi-C R2 file')
@click.option('--prefix', '-p',
              default="genome_sample",
              type=str,
              help='📝 样本前缀 | Sample prefix (default: genome_sample)')
@click.option('--threads', '-t',
              default=88,
              type=int,
              help='🧵 线程数 | Number of threads (default: 88)')
@click.option('--genome-size', '-g',
              default="1.45g",
              type=str,
              help='🧬 预估基因组大小 | Estimated genome size (default: 1.45g)')
@click.option('--n-hap',
              default=2,
              type=int,
              help='🔢 倍性 | Ploidy (default: 2)')
@click.option('--output', '-o',
              default="./assembly_output",
              type=click.Path(),
              help='📁 输出目录 | Output directory (default: ./assembly_output)')
@click.option('--verbose', '-v',
              is_flag=True,
              help='📝 详细输出 | Verbose output')
def hifi_hic_sub(hifi, hic_r1, hic_r2, prefix, threads, genome_size, 
                        n_hap, output, verbose):
    """HiFi + Hi-C基因组组装 | HiFi + Hi-C Genome Assembly"""
    hifi_hic.callback(hifi, hic_r1, hic_r2, prefix, threads, 
                            genome_size, n_hap, output, verbose)


# 如果直接运行此脚本，提供独立接口
if __name__ == '__main__':
    hifi_hic()