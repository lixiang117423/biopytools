"""
ALLHiC流水线CLI包装器 | ALLHiC Pipeline CLI Wrapper
"""

import click
import sys
import os
from pathlib import Path

# 懒加载导入 | Lazy import
def get_allhic_main():
    """懒加载获取ALLHiC主函数 | Lazy load ALLHiC main function"""
    try:
        from ...allhic.main import main as allhic_main
        return allhic_main
    except ImportError:
        try:
            from main import main as allhic_main
            return allhic_main
        except ImportError as e:
            def error_func():
                click.echo(f"❌ 导入错误 | Import Error: {e}", err=True)
                click.echo("请确保allhic模块已正确安装", err=True)
                sys.exit(1)
            return error_func


@click.command(short_help='ALLHiC Hi-C基因组支架构建流水线',
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('-r', '--reference', required=True,
              type=click.Path(exists=True),
              help='📄 参考基因组文件 | Reference genome file')
@click.option('-1', '--read1', required=True,
              type=click.Path(exists=True),
              help='📄 Hi-C读段1文件 | Hi-C read 1 file')
@click.option('-2', '--read2', required=True,
              type=click.Path(exists=True),
              help='📄 Hi-C读段2文件 | Hi-C read 2 file')
@click.option('-k', '--chr-num', required=True,
              type=int,
              help='🔢 染色体数量 | Number of chromosomes')
@click.option('-e', '--enzyme',
              default="GATC",
              type=str,
              help='✂️ ALLHiC酶切位点 | ALLHiC enzyme motif (default: GATC)')
@click.option('-t', '--threads',
              default=88,
              type=int,
              help='⚡ CPU线程数 | CPU threads (default: 88)')
@click.option('-w', '--workdir',
              default="./allhic_output",
              type=click.Path(),
              help='📂 工作目录 | Working directory (default: ./allhic_output)')
@click.option('--mapq-step1',
              default=1,
              type=int,
              help='📏 Step 1 MapQ阈值 | Step 1 MapQ threshold (default: 1)')
@click.option('--bin-size',
              default="500k",
              type=str,
              help='📊 二进制大小 | Bin size (default: 500k)')
@click.option('--min-bin-size',
              default="50k",
              type=str,
              help='📊 最小二进制大小 | Minimum bin size (default: 50k)')
@click.option('--skip-mapping',
              is_flag=True,
              help='⏭️ 跳过步骤1: 比对 | Skip Step 1: Mapping')
@click.option('--skip-allele',
              is_flag=True,
              help='⏭️ 跳过步骤1.5: 等位基因检测 | Skip Step 1.5: Allele Detection')
@click.option('--skip-prune',
              is_flag=True,
              help='⏭️ 跳过步骤2: 修剪 | Skip Step 2: Pruning')
@click.option('--skip-partition',
              is_flag=True,
              help='⏭️ 跳过步骤3: 分区 | Skip Step 3: Partition')
@click.option('--skip-extract',
              is_flag=True,
              help='⏭️ 跳过步骤3.5: 矩阵提取 | Skip Step 3.5: Extract Matrix')
@click.option('--skip-rescue',
              is_flag=True,
              help='⏭️ 跳过步骤4: 拯救 | Skip Step 4: Rescue')
@click.option('--skip-optimize',
              is_flag=True,
              help='⏭️ 跳过步骤5: 优化 | Skip Step 5: Optimization')
@click.option('--skip-build',
              is_flag=True,
              help='⏭️ 跳过步骤6: 构建 | Skip Step 6: Build FASTA')
@click.option('--skip-plot',
              is_flag=True,
              help='⏭️ 跳过步骤7: 绘图 | Skip Step 7: Plot Heatmap')
@click.option('--skip-asmkit',
              is_flag=True,
              help='⏭️ 跳过步骤8: JBAT生成 | Skip Step 8: JBAT Generation')
@click.option('--diagnose',
              is_flag=True,
              help='🔍 运行诊断模式 | Run diagnostic mode')
@click.option('--verbose', '-v',
              is_flag=True,
              help='📝 详细输出 | Verbose output')
def allhic(reference, read1, read2, chr_num, enzyme, threads, workdir,
                    mapq_step1, bin_size, min_bin_size, skip_mapping, skip_allele,
                    skip_prune, skip_partition, skip_extract, skip_rescue,
                    skip_optimize, skip_build, skip_plot, skip_asmkit,
                    diagnose, verbose):
    """
    🧬 ALLHiC Hi-C基因组支架构建流水线
    
    使用ALLHiC进行高效的Hi-C基因组支架构建，支持二倍体和多倍体基因组，
    包含完整的从比对到最终组装的流水线，以及asmkit JBAT文件生成。
    
    功能特点 | Features:
    - 基于ALLHiC的完整流水线
    - Hi-C数据支架构建
    - 等位基因检测和单倍型分离
    - 嵌合contig修剪
    - 接触矩阵分析和优化
    - asmkit JBAT文件生成
    - 支持跳过特定步骤
    
    流水线步骤 | Pipeline Steps:
    - Step 0: 数据准备和索引构建
    - Step 1: Hi-C读段比对
    - Step 1.5: 等位基因检测
    - Step 2: 嵌合contig修剪
    - Step 3: Contig分区
    - Step 3.5: 接触矩阵提取
    - Step 4: 未定位contig拯救
    - Step 5: Contig顺序优化
    - Step 6: 最终组装构建
    - Step 7: 接触图谱绘制
    - Step 8: asmkit JBAT生成
    
    输出目录结构 | Output Directory Structure:
    - 00_data/: 输入数据和索引文件
    - 01_mapping/: 比对结果
    - 02_allele_table/: 等位基因检测结果
    - 03_pruning/: 修剪结果
    - 04_partition/: 分区结果
    - 05_extract_matrix/: 接触矩阵
    - 06_rescue/: 拯救结果
    - 07_optimize/: 优化结果
    - 08_build/: 最终组装
    - 09_plot/: 绘图结果
    - 10_jbat_asmkit/: asmkit JBAT文件
    
    示例 | Examples:
    
    \b
    # 基本流水线运行
    biopytools allhic -r draft.asm.fasta \\
        -1 reads_R1.fastq.gz -2 reads_R2.fastq.gz -k 12
    
    \b
    # 自定义参数
    biopytools allhic --reference genome.fa \\
        --read1 hic_R1.fq.gz --read2 hic_R2.fq.gz \\
        --chr-num 24 --threads 88 --enzyme GATC
    
    \b
    # 跳过某些步骤
    biopytools allhic -r genome.fa -1 hic_R1.fq.gz -2 hic_R2.fq.gz \\
        -k 12 --skip-prune --skip-rescue
    
    \b
    # 诊断模式
    biopytools allhic -r genome.fa -1 hic_R1.fq.gz -2 hic_R2.fq.gz \\
        -k 12 --diagnose
    
    \b
    # 高性能配置
    biopytools allhic -r large_genome.fa \\
        -1 hic_R1.fastq.gz -2 hic_R2.fastq.gz \\
        -k 32 --threads 128 --workdir /fast_storage/allhic/
    
    输入文件要求 | Input File Requirements:
    
    参考基因组:
    - FASTA格式
    - Contig级别组装
    - 建议包含足够的N50
    
    Hi-C数据:
    - 双端FASTQ格式
    - 建议150bp或更长读长
    - 数据量建议 > 100X有效接触
    
    性能建议 | Performance Tips:
    - 大基因组推荐使用高线程数 (-t 88或更多)
    - 使用高速存储作为工作目录可显著提升性能
    - 确保有足够的RAM（建议 > 200GB）
    - 对于复杂基因组，考虑增加MapQ过滤阈值
    
    故障排除 | Troubleshooting:
    - 步骤失败：检查输入文件格式和质量
    - 内存不足：减少线程数或增加RAM
    - 磁盘空间不足：清理临时文件或使用更大存储
    - 软件依赖缺失：检查PATH中的软件安装
    
    输出解读 | Output Interpretation:
    
    最终组装:
    - groups.asm.fasta: 支架构建后的基因组
    - groups.agp: 组装描述文件
    
    质量评估:
    - 接触图谱显示染色体级别结构
    - N50和基因组大小统计
    - 装配完整性指标
    
    后续分析:
    - groups.assembly: JBAT格式（用于Juicebox）
    - groups.hic: Hi-C矩阵文件（用于可视化）
    """
    
    # 懒加载主函数
    allhic_main = get_allhic_main()
    
    # 构建参数列表传递给原始main函数
    args = ['allhic.py']
    
    # 必需参数
    args.extend(['-r', reference])
    args.extend(['-1', read1])
    args.extend(['-2', read2])
    args.extend(['-k', str(chr_num)])
    
    # 基本参数
    if enzyme != "GATC":
        args.extend(['-e', enzyme])
    
    if threads != 88:
        args.extend(['-t', str(threads)])
    
    if workdir != "./allhic_output":
        args.extend(['-w', workdir])
    
    if mapq_step1 != 1:
        args.extend(['--mapq-step1', str(mapq_step1)])
    
    if bin_size != "500k":
        args.extend(['--bin-size', bin_size])
    
    if min_bin_size != "50k":
        args.extend(['--min-bin-size', min_bin_size])
    
    # 跳过步骤选项
    if skip_mapping:
        args.append('--skip-mapping')
    
    if skip_allele:
        args.append('--skip-allele')
    
    if skip_prune:
        args.append('--skip-prune')
    
    if skip_partition:
        args.append('--skip-partition')
    
    if skip_extract:
        args.append('--skip-extract')
    
    if skip_rescue:
        args.append('--skip-rescue')
    
    if skip_optimize:
        args.append('--skip-optimize')
    
    if skip_build:
        args.append('--skip-build')
    
    if skip_plot:
        args.append('--skip-plot')
    
    if skip_asmkit:
        args.append('--skip-asmkit')
    
    # 其他选项
    if diagnose:
        args.append('--diagnose')
    
    if verbose:
        args.append('--verbose')
    
    # 保存并恢复sys.argv
    original_argv = sys.argv
    
    try:
        # 设置新的sys.argv
        sys.argv = args
        
        # 调用原始的main函数
        allhic_main()
        
    except SystemExit as e:
        # 处理程序正常退出
        if e.code != 0:
            sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n流水线被用户中断 | Pipeline interrupted by user", err=True)
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
    """ALLHiC流水线CLI工具组 | ALLHiC Pipeline CLI Tools"""
    if ctx.invoked_subcommand is None:
        click.echo(ctx.get_help())


@cli.command('allhic')
@click.option('-r', '--reference', required=True,
              type=click.Path(exists=True),
              help='📄 参考基因组文件 | Reference genome file')
@click.option('-1', '--read1', required=True,
              type=click.Path(exists=True),
              help='📄 Hi-C读段1文件 | Hi-C read 1 file')
@click.option('-2', '--read2', required=True,
              type=click.Path(exists=True),
              help='📄 Hi-C读段2文件 | Hi-C read 2 file')
@click.option('-k', '--chr-num', required=True,
              type=int,
              help='🔢 染色体数量 | Number of chromosomes')
@click.option('-e', '--enzyme',
              default="GATC",
              type=str,
              help='✂️ ALLHiC酶切位点 | ALLHiC enzyme motif (default: GATC)')
@click.option('-t', '--threads',
              default=88,
              type=int,
              help='⚡ CPU线程数 | CPU threads (default: 88)')
@click.option('-w', '--workdir',
              default="./allhic_output",
              type=click.Path(),
              help='📂 工作目录 | Working directory (default: ./allhic_output)')
@click.option('--diagnose',
              is_flag=True,
              help='🔍 运行诊断模式 | Run diagnostic mode')
def allhic_sub(reference, read1, read2, chr_num, enzyme, threads, workdir,
               diagnose):
    """ALLHiC Hi-C基因组支架构建流水线 | ALLHiC Hi-C Genome Scaffolding Pipeline"""
    allhic.callback(reference, read1, read2, chr_num, enzyme, threads,
                            workdir, 1, "500k", "50k", False, False, False,
                            False, False, False, False, False, False,
                            False, diagnose, False)


# 如果直接运行此脚本，提供独立接口
if __name__ == '__main__':
    allhic()