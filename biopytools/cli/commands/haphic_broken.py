"""
HapHiC基因组scaffolding命令 | HapHiC Genome Scaffolding Command
"""

import click
import sys
import os


def _lazy_import_haphic_main():
    """懒加载HapHiC main函数 | Lazy load HapHiC main function"""
    try:
        from ...haphic.main import HapHiCProcessor
        return HapHiCProcessor
    except ImportError as e:
        click.echo(f"❌ 导入错误 | Import Error: {e}", err=True)
        sys.exit(1)


def _validate_file_exists(file_path):
    """验证文件是否存在 | Validate file existence"""
    if not file_path or not os.path.exists(file_path):
        raise click.BadParameter(f"文件不存在 | File does not exist: {file_path}")
    return file_path


@click.command(
    short_help='HapHiC基因组scaffolding工具',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.argument('asm_file',
              type=click.Path(exists=True, file_okay=True, dir_okay=False, resolve_path=True))
@click.argument('hic_file',
              type=click.Path(exists=True, file_okay=True, dir_okay=True, resolve_path=True))
@click.argument('nchrs',
              type=int)
@click.option('--hic-file-type',
              type=click.Choice(['bam', 'fastq']),
              default='bam',
              help='📊 Hi-C文件类型: bam=已比对文件, fastq=原始测序文件 (默认: bam) | Hi-C file type: bam=aligned file, fastq=raw sequencing files (default: bam)')
@click.option('--haphic-bin',
              default="haphic",
              help='🛠️ HapHiC可执行文件路径 (默认: haphic) | HapHiC executable path (default: haphic)')
@click.option('-o', '--output-dir',
              type=click.Path(file_okay=False, dir_okay=True, resolve_path=True),
              help='📂 输出目录路径 | Output directory path')
@click.option('--prefix',
              help='📝 输出文件前缀 | Output file prefix')
@click.option('--mapq-threshold', type=int, default=1,
              help='🎯 MAPQ阈值 (默认: 1) | MAPQ threshold (default: 1)')
@click.option('--edit-distance', type=int, default=3,
              help='🎯 编辑距离阈值 (默认: 3) | Edit distance threshold (default: 3)')
@click.option('--min-RE-sites', type=int, default=25,
              help='🎯 最小RE位点数 (默认: 25) | Minimum RE sites (default: 25)')
@click.option('--min-links', type=int, default=25,
              help='🎯 最小连接数 (默认: 25) | Minimum links (default: 25)')
@click.option('--min-link-density', type=float, default=0.0001,
              help='🎯 最小连接密度 (默认: 0.0001) | Minimum link density (default: 0.0001)')
@click.option('--min-inflation', type=float, default=1.0,
              help='🎯 最小膨胀参数 (默认: 1.0) | Min inflation (default: 1.0)')
@click.option('--max-inflation', type=float, default=3.0,
              help='🎯 最大膨胀参数 (默认: 3.0) | Max inflation (default: 3.0)')
@click.option('--inflation-step', type=float, default=0.2,
              help='🎯 膨胀参数步长 (默认: 0.2) | Inflation step (default: 0.2)')
@click.option('--Nx', type=int, default=80,
              help='🎯 Nx参数 (默认: 80) | Nx parameter (default: 80)')
@click.option('--min-group-len', type=int, default=0,
              help='🎯 最小分组长度 (默认: 0) | Min group length (default: 0)')
@click.option('--processes', type=int, default=8,
              help='🔧 并行进程数 (默认: 8) | Number of parallel processes (default: 8)')
@click.option('--threads', type=int, default=8,
              help='🔧 线程数 (默认: 8) | Number of threads (default: 8)')
@click.option('--correct-nrounds', type=int, default=2,
              help='🔧 组装校正轮数 (默认: 2) | Assembly correction rounds (default: 2)')
@click.option('--correct-min-coverage', type=float, default=10.0,
              help='🔧 校正最小覆盖度 (默认: 10.0) | Correction min coverage (default: 10.0)')
@click.option('--remove-allelic-links', type=int,
              help='🧬 移除等位基因连接数 | Remove allelic links count')
@click.option('--phasing-weight', type=float, default=1.0,
              help='🧬 分相权重 (默认: 1.0) | Phasing weight (default: 1.0)')
@click.option('--gfa-files',
              help='🧬 GFA文件路径(逗号分隔) | GFA files path (comma-separated)')
@click.option('--RE', default="GATC",
              help='🧬 限制性内切酶位点 (默认: GATC) | Restriction enzyme sites (default: GATC)')
@click.option('--bin-size', type=int, default=500,
              help='📊 接触图装箱大小 (默认: 500) | Contact map bin size (default: 500)')
@click.option('--min-len', type=float, default=1.0,
              help='📊 最小scaffold长度(默认: 1.0) | Min scaffold length (default: 1.0)')
@click.option('--quick-view', is_flag=True,
              help='🔍 快速查看模式 | Quick view mode')
@click.option('--generate-plots', is_flag=True,
              help='📈 生成可视化图表 | Generate visualization plots')
@click.option('--separate-plots', is_flag=True,
              help='📈 生成单独图表 | Generate separate plots')
@click.option('--no-fast-sorting', is_flag=True,
              help='⏭️ 禁用快速排序 | Disable fast sorting')
@click.option('--no-allhic-optimization', is_flag=True,
              help='⏭️ 禁用ALLHiC优化 | Disable ALLHiC optimization')
@click.option('--no-agp', is_flag=True,
              help='⏭️ 不输出AGP文件 | Don\'t output AGP file')
@click.option('--no-fasta', is_flag=True,
              help='⏭️ 不输出FASTA文件 | Don\'t output FASTA file')
@click.option('--no-juicebox', is_flag=True,
              help='⏭️ 不生成Juicebox脚本 | Don\'t generate Juicebox script')
@click.option('--continue-from',
              type=click.Choice(['cluster', 'reassign', 'sort', 'build']),
              help='🔄 从指定步骤继续 | Continue from specified step')
@click.option('--verbose', '-v', is_flag=True,
              help='📝 详细输出模式 | Verbose output mode')
@click.option('--log-file',
              type=click.Path(file_okay=True, dir_okay=False, resolve_path=True),
              help='📊 日志文件路径 | Log file path')
@click.option('--bwa-bin',
              default="bwa",
              help='🔗 BWA可执行文件路径 (默认: bwa) | BWA executable path (default: bwa)')
@click.option('--samtools-bin',
              default="samtools",
              help='🔧 Samtools可执行文件路径 (默认: samtools) | Samtools executable path (default: samtools)')
@click.option('--read1-pattern',
              default="_R1.fastq.gz",
              help='📄 Read1文件名模式 (默认: _R1.fastq.gz) | Read1 filename pattern (default: _R1.fastq.gz)')
@click.option('--read2-pattern',
              default="_R2.fastq.gz",
              help='📄 Read2文件名模式 (默认: _R2.fastq.gz) | Read2 filename pattern (default: _R2.fastq.gz)')
@click.option('--memory-limit',
              help='💾 内存限制 | Memory limit (e.g., 64G)')
@click.option('--dry-run', is_flag=True,
              help='🧪 测试模式，不执行实际命令 | Test mode, do not execute actual commands')
def haphic(asm_file, hic_file, nchrs, hic_file_type, haphic_bin, output_dir, prefix,
           mapq_threshold, edit_distance, min_re_sites, min_links, min_link_density,
           min_inflation, max_inflation, inflation_step, nx, min_group_len,
           processes, threads, correct_nrounds, correct_min_coverage,
           remove_allelic_links, phasing_weight, gfa_files, re, re_sites,
           bin_size, min_len, quick_view, generate_plots, separate_plots,
           no_fast_sorting, no_allhic_optimization, no_agp, no_fasta, no_juicebox,
           continue_from, verbose, log_file, bwa_bin, samtools_bin,
           read1_pattern, read2_pattern, memory_limit, dry_run, **kwargs):
    """
    HapHiC基因组scaffolding工具

    基于Hi-C数据的快速、参考基因组独立的等位基因感知scaffolding工具。
    支持单倍型分相基因组组装、二倍体和多倍体基因组组装，
    能够在1小时内完成大多数基因组的scaffolding工作。

    This is a fast, reference-independent, allele-aware scaffolding tool
    based on Hi-C data. It supports haplotype-phased, haplotype-collapsed
    diploid and allopolyploid genome assemblies, and can scaffold most
    genomes within 1 hour.
    """

    # 调试输出 - 显示额外的参数
    if kwargs:
        print(f"Debug - Extra kwargs: {kwargs}", file=sys.stderr)
        for key, value in kwargs.items():
            print(f"Debug - {key}: {value}", file=sys.stderr)

    参数 | Arguments:
        ASM_FILE          🧬 基因组组装文件路径 | Genome assembly file path (FASTA)
        HIC_FILE          📊 Hi-C文件路径 | Hi-C file path (BAM或FASTQ格式 | BAM or FASTQ format)
        NCHRS             🔢 染色体数量 | Number of chromosomes

    🧬 核心功能 | Core Features:
    - 染色体级别scaffolding，无需参考基因组
    - 高效校正嵌合contigs
    - 对嵌合、折叠和交换错误具有高容错性
    - 超快速且内存高效
    - 支持hifiasm分相信息
    - 生成可视化图表和Juicebox脚本

    🔄 四步流程 | Four-step Pipeline:
    1. 🔍 聚类分析 | Clustering: 使用马尔可夫聚类算法将contigs分组
    2. 🔄 重新分配 | Reassignment: 重新分配和重组contigs
    3. 📊 排序和定向 | Ordering & Orientation: 排序和定向contigs
    4. 🏗️ 构建scaffolds | Building: 构建最终的scaffold序列

    💡 使用场景 | Use Cases:
    - 单倍型分相基因组组装
    - 二倍体和多倍体基因组组装
    - 从contigs到染色体级别组装
    - 基因组错误校正
    - Hi-C接触图可视化

    📋 输入要求 | Input Requirements:
    - 基因组组装文件(FASTA格式)
    - Hi-C文件(BAM格式 - 按read name排序，或FASTQ格式 - 原始测序数据)
    - 染色体数量估计

  🔧 自动处理流程 | Automatic Processing:
    - FASTQ输入: 自动执行BWA比对 → 按read name排序 → HapHiC分析
    - BAM输入: 直接进行HapHiC分析(需按read name排序)

    📂 输出文件 | Output Files:
    - scaffolds.fa: 最终scaffold序列
    - scaffolds.agp: SALSA格式AGP文件
    - scaffolds.raw.agp: YaHS格式AGP文件
    - juicebox.sh: Juicebox可视化脚本
    - contact_map.pdf: Hi-C接触图可视化

    🎯 参数优化建议 | Parameter Optimization:
    - 小基因组: 降低--Nx值，减少--min_group_len
    - 大基因组: 增加线程数和内存限制
    - 高质量Hi-C: 使用更严格的过滤参数
    - 低质量Hi-C: 使用更宽松的参数

    ⚡ 性能优势 | Performance Advantages:
    - 比其他工具快10-100倍
    - 内存使用效率高
    - 自动参数优化
    - 支持断点续传

    示例 | Examples:

    \b
    # 基本scaffolding (使用BAM文件)
    biopytools haphic assembly.fa hic.bam 24

    \b
    # 使用FASTQ文件 (自动执行BWA比对)
    biopytools haphic assembly.fa hic_reads/ 24 \\
        --hic-file-type fastq

    \b
    # 单倍型分相组装 (使用FASTQ)
    biopytools haphic phased.fa hic_fastq/ 24 \\
        --hic-file-type fastq \\
        --gfa-files "hap1.p_ctg.gfa,hap2.p_ctg.gfa" \\
        --remove-allelic-links 2

    \b
    # 包含组装校正和可视化 (默认启用校正)
    biopytools haphic assembly.fa hic.bam 24 \\
        --generate-plots --bin-size 1000

    \b
    # 自定义酶切位点和BWA路径
    biopytools haphic assembly.fa hic.bam 24 \\
        --RE "AAGCTT" \\
        --bwa-bin /path/to/bwa

    \b
    # 从指定步骤继续
    biopytools haphic assembly.fa hic.bam 24 \\
        --continue-from sort

    \b
    # 大基因组高性能配置
    biopytools haphic assembly.fa hic.bam 24 \\
        --threads 32 --processes 16 \\
        --memory-limit 128G

    ⚠️ 重要提示 | Important Notes:
    - BAM文件必须按read name排序，不是坐标排序
    - 基因组文件可以是contigs或scaffolds
    - 快速查看模式适用于染色体数未知的情况
    - 组装校正现在通常不必要
    - 推荐使用YaHS-style AGP文件进行Juicebox可视化

    🔧 技术细节 | Technical Details:
    - 使用马尔可夫聚类算法(MCL)进行contig聚类
    - 集成3D-DNA和ALLHiC算法进行排序优化
    - 支持多种限制性内切酶
    - 自动处理分相信息
    - 内置错误检测和恢复机制
    """

    try:
        # 懒加载 | Lazy loading
        HapHiCProcessor = _lazy_import_haphic_main()

        # 创建处理器 | Create processor
        processor = HapHiCProcessor(
            asm_file=str(asm_file),
            hic_file=str(hic_file),
            nchrs=nchrs,
            hic_file_type=hic_file_type,
            haphic_bin=str(haphic_bin),
            output_dir=str(output_dir) if output_dir else None,
            prefix=prefix,
            mapq_threshold=mapq_threshold,
            edit_distance=edit_distance,
            min_re_sites=min_re_sites,
            min_links=min_links,
            min_link_density=min_link_density,
            min_inflation=min_inflation,
            max_inflation=max_inflation,
            inflation_step=inflation_step,
            nx=nx,
            min_group_len=min_group_len,
            processes=processes,
            threads=threads,
            fast_sorting=not no_fast_sorting,
            allhic_optimization=not no_allhic_optimization,
            correct_nrounds=correct_nrounds,
            correct_min_coverage=correct_min_coverage,
            remove_allelic_links=remove_allelic_links,
            phasing_weight=phasing_weight,
            gfa_files=gfa_files,
            re_sites=re,
            bin_size=bin_size,
            min_len=min_len,
            quick_view=quick_view,
            generate_plots=generate_plots,
            separate_plots=separate_plots,
            output_agp=not no_agp,
            output_fasta=not no_fasta,
            output_juicebox=not no_juicebox,
            continue_from_step=continue_from,
            verbose=verbose,
            log_file=str(log_file) if log_file else None,
            bwa_bin=bwa_bin,
            samtools_bin=samtools_bin,
            read1_pattern=read1_pattern,
            read2_pattern=read2_pattern,
            memory_limit=memory_limit
        )

        # 设置dry run模式
        if dry_run:
            processor.config.dry_run = True
            click.echo("🧪 测试模式: 不会执行实际的HapHiC命令")

        # 运行分析 | Run analysis
        if continue_from:
            success = processor.continue_from_step(continue_from)
        else:
            success = processor.run_pipeline()

        if success:
            click.echo("✅ HapHiC scaffolding完成!")
        else:
            click.echo("❌ HapHiC scaffolding失败!", err=True)
            sys.exit(1)

    except KeyboardInterrupt:
        click.echo("\n⚠️ 用户中断操作 | Operation interrupted by user", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"💥 发生错误 | Error occurred: {e}", err=True)
        sys.exit(1)