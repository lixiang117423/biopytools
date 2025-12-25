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
    short_help='HapHiC基因组scaffolding工具 - Pipeline模式',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)

# 主要输入参数 | Main input parameters
@click.option('--assembly', '-a',
              required=True,
              type=click.Path(exists=True, file_okay=True, dir_okay=False, resolve_path=True),
              help='🧬 基因组组装文件路径 (FASTA格式) | Genome assembly file path (FASTA format)')
@click.option('--bam', '-b',
              type=click.Path(exists=True, file_okay=True, dir_okay=False, resolve_path=True),
              help='📊 Hi-C BAM文件路径 (按read name排序) | Hi-C BAM file path (sorted by read name)')
@click.option('--hic1', '-1',
              type=click.Path(exists=True, file_okay=True, dir_okay=False, resolve_path=True),
              help='📄 Hi-C Read1文件路径 (FASTQ格式) | Hi-C Read1 file path (FASTQ format)')
@click.option('--hic2', '-2',
              type=click.Path(exists=True, file_okay=True, dir_okay=False, resolve_path=True),
              help='📄 Hi-C Read2文件路径 (FASTQ格式) | Hi-C Read2 file path (FASTQ format)')

# 必需参数 | Required parameters
@click.option('--chr-number', '-c', type=int, default=12,
              help='🔢 染色体数量 (默认: 12) | Number of chromosomes (default: 12)')

# 输出配置 | Output configuration
@click.option('--output-dir', '-o',
              type=click.Path(file_okay=False, dir_okay=True, resolve_path=True),
              help='📂 输出目录路径 | Output directory path')
@click.option('--prefix',
              help='📝 输出文件前缀 | Output file prefix')
@click.option('--force-rerun', is_flag=True,
              help='🔄 强制重新运行所有步骤(禁用断点续传) | Force rerun all steps (disable resume mode)')

# Hi-C数据处理参数 | Hi-C data processing parameters
@click.option('--mapq-threshold', type=int, default=1,
              help='🎯 MAPQ阈值 (默认: 1) | MAPQ threshold (default: 1)')
@click.option('--edit-distance', type=int, default=3,
              help='🎯 编辑距离阈值 (默认: 3) | Edit distance threshold (default: 3)')
@click.option('--min-RE-sites', type=int, default=25,
              help='🎯 最小RE位点数 (默认: 25) | Minimum RE sites (default: 25)')

# 聚类参数 | Clustering parameters
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

# 性能参数 | Performance parameters
@click.option('--processes', type=int, default=8,
              help='🔧 并行进程数 (默认: 8) | Number of parallel processes (default: 8)')
@click.option('--threads', type=int, default=8,
              help='🔧 线程数 (默认: 8) | Number of threads (default: 8)')
@click.option('--memory-limit',
              help='💾 内存限制 | Memory limit (e.g., 64G)')

# 组装校正参数 | Assembly correction parameters
@click.option('--correct-nrounds', type=int, default=2,
              help='🔧 组装校正轮数 (默认: 2) | Assembly correction rounds (default: 2)')
@click.option('--correct-min-coverage', type=float, default=10.0,
              help='🔧 校正最小覆盖度 (默认: 10.0) | Correction min coverage (default: 10.0)')

# ALLHiC优化参数 | ALLHiC optimization parameters
@click.option('--allhic-optimization/--no-allhic-optimization', default=False,
              help='🚀 启用ALLHiC优化 (默认: 禁用) | Enable ALLHiC optimization (default: disabled)')

# 单倍型分相参数 | Haplotype phasing parameters
@click.option('--remove-allelic-links', type=int,
              help='🧬 移除等位基因连接数 | Remove allelic links count')
@click.option('--phasing-weight', type=float, default=1.0,
              help='🧬 分相权重 (默认: 1.0) | Phasing weight (default: 1.0)')
@click.option('--gfa-files',
              help='🧬 GFA文件路径(逗号分隔) | GFA files path (comma-separated)')

# 可视化参数 | Visualization parameters
@click.option('--generate-plots', is_flag=True,
              help='📈 生成可视化图表 | Generate visualization plots')
@click.option('--bin-size', type=int, default=500,
              help='📊 接触图装箱大小 (默认: 500) | Contact map bin size (default: 500)')
@click.option('--min-len', type=float, default=1.0,
              help='📊 最小scaffold长度(默认: 1.0) | Min scaffold length (default: 1.0)')
@click.option('--separate-plots', is_flag=True,
              help='📈 生成单独图表 | Generate separate plots')

# 工具配置 | Tool configuration
@click.option('--haphic-bin',
              default="/share/org/YZWL/yzwl_lixg/miniforge3/envs/haphic/bin/haphic",
              help='🛠️ HapHiC可执行文件路径 (默认: /share/org/YZWL/yzwl_lixg/miniforge3/envs/haphic/bin/haphic) | HapHiC executable path (default: /share/org/YZWL/yzwl_lixg/miniforge3/envs/haphic/bin/haphic)')
@click.option('--bwa-bin',
              default="bwa",
              help='🔗 BWA可执行文件路径 (默认: bwa) | BWA executable path (default: bwa)')
@click.option('--samtools-bin',
              default="samtools",
              help='🔧 Samtools可执行文件路径 (默认: samtools) | Samtools executable path (default: samtools)')

# BWA比对配置 | BWA alignment configuration
@click.option('--samblaster-bin',
              default="samblaster",
              help='🔧 Samblaster可执行文件路径 (默认: samblaster) | Samblaster executable path (default: samblaster)')
@click.option('--haphic-filter-bam-bin',
              default="/share/org/YZWL/yzwl_lixg/miniforge3/envs/haphic/bin/filter_bam",
              help='🔍 HapHiC filter_bam工具路径 | HapHiC filter_bam tool path')
@click.option('--use-samblaster/--no-use-samblaster', default=True,
              help='🔧 使用samblaster去重 (默认: 启用) | Use samblaster deduplication (default: enabled)')
@click.option('--use-haphic-filter/--no-use-haphic-filter', default=True,
              help='🔍 使用HapHiC过滤 (默认: 启用) | Use HapHiC filtering (default: enabled)')

# Juicebox配置 | Juicebox configuration
@click.option('--generate-juicebox/--no-generate-juicebox', default=True,
              help='🥤 生成Juicebox兼容文件 (默认: 启用) | Generate Juicebox compatible files (default: enabled)')
@click.option('--matlock-bin',
              default="matlock",
              help='🔧 Matlock可执行文件路径 (默认: matlock) | Matlock executable path (default: matlock)')
@click.option('--three-d-dna-dir',
              default="/share/org/YZWL/yzwl_lixg/software/3d-dna",
              help='🧬 3D-DNA目录路径 (默认: /share/org/YZWL/yzwl_lixg/software/3d-dna) | 3D-DNA directory path (default: /share/org/YZWL/yzwl_lixg/software/3d-dna)')
@click.option('--agp2assembly-script',
              default="/share/org/YZWL/yzwl_lixg/software/3d-dna/utils/agp2assembly.py",
              help='📄 agp2assembly脚本路径 | agp2assembly script path')
@click.option('--asm-visualizer-script',
              default="/share/org/YZWL/yzwl_lixg/software/3d-dna/visualize/run-asm-visualizer.sh",
              help='📊 asm-visualizer脚本路径 | asm-visualizer script path')

# 高级选项 | Advanced options
@click.option('--RE', default="GATC",
              help='🧬 限制性内切酶位点 (默认: GATC) | Restriction enzyme sites (default: GATC)')
@click.option('--quick-view', is_flag=True,
              help='🔍 快速查看模式 | Quick view mode')

# 输出格式选项 | Output format options
@click.option('--no-agp', is_flag=True,
              help='⏭️ 不输出AGP文件 | Don\'t output AGP file')
@click.option('--no-fasta', is_flag=True,
              help='⏭️ 不输出FASTA文件 | Don\'t output FASTA file')
@click.option('--no-generate-plots', is_flag=True,
              help='⏭️ 不生成可视化图表 | Don\'t generate visualization plots')
@click.option('--no-juicebox', is_flag=True,
              help='⏭️ 不生成Juicebox脚本 | Don\'t generate Juicebox script')

# 日志配置 | Logging configuration
@click.option('--verbose', '-v', is_flag=True,
              help='📝 详细输出模式 | Verbose output mode')
@click.option('--log-file',
              type=click.Path(file_okay=True, dir_okay=False, resolve_path=True),
              help='📊 日志文件路径 | Log file path')

@click.option('--dry-run', is_flag=True,
              help='🧪 测试模式，不执行实际命令 | Test mode, do not execute actual commands')
def haphic(assembly, bam, hic1, hic2, chr_number, output_dir, prefix, force_rerun,
           mapq_threshold, edit_distance, min_re_sites,
           min_inflation, max_inflation, inflation_step, nx, min_group_len,
           processes, threads, memory_limit, correct_nrounds, correct_min_coverage,
           allhic_optimization,
           remove_allelic_links, phasing_weight, gfa_files, re,
           bin_size, min_len, generate_plots, separate_plots,
           haphic_bin, bwa_bin, samtools_bin,
           samblaster_bin, haphic_filter_bam_bin, use_samblaster, use_haphic_filter,
           generate_juicebox, matlock_bin, three_d_dna_dir, agp2assembly_script, asm_visualizer_script,
           quick_view,
           no_agp, no_fasta, no_generate_plots, no_juicebox,
           verbose, log_file, dry_run):
    """
    HapHiC基因组scaffolding工具 - Pipeline模式

    基于Hi-C数据的快速、参考基因组独立的等位基因感知scaffolding工具。
    现在使用HapHiC Pipeline模式，一步完成所有scaffolding步骤：
    1. BWA比对 (如果输入为FASTQ)
    2. HapHiC聚类 (cluster)
    3. 重新分配 (reassign)
    4. 排序和定向 (sort)
    5. 构建scaffolds (build)
    6. 生成可视化图表 (默认执行)
    7. 生成Juicebox文件 (可选)

    Pipeline模式优势 | Pipeline Mode Advantages:
        ✅ 更高效率: 减少文件I/O开销和进程启动时间
        ✅ 更好一致性: 使用HapHiC原生pipeline确保参数传递一致性
        ✅ 更简操作: 无需关心中间步骤，一键完成整个流程
        ✅ 更强鲁棒性: 原生pipeline有更好的错误处理和参数验证
        ✅ 断点续传: 自动检测已完成的步骤，支持中断恢复

    输入选项 | Input Options:
        必需提供以下输入之一:
        1. --bam: 已比对的BAM文件
        2. --hic1 + --hic2: Hi-C原始测序文件(FASTQ格式)

    自动处理流程 | Automatic Processing:
        BAM输入: 直接进行HapHiC Pipeline分析
        FASTQ输入: 自动执行BWA比对 → HapHiC Pipeline分析
          • BWA mem (-5SP) + samblaster + samtools view (-F 3340)
          • HapHiC filter_bam (MAPQ≥1, NM<3)
          • HapHiC Pipeline一步完成所有scaffolding步骤

    示例 | Examples:
        # 使用BAM文件
        biopytools haphic -a assembly.fa -b hic.bam -c 24

        # 使用FASTQ文件
        biopytools haphic -a assembly.fa -1 hic_R1.fastq.gz -2 hic_R2.fastq.gz -c 24

        # 高性能配置
        biopytools haphic -a assembly.fa -1 hic_R1.fastq.gz -2 hic_R2.fastq.gz -c 24 \\
            --threads 32 --processes 16 --correct-nrounds 2

        # 调整聚类参数
        biopytools haphic -a assembly.fa -b hic.bam -c 24 \\
            --min-inflation 1.0 --max-inflation 3.0 --inflation-step 0.2

        # 断点续传示例 (默认启用)
        biopytools haphic -a assembly.fa -b hic.bam -c 24
        # 如果中断后再次运行，会自动跳过已完成的步骤

        # 强制重新运行所有步骤
        biopytools haphic -a assembly.fa -b hic.bam -c 24 --force-rerun

    🚀 性能优化 (默认配置):
        - 默认跳过ALLHiC优化 (避免文件命名问题)
        - 使用FAST sorting算法 (高效且可靠)
        - 可通过--allhic-optimization启用ALLHiC优化

    🔄 断点续传功能 (默认启用):
        - 自动检测各步骤输出文件，判断是否已完成
        - HapHiC Pipeline步骤智能跳过
        - 可视化和Juicebox文件生成也支持断点续传
        - 使用--force-rerun禁用断点续传，强制重新运行

    📊 可视化图表生成 (默认执行):
        - 自动生成Hi-C接触图可视化
        - 存储在05.plots目录中
        - 可通过--no-generate-plots禁用

    🥤 Juicebox文件生成 (自动启用):
        - 自动生成Juicebox兼容的.hic和.assembly文件
        - 存储在06.juicebox目录中
        - 使用matlock、agp2assembly和asm-visualizer工具
        - 可通过--no-generate-juicebox禁用

    📂 输出文件 | Output Files:
        - 04.build/{prefix}.fa: 最终scaffold序列
        - 04.build/{prefix}.agp: SALSA格式AGP文件
        - 04.build/{prefix}.raw.agp: YaHS格式AGP文件
        - 05.plots/: 可视化图表目录 (PDF/PNG格式)
        - 06.juicebox/{prefix}.hic: Juicebox兼容的Hi-C文件
        - 06.juicebox/{prefix}.assembly: 3D-DNA assembly文件
        - {prefix}_haphic.log: 完整的运行日志
    """

    try:
        # 参数验证 | Parameter validation
        input_count = sum(1 for x in [bam, hic1, hic2] if x is not None)

        if input_count == 0:
            raise click.BadParameter("必须提供输入文件: 使用 --bam 或 --hic1/--hic2")
        elif bam and (hic1 or hic2):
            raise click.BadParameter("不能同时使用 --bam 和 --hic1/--hic2")
        elif (hic1 and not hic2) or (hic2 and not hic1):
            raise click.BadParameter("使用FASTQ输入时，必须同时提供 --hic1 和 --hic2")

        # 确定输入类型 | Determine input type
        if bam:
            hic_file = bam
            hic_file_type = "bam"
        else:
            hic_file = hic1  # 使用hic1作为主要文件
            hic_file_type = "fastq"

        # 懒加载 | Lazy loading
        HapHiCProcessor = _lazy_import_haphic_main()

        # 创建处理器 | Create processor
        processor = HapHiCProcessor(
            asm_file=str(assembly),
            hic_file=str(hic_file),
            nchrs=chr_number,
            hic_file_type=hic_file_type,
            haphic_bin=haphic_bin,
            output_dir=str(output_dir) if output_dir else None,
            prefix=prefix,
            force_rerun=force_rerun,
            mapq_threshold=mapq_threshold,
            edit_distance=edit_distance,
            min_re_sites=min_re_sites,
            min_inflation=min_inflation,
            max_inflation=max_inflation,
            inflation_step=inflation_step,
            nx=nx,
            min_group_len=min_group_len,
            processes=processes,
            threads=threads,
            correct_nrounds=correct_nrounds,
            correct_min_coverage=correct_min_coverage,
            allhic_optimization=allhic_optimization,
            remove_allelic_links=remove_allelic_links,
            phasing_weight=phasing_weight,
            gfa_files=gfa_files,
            re_sites=re,
            bin_size=bin_size,
            min_len=min_len,
            quick_view=quick_view,
            generate_plots=not no_generate_plots,
            separate_plots=separate_plots,
            output_agp=not no_agp,
            output_fasta=not no_fasta,
            output_juicebox=not no_juicebox,
            generate_juicebox=generate_juicebox,
            matlock_bin=matlock_bin,
            three_d_dna_dir=three_d_dna_dir,
            agp2assembly_script=agp2assembly_script,
            asm_visualizer_script=asm_visualizer_script,
                        verbose=verbose,
            log_file=str(log_file) if log_file else None,
            bwa_bin=bwa_bin,
            samtools_bin=samtools_bin,
            samblaster_bin=samblaster_bin,
            haphic_filter_bam_bin=haphic_filter_bam_bin,
            use_samblaster=use_samblaster,
            use_haphic_filter=use_haphic_filter,
            dry_run=dry_run,
            memory_limit=memory_limit
        )

        # 如果是FASTQ输入，需要传递hic2文件
        if hic_file_type == "fastq":
            processor.config.hic2_file = str(hic2)

        # 设置dry run模式
        if dry_run:
            processor.config.dry_run = True
            click.echo("🧪 测试模式: 不会执行实际的HapHiC命令")

        # 运行分析 | Run analysis
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