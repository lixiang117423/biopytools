"""
HapHiC基因组支架|HapHiC Genome Scaffolding Command
"""

import click
import sys
import os


def _lazy_import_haphic_main():
    """延迟加载HapHiC主函数|Lazy load HapHiC main function"""
    try:
        from ...haphic.main import HapHiCProcessor
        return HapHiCProcessor
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _validate_file_exists(file_path):
    """验证文件存在|Validate file existence"""
    if not file_path or not os.path.exists(file_path):
        raise click.BadParameter(f"文件不存在|File does not exist: {file_path}")
    return file_path


@click.command(
    short_help='HapHiC基因组支架流程|HapHiC Genome Scaffolding Pipeline',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)

# 主要输入参数|Main input parameters
@click.option('--input', '-i',
              required=True,
              type=click.Path(exists=True, file_okay=True, dir_okay=False, resolve_path=True),
              help='基因组组装文件(FASTA)|Genome assembly file path (FASTA format)')
@click.option('--bam', '-b',
              type=click.Path(exists=True, file_okay=True, dir_okay=False, resolve_path=True),
              help='Hi-C BAM文件(按read名排序)|Hi-C BAM file path (sorted by read name)')
@click.option('--hic1', '-1',
              type=click.Path(exists=True, file_okay=True, dir_okay=False, resolve_path=True),
              help='Hi-C Read1文件(FASTQ)|Hi-C Read1 file path (FASTQ format)')
@click.option('--hic2', '-2',
              type=click.Path(exists=True, file_okay=True, dir_okay=False, resolve_path=True),
              help='Hi-C Read2文件(FASTQ)|Hi-C Read2 file path (FASTQ format)')

# 必需参数|Required parameters
@click.option('--chr-number', '-c', type=int, default=12, show_default=True,
              help='染色体数量|Number of chromosomes')

# 输出配置|Output configuration
@click.option('--output-dir', '-o',
              type=click.Path(file_okay=False, dir_okay=True, resolve_path=True),
              default='.',
              show_default=True,
              help='输出目录路径|Output directory path')
@click.option('--prefix',
              help='输出文件前缀|Output file prefix')
@click.option('--force-rerun', is_flag=True,
              help='强制重新运行所有步骤|Force rerun all steps (disable resume mode)')

# Hi-C数据处理参数|Hi-C data processing parameters
@click.option('--mapq-threshold', type=int, default=1, show_default=True,
              help='MAPQ阈值|MAPQ threshold')
@click.option('--edit-distance', type=int, default=3, show_default=True,
              help='编辑距离阈值|Edit distance threshold')
@click.option('--re-site-cutoff', type=int, default=5, show_default=True,
              help='Step1 RE位点过滤阈值|Step1 RE site filtering threshold')
@click.option('--min-RE-sites', type=int, default=25, show_default=True,
              help='Step2重分配最小RE位点数|Step2 reassignment min RE sites')
@click.option('--aln-format', type=click.Choice(['auto', 'bam', 'pairs']), default='auto', show_default=True,
              help='比对文件格式|Alignment file format')

# 聚类参数|Clustering parameters
@click.option('--min-inflation', type=float, default=1.1, show_default=True,
              help='最小膨胀值|Min inflation')
@click.option('--max-inflation', type=float, default=3.0, show_default=True,
              help='最大膨胀值|Max inflation')
@click.option('--inflation-step', type=float, default=0.1, show_default=True,
              help='膨胀值步长|Inflation step')
@click.option('--Nx', type=int, default=80, show_default=True,
              help='Nx参数|Nx parameter')
@click.option('--min-group-len', type=float, default=5.0, show_default=True,
              help='最小组长度(Mbp)|Min group length (Mbp)')
@click.option('--flank', type=int, default=500, show_default=True,
              help='邻接矩阵侧翼区域(kbp)|Adjacency matrix flank region (kbp)')
@click.option('--bin-size-kbp', type=int, default=-1, show_default=True,
              help='聚类分箱大小(kbp),-1=自动|Clustering bin size (kbp), -1=auto')

# 性能参数|Performance parameters
@click.option('--processes', type=int, default=8, show_default=True,
              help='并行进程数|Number of parallel processes')
@click.option('--threads', '-t', type=int, default=12, show_default=True,
              help='线程数|Number of threads')
@click.option('--memory-limit', default='100G', show_default=True,
              help='内存限制|Memory limit (e.g., 64G, 300G)')

# 组装修正参数|Assembly correction parameters
@click.option('--correct-nrounds', type=int, default=2, show_default=True,
              help='组装修正轮数(0=禁用)|Assembly correction rounds (0=disabled)')
@click.option('--correct-min-coverage', type=float, default=10.0, show_default=True,
              help='修正最小覆盖度|Correction min coverage')
@click.option('--median-cov-ratio', type=float, default=0.2, show_default=True,
              help='覆盖率截断乘数|Coverage cutoff multiplier')
@click.option('--region-len-ratio', type=float, default=0.1, show_default=True,
              help='高覆盖区域长度比|High-coverage region length ratio')
@click.option('--min-region-cutoff', type=int, default=5000, show_default=True,
              help='高覆盖区域最小长度(bp)|Min high-coverage region length (bp)')

# ALLHiC优化参数|ALLHiC optimization parameters
@click.option('--skip-fast-sort', is_flag=True,
              help='跳过快速排序|Skip fast sorting')
@click.option('--skip-allhic', is_flag=True,
              help='跳过ALLHiC优化|Skip ALLHiC optimization')
@click.option('--skip-ga', is_flag=True,
              help='跳过ALLHiC遗传算法|Skip ALLHiC genetic algorithm')
@click.option('--sort-by-input', is_flag=True,
              help='按输入顺序排序|Sort output by input order')
@click.option('--no-additional-rescue', is_flag=True,
              help='跳过额外救援轮|Skip additional rescue round')
@click.option('--remove-concentrated-links', is_flag=True,
              help='移除高密度集中链接|Remove concentrated links')
@click.option('--normalize-by-nlinks', is_flag=True,
              help='按链接数归一化|Normalize by number of links')
@click.option('--dense-matrix', is_flag=True,
              help='使用稠密矩阵|Use dense matrix')

# 单倍型相位参数|Haplotype phasing parameters
@click.option('--remove-allelic-links', type=int,
              help='移除等位基因连锁数|Remove allelic links count')
@click.option('--phasing-weight', type=float, default=1.0, show_default=True,
              help='相位权重|Phasing weight')
@click.option('--gfa-files',
              help='GFA文件路径(逗号分隔)|GFA files path (comma-separated)')

# 可视化参数|Visualization parameters
@click.option('--generate-plots', is_flag=True,
              help='生成可视化图表|Generate visualization plots')
@click.option('--bin-size', type=int, default=500, show_default=True,
              help='接触图bin大小|Contact map bin size')
@click.option('--min-len', type=float, default=1.0, show_default=True,
              help='最小scaffold长度|Min scaffold length')
@click.option('--separate-plots', is_flag=True,
              help='生成单独的图表|Generate separate plots')

# 工具配置|Tool configuration
@click.option('--haphic-bin',
              default="~/miniforge3/envs/haphic/bin/haphic",
              show_default=True,
              help='HapHiC可执行文件路径|HapHiC executable path')
@click.option('--bwa-bin',
              default="~/miniforge3/envs/Population_genetics/bin/bwa",
              show_default=True,
              help='BWA可执行文件路径|BWA executable path')
@click.option('--samtools-bin',
              default="~/miniforge3/envs/GATK_v.4.6.2.0/bin/samtools",
              show_default=True,
              help='Samtools可执行文件路径|Samtools executable path')

# BWA比对配置|BWA alignment configuration
@click.option('--samblaster-bin',
              default="~/miniforge3/envs/Population_genetics/bin/samblaster",
              show_default=True,
              help='Samblaster可执行文件路径|Samblaster executable path')
@click.option('--haphic-filter-bam-bin',
              default="~/miniforge3/envs/haphic/bin/filter_bam",
              show_default=True,
              help='HapHiC filter_bam工具路径|HapHiC filter_bam tool path')
@click.option('--use-samblaster/--no-use-samblaster', default=True,
              help='使用samblaster去重|Use samblaster deduplication')
@click.option('--use-haphic-filter/--no-use-haphic-filter', default=True,
              help='使用HapHiC过滤|Use HapHiC filtering')

# Juicebox配置|Juicebox configuration
@click.option('--generate-juicebox/--no-generate-juicebox', default=True,
              help='生成Juicebox兼容文件|Generate Juicebox compatible files')
@click.option('--matlock-bin',
              default="matlock",
              show_default=True,
              help='Matlock可执行文件路径|Matlock executable path')
@click.option('--three-d-dna-dir',
              default="~/software/3d-dna",
              show_default=True,
              help='3D-DNA目录路径|3D-DNA directory path')
@click.option('--agp2assembly-script',
              default="~/software/3d-dna/utils/agp2assembly.py",
              show_default=True,
              help='agp2assembly脚本路径|agp2assembly script path')
@click.option('--asm-visualizer-script',
              default="~/software/3d-dna/visualize/run-assembly-visualizer.sh",
              show_default=True,
              help='asm-visualizer脚本路径|asm-visualizer script path')

# 高级选项|Advanced options
@click.option('--RE', default="GATC",
              show_default=True,
              help='限制性酶切位点|Restriction enzyme sites')
@click.option('--quick-view', is_flag=True,
              help='快速查看模式|Quick view mode')

# 输出格式选项|Output format options
@click.option('--no-agp', is_flag=True,
              help='不输出AGP文件|Don\'t output AGP file')
@click.option('--no-fasta', is_flag=True,
              help='不输出FASTA文件|Don\'t output FASTA file')
@click.option('--no-generate-plots', is_flag=True,
              help='不生成可视化图表|Don\'t generate visualization plots')
@click.option('--no-juicebox', is_flag=True,
              help='不生成Juicebox脚本|Don\'t generate Juicebox script')

# 日志配置|Logging configuration
@click.option('--verbose', '-v', is_flag=True,
              help='详细输出模式|Verbose output mode')
@click.option('--log-file',
              type=click.Path(file_okay=True, dir_okay=False, resolve_path=True),
              help='日志文件路径|Log file path')

@click.option('--dry-run', is_flag=True,
              help='测试模式,不执行实际命令|Test mode, do not execute actual commands')
def haphic(input, bam, hic1, hic2, chr_number, output_dir, prefix, force_rerun,
           mapq_threshold, edit_distance, re_site_cutoff, min_re_sites, aln_format,
           min_inflation, max_inflation, inflation_step, nx, min_group_len, flank, bin_size_kbp,
           skip_fast_sort, skip_allhic, skip_ga, sort_by_input, no_additional_rescue,
           remove_concentrated_links, normalize_by_nlinks, dense_matrix,
           processes, threads, memory_limit,
           correct_nrounds, correct_min_coverage, median_cov_ratio, region_len_ratio, min_region_cutoff,
           remove_allelic_links, phasing_weight, gfa_files, re,
           bin_size, min_len, generate_plots, separate_plots,
           haphic_bin, bwa_bin, samtools_bin,
           samblaster_bin, haphic_filter_bam_bin, use_samblaster, use_haphic_filter,
           generate_juicebox, matlock_bin, three_d_dna_dir, agp2assembly_script, asm_visualizer_script,
           quick_view,
           no_agp, no_fasta, no_generate_plots, no_juicebox,
           verbose, log_file, dry_run):
    """
    HapHiC基因组支架流程|HapHiC Genome Scaffolding Pipeline

    基于Hi-C数据的基因组支架构建自动化流程
    Automated genome scaffolding pipeline based on Hi-C data

    示例|Examples: biopytools haphic -i assembly.fa -b hic.bam -c 24
    """

    try:
        # 参数验证|Parameter validation
        input_count = sum(1 for x in [bam, hic1, hic2] if x is not None)

        if input_count == 0:
            raise click.BadParameter("必须提供输入: 指定--bam或--hic1/--hic2|Must provide input: --bam or --hic1/--hic2")
        elif bam and (hic1 or hic2):
            raise click.BadParameter("不能同时使用--bam和--hic1/--hic2|Cannot use both --bam and --hic1/--hic2")
        elif (hic1 and not hic2) or (hic2 and not hic1):
            raise click.BadParameter("FASTQ模式需要同时提供--hic1和--hic2|FASTQ mode requires both --hic1 and --hic2")

        # 确定输入类型|Determine input type
        if bam:
            hic_file = bam
            hic_file_type = "bam"
        else:
            hic_file = hic1  # 使用hic1作为基准|Use hic1 as reference
            hic_file_type = "fastq"

        # 延迟加载|Lazy loading
        HapHiCProcessor = _lazy_import_haphic_main()

        # 创建处理器|Create processor
        processor = HapHiCProcessor(
            asm_file=str(input),
            hic_file=str(hic_file),
            nchrs=chr_number,
            hic_file_type=hic_file_type,
            haphic_bin=haphic_bin,
            output_dir=str(output_dir) if output_dir else None,
            prefix=prefix,
            force_rerun=force_rerun,
            mapq_threshold=mapq_threshold,
            edit_distance=edit_distance,
            re_site_cutoff=re_site_cutoff,
            min_re_sites=min_re_sites,
            aln_format=aln_format,
            min_inflation=min_inflation,
            max_inflation=max_inflation,
            inflation_step=inflation_step,
            nx=nx,
            min_group_len=min_group_len,
            flank=flank,
            bin_size_kbp=bin_size_kbp,
            skip_fast_sort=skip_fast_sort,
            skip_allhic=skip_allhic,
            skip_ga=skip_ga,
            sort_by_input=sort_by_input,
            no_additional_rescue=no_additional_rescue,
            remove_concentrated_links=remove_concentrated_links,
            normalize_by_nlinks=normalize_by_nlinks,
            dense_matrix=dense_matrix,
            processes=processes,
            threads=threads,
            correct_nrounds=correct_nrounds,
            correct_min_coverage=correct_min_coverage,
            median_cov_ratio=median_cov_ratio,
            region_len_ratio=region_len_ratio,
            min_region_cutoff=min_region_cutoff,
            remove_allelic_links=remove_allelic_links,
            phasing_weight=phasing_weight,
            gfa_files=gfa_files,
            re_sites=re,
            bin_size=bin_size,
            min_len=min_len,
            quick_view=quick_view,
            generate_plots=generate_plots and not no_generate_plots,
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

        # FASTQ模式需要设置hic2|FASTQ mode requires setting hic2
        if hic_file_type == "fastq":
            processor.config.hic2_file = str(hic2)

        # 测试运行|Dry run
        if dry_run:
            processor.config.dry_run = True
            click.echo("测试模式|Test mode: HapHiC")

        # 运行分析|Run analysis
        success = processor.run_pipeline()

        if success:
            click.echo("HapHiC scaffolding完成|HapHiC scaffolding completed!")
        else:
            click.echo("HapHiC scaffolding失败|HapHiC scaffolding failed!", err=True)
            sys.exit(1)

    except KeyboardInterrupt:
        click.echo("\n操作被用户中断|Operation interrupted by user", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"发生错误|Error occurred: {e}", err=True)
        sys.exit(1)