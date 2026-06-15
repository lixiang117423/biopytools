"""Hi-C热图分析CLI包装器|Hi-C heatmap analysis CLI wrapper (HiCPro + PlotHiC)"""

import click
import sys
import os


def _lazy_import_main():
    """延迟加载主函数|Lazy load main function"""
    try:
        from biopytools.hic_heatmap.main import main as hic_heatmap_main
        return hic_heatmap_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_path_exists(path):
    """验证路径存在|Validate path exists"""
    if not _is_help_request() and path and not os.path.exists(path):
        raise click.BadParameter(f"路径不存在|Path does not exist: {path}")
    return path


@click.command(short_help="Hi-C热图分析 (HiCPro + PlotHiC)|Hi-C heatmap analysis (HiCPro + PlotHiC)")
@click.option('-i', '--genome',
              required=True,
              callback=lambda ctx, param, value: _validate_path_exists(value),
              help='基因组FASTA文件|Genome FASTA file')
@click.option('-g', '--genome-id',
              required=True,
              help='基因组ID（用于输出文件命名，如hg19, mm10）|Genome ID (for output file naming, e.g., hg19, mm10)')
@click.option('-o', '--output-dir',
              default='./hic_output',
              show_default=True,
              help='输出目录|Output directory')
@click.option('-1', '--fastq-r1',
              required=True,
              callback=lambda ctx, param, value: _validate_path_exists(value),
              help='R1测序文件|R1 sequencing file')
@click.option('-2', '--fastq-r2',
              required=True,
              callback=lambda ctx, param, value: _validate_path_exists(value),
              help='R2测序文件|R2 sequencing file')
@click.option('-t', '--threads',
              default=64,
              show_default=True,
              help='线程数|Threads')
@click.option('--max-memory',
              default=200,
              type=int,
              show_default=True,
              help='HiC-Pro最大内存限制（GB）|Maximum memory limit for HiC-Pro in GB')
@click.option('--restriction-enzyme',
              default='MboI',
              show_default=True,
              help='限制性内切酶|Restriction enzyme (default: MboI). Options: MboI, HindIII, NcoI, EcoRI, BamHI')
@click.option('--bowtie2-idx',
              help='Bowtie2索引路径（默认自动生成）|Bowtie2 index path (auto-generated if not specified)')
@click.option('--bin-sizes',
              default='20000 40000 150000 500000 1000000',
              show_default=True,
              help='Contact map bin大小（空格分隔）|Contact map bin sizes, space-separated')
@click.option('--resolution',
              type=int,
              default=100000,
              show_default=True,
              help='热图分辨率（bp）|Heatmap resolution in bp (default: 100000, 100kb)')
@click.option('--color-map',
              default='YlOrRd',
              show_default=True,
              help='颜色方案|Color scheme (PlotHiC default: YlOrRd)')
@click.option('--dpi',
              type=int,
              default=300,
              show_default=True,
              help='图像分辨率|Image DPI')
@click.option('--format',
              default='pdf',
              show_default=True,
              help='输出格式|Output format (pdf, png, svg, etc.)')
@click.option('--bar-max',
              type=int,
              default=1,
              show_default=True,
              help='颜色条最大值|Color bar maximum value (after log transform)')
@click.option('--hicpro-sif',
              default='',
              help='HiCPro singularity镜像路径（留空则直接使用HiC-Pro）|HiCPro singularity image path (leave empty to use HiC-Pro directly)')
@click.option('--plothic-path',
              default='~/miniforge3/envs/plothic_v.1.0.0/bin/plothic',
              show_default=True,
              help='PlotHiC可执行文件路径|PlotHiC executable path')
@click.option('--force',
              is_flag=True,
              help='强制重新运行|Force rerun all steps')
@click.option('--verbose',
              is_flag=True,
              help='显示详细日志|Verbose logging')
@click.option('--quiet',
              is_flag=True,
              help='仅显示错误|Errors only')
def hic_heatmap(genome, genome_id, output_dir, fastq_r1, fastq_r2,
                threads, max_memory, restriction_enzyme, bowtie2_idx, bin_sizes,
                resolution, color_map, dpi, format, bar_max,
                hicpro_sif, plothic_path, force, verbose, quiet):
    """
    Hi-C全基因组热图分析流程 (HiCPro + PlotHiC)|Hi-C whole genome heatmap analysis pipeline (HiCPro + PlotHiC)

    使用HiCPro进行比对和矩阵生成，使用PlotHiC绘制全基因组热图|Use HiCPro for alignment and matrix generation, PlotHiC for plotting.

    示例|Example: biopytools hic-heatmap -i genome.fa -g EcA -1 R1.fq.gz -2 R2.fq.gz -o output
    """
    hic_heatmap_main = _lazy_import_main()

    # 构建参数列表|Build argument list
    args = ['hic_heatmap.py']
    args.extend(['-i', genome])
    args.extend(['-g', genome_id])
    args.extend(['-o', output_dir])
    args.extend(['-1', fastq_r1])
    args.extend(['-2', fastq_r2])
    args.extend(['-t', str(threads)])
    args.extend(['--max-memory', str(max_memory)])

    # HiCPro参数|HiCPro parameters
    args.extend(['--restriction-enzyme', restriction_enzyme])
    if bowtie2_idx:
        args.extend(['--bowtie2-idx', bowtie2_idx])
    args.extend(['--bin-sizes', bin_sizes])

    # PlotHiC参数|PlotHiC parameters
    args.extend(['--resolution', str(resolution)])
    args.extend(['--color-map', color_map])
    args.extend(['--dpi', str(dpi)])
    args.extend(['--format', format])
    args.extend(['--bar-max', str(bar_max)])

    # 工具路径|Tool paths
    args.extend(['--hicpro-sif', hicpro_sif])
    args.extend(['--plothic-path', plothic_path])

    # 流程控制|Process control
    if force:
        args.append('--force')
    if verbose:
        args.append('--verbose')
    if quiet:
        args.append('--quiet')

    original_argv = sys.argv
    sys.argv = args

    try:
        hic_heatmap_main()
    except SystemExit as e:
        sys.exit(e.code)
    finally:
        sys.argv = original_argv
