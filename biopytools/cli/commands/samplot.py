"""
Samplot命令|Samplot Command
"""

import click
import sys
import os


def _lazy_import_plot_main():
    """延迟加载samplot plot主函数|Lazy load samplot plot main function"""
    try:
        from ...samplot.main import main_plot
        return main_plot
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _lazy_import_vcf_main():
    """延迟加载samplot vcf主函数|Lazy load samplot vcf main function"""
    try:
        from ...samplot.main import main_vcf
        return main_vcf
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_file_exists(file_path):
    """验证文件存在|Validate file existence"""
    if _is_help_request():
        return file_path
    if file_path and not os.path.exists(file_path):
        raise click.BadParameter(f"文件不存在|File does not exist: {file_path}")
    return file_path


@click.group(short_help='Samplot结构变异可视化|Samplot SV visualization',
             context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120),
             invoke_without_command=True)
@click.pass_context
def samplot(ctx):
    """Samplot结构变异可视化工具|Samplot SV visualization tool

    基因组结构变异(SV)可视化，支持单区域绘图和VCF批量绘图
    |Genomic structural variant visualization, supports single region and VCF batch plotting

    子命令|Subcommands: plot, vcf

    示例|Examples: biopytools samplot plot -b sample.bam -c chr1 -s 1000 -e 5000
    """
    if ctx.invoked_subcommand is None:
        click.echo(ctx.get_help())


@samplot.command('plot', short_help='绘制单区域SV图|Plot single SV region',
                 context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('-b', '--bams',
              required=True,
              multiple=True,
              callback=lambda ctx, param, value: [_validate_file_exists(f) for f in value] if value else None,
              help='BAM/CRAM文件路径|BAM/CRAM file paths')
@click.option('-c', '--chrom',
              required=True,
              help='染色体名称|Chromosome name')
@click.option('-s', '--start',
              required=True,
              type=int,
              help='起始位置|Start position')
@click.option('-e', '--end',
              required=True,
              type=int,
              help='结束位置|End position')
@click.option('-t', '--sv-type',
              default=None,
              help='SV类型(DEL/DUP/INV/BND)|SV type')
@click.option('-o', '--output-file',
              default=None,
              help='输出文件名|Output file name')
@click.option('--output-dir',
              default='.',
              show_default=True,
              help='输出目录|Output directory')
@click.option('-r', '--reference',
              default=None,
              help='参考基因组(CRAM必须)|Reference genome (required for CRAM)')
@click.option('-d', '--max-depth',
              type=int,
              default=1,
              show_default=True,
              help='最大正常pair数|Max normal pairs to plot')
@click.option('-w', '--window',
              type=int,
              default=None,
              help='窗口大小|Window size')
@click.option('-z', '--z',
              type=int,
              default=4,
              show_default=True,
              help='标准差倍数|Number of stdevs from mean')
@click.option('-H', '--plot-height',
              type=int,
              default=None,
              help='图高|Plot height')
@click.option('-W', '--plot-width',
              type=int,
              default=8,
              show_default=True,
              help='图宽|Plot width')
@click.option('--dpi',
              type=int,
              default=300,
              show_default=True,
              help='DPI|Dots per inch')
@click.option('--long-read',
              type=int,
              default=1000,
              show_default=True,
              help='长读长最小长度|Min length for long-read')
@click.option('--coverage-only',
              is_flag=True,
              default=False,
              help='仅显示覆盖度|Show only coverage')
@click.option('--same-yaxis-scales',
              is_flag=True,
              default=False,
              help='统一Y轴|Use same Y-axis scales')
@click.option('-n', '--titles',
              multiple=True,
              default=None,
              help='样本标题|Sample titles')
@click.option('--samplot-path',
              default='~/miniforge3/envs/samplot_v.1.3.0/bin/samplot',
              show_default=True,
              help='samplot路径|samplot binary path')
def plot(bams, chrom, start, end, sv_type, output_file, output_dir, reference,
         max_depth, window, z, plot_height, plot_width, dpi, long_read,
         coverage_only, same_yaxis_scales, titles, samplot_path):
    """
    绘制单区域结构变异图|Plot single SV region

    示例|Examples: biopytools samplot plot -b sample.bam -c chr1 -s 1000 -e 5000 -t DEL
    """
    main_plot = _lazy_import_plot_main()

    args = ['samplot_plot.py']
    args.extend(['-b'] + list(bams))
    args.extend(['-c', chrom])
    args.extend(['-s', str(start)])
    args.extend(['-e', str(end)])

    if sv_type:
        args.extend(['-t', sv_type])
    if output_file:
        args.extend(['-o', output_file])
    if output_dir != '.':
        args.extend(['--output-dir', output_dir])
    if reference:
        args.extend(['-r', reference])
    if max_depth != 1:
        args.extend(['-d', str(max_depth)])
    if window is not None:
        args.extend(['-w', str(window)])
    if z != 4:
        args.extend(['-z', str(z)])
    if plot_height is not None:
        args.extend(['-H', str(plot_height)])
    if plot_width != 8:
        args.extend(['-W', str(plot_width)])
    if dpi != 300:
        args.extend(['--dpi', str(dpi)])
    if long_read != 1000:
        args.extend(['--long-read', str(long_read)])
    if coverage_only:
        args.append('--coverage-only')
    if same_yaxis_scales:
        args.append('--same-yaxis-scales')
    if titles:
        args.extend(['-n'] + list(titles))
    if samplot_path != '~/miniforge3/envs/samplot_v.1.3.0/bin/samplot':
        args.extend(['--samplot-path', samplot_path])

    original_argv = sys.argv
    sys.argv = args

    try:
        main_plot()
    except SystemExit as e:
        sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n用户中断|User interrupted", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv


@samplot.command('vcf', short_help='批量绘制VCF中的SV|Batch plot SVs from VCF',
                 context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('-b', '--bams',
              required=True,
              multiple=True,
              callback=lambda ctx, param, value: [_validate_file_exists(f) for f in value] if value else None,
              help='BAM/CRAM文件路径|BAM/CRAM file paths')
@click.option('--vcf',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='VCF文件路径|VCF file path')
@click.option('-d', '--output-dir',
              default='samplot-out',
              show_default=True,
              help='输出目录|Output directory')
@click.option('-O', '--output-type',
              default='png',
              show_default=True,
              type=click.Choice(['png', 'pdf', 'eps', 'jpg']),
              help='输出格式|Output format')
@click.option('-t', '--threads',
              type=int,
              default=1,
              show_default=True,
              help='线程数|Number of threads')
@click.option('--downsample',
              type=int,
              default=1,
              show_default=True,
              help='下采样数|Downsample count')
@click.option('--min-bp',
              type=int,
              default=20,
              show_default=True,
              help='最小SV长度(bp)|Min SV length in bp')
@click.option('--max-mb',
              type=int,
              default=None,
              help='最大SV长度(MB)|Max SV length in MB')
@click.option('--sample-ids',
              multiple=True,
              default=None,
              help='样本ID列表|Sample ID list')
@click.option('--plot-all',
              is_flag=True,
              default=False,
              help='绘制所有样本|Plot all samples')
@click.option('--min-call-rate',
              type=float,
              default=None,
              help='最小call rate|Min call rate')
@click.option('--max-hets',
              type=int,
              default=None,
              help='最大杂合数|Max heterozygotes')
@click.option('--min-entries',
              type=int,
              default=6,
              show_default=True,
              help='最小样本数|Min entries to plot')
@click.option('--max-entries',
              type=int,
              default=10,
              show_default=True,
              help='最大样本数|Max entries to plot')
@click.option('--gff3',
              default=None,
              help='GFF3注释文件|GFF3 annotation file')
@click.option('--samplot-path',
              default='~/miniforge3/envs/samplot_v.1.3.0/bin/samplot',
              show_default=True,
              help='samplot路径|samplot binary path')
def vcf(bams, vcf, output_dir, output_type, threads, downsample, min_bp, max_mb,
        sample_ids, plot_all, min_call_rate, max_hets, min_entries, max_entries,
        gff3, samplot_path):
    """
    批量绘制VCF中的结构变异|Batch plot SVs from VCF

    示例|Examples: biopytools samplot vcf --vcf variants.vcf -b sample.bam -d output/
    """
    main_vcf = _lazy_import_vcf_main()

    args = ['samplot_vcf.py']
    args.extend(['-b'] + list(bams))
    args.extend(['--vcf', vcf])

    if output_dir != 'samplot-out':
        args.extend(['-d', output_dir])
    if output_type != 'png':
        args.extend(['-O', output_type])
    if threads != 1:
        args.extend(['-t', str(threads)])
    if downsample != 1:
        args.extend(['--downsample', str(downsample)])
    if min_bp != 20:
        args.extend(['--min-bp', str(min_bp)])
    if max_mb is not None:
        args.extend(['--max-mb', str(max_mb)])
    if sample_ids:
        args.extend(['--sample-ids'] + list(sample_ids))
    if plot_all:
        args.append('--plot-all')
    if min_call_rate is not None:
        args.extend(['--min-call-rate', str(min_call_rate)])
    if max_hets is not None:
        args.extend(['--max-hets', str(max_hets)])
    if min_entries != 6:
        args.extend(['--min-entries', str(min_entries)])
    if max_entries != 10:
        args.extend(['--max-entries', str(max_entries)])
    if gff3:
        args.extend(['--gff3', gff3])
    if samplot_path != '~/miniforge3/envs/samplot_v.1.3.0/bin/samplot':
        args.extend(['--samplot-path', samplot_path])

    original_argv = sys.argv
    sys.argv = args

    try:
        main_vcf()
    except SystemExit as e:
        sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n用户中断|User interrupted", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
