"""
GenomeSyn2比较基因组学可视化命令|GenomeSyn2 Comparative Genomics Visualization Command
"""

import click
import sys
import os


def _lazy_import_genomesyn2_main():
    """延迟加载genomesyn2主函数|Lazy load genomesyn2 main function"""
    try:
        from ...genomesyn2.main import main as genomesyn2_main
        return genomesyn2_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_path_exists(path):
    """验证路径存在(仅在非帮助模式)|Validate path exists (only in non-help mode)"""
    if path and not _is_help_request() and not os.path.exists(path):
        raise click.BadParameter(f"路径不存在|Path does not exist: {path}")
    return path


@click.command(
    short_help='GenomeSyn2比较基因组学可视化工具|GenomeSyn2 comparative genomics visualization',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)

# ========== 比对模式参数 | Alignment Mode Parameters ==========
@click.option('--align',
              type=click.Choice(['mummer', 'minimap2', 'blastp', 'mmseqs', 'diamond']),
              help='比对软件类型|Alignment software type')
@click.option('--genome',
              callback=lambda ctx, param, value: _validate_path_exists(value) if value else None,
              help='基因组文件目录|Genome files directory')
@click.option('--gene',
              callback=lambda ctx, param, value: _validate_path_exists(value) if value else None,
              help='基因注释文件目录|Gene annotation files directory')
@click.option('--outdir',
              help='输出目录|Output directory')

# ========== VCF模式参数 | VCF Mode Parameters ==========
@click.option('--vcf',
              callback=lambda ctx, param, value: _validate_path_exists(value) if value else None,
              help='VCF文件路径|VCF file path for SNP analysis')
@click.option('--bin',
              type=int,
              default=50000,
              show_default=True,
              help='Bin大小(用于SNP分析)|Bin size for SNP analysis')
@click.option('--identity',
              callback=lambda ctx, param, value: _validate_path_exists(value) if value else None,
              help='SNP一致性文件|SNP identity BED file')
@click.option('--density',
              callback=lambda ctx, param, value: _validate_path_exists(value) if value else None,
              help='SNP密度文件|SNP density BED file')

# ========== 绘图模式参数 | Plotting Mode Parameters ==========
@click.option('--conf',
              callback=lambda ctx, param, value: _validate_path_exists(value) if value else None,
              help='配置文件路径|Configuration file path')
@click.option('--anno',
              is_flag=True,
              help='显示注释配置选项|Show annotation configuration options')

# ========== 文件生成模式参数 | File Generation Mode Parameters ==========
@click.option('--type',
              type=click.Choice(['fa', 'prot', 'anno']),
              help='文件类型(用于生成文件列表)|File type for generating file list')
@click.option('--path',
              callback=lambda ctx, param, value: _validate_path_exists(value) if value else None,
              help='文件路径|File path for generating list')
@click.option('--out',
              help='输出文件名|Output file name')

# ========== 通用参数 | Common Parameters ==========
@click.option('-t', '--threads',
              type=int,
              default=12,
              show_default=True,
              help='线程数|Number of threads')

def genomesyn2(align, genome, gene, outdir, vcf, bin, identity, density,
               conf, anno, type, path, out, threads):
    """
    GenomeSyn2比较基因组学可视化工具|GenomeSyn2 Comparative Genomics Visualization Tool

    示例|Examples: biopytools genomesyn2 --align mummer --genome ./genome_dir/ --outdir ./output/ -t 12
    """

    # 延迟加载|Lazy loading
    genomesyn2_main = _lazy_import_genomesyn2_main()

    # 构建参数列表|Build argument list
    args = ['genomesyn2.py']

    # 比对模式参数|Alignment mode parameters
    if align:
        args.extend(['--align', align])
    if genome:
        args.extend(['--genome', genome])
    if gene:
        args.extend(['--gene', gene])
    if outdir:
        args.extend(['--outdir', outdir])

    # VCF模式参数|VCF mode parameters
    if vcf:
        args.extend(['--vcf', vcf])
    if bin != 50000:
        args.extend(['--bin', str(bin)])
    if identity:
        args.extend(['--identity', identity])
    if density:
        args.extend(['--density', density])

    # 绘图模式参数|Plotting mode parameters
    if conf:
        args.extend(['--conf', conf])
    if anno:
        args.append('--anno')

    # 文件生成模式参数|File generation mode parameters
    if type:
        args.extend(['--type', type])
    if path:
        args.extend(['--path', path])
    if out:
        args.extend(['--out', out])

    # 通用参数|Common parameters
    if threads != 12:
        args.extend(['--threads', str(threads)])

    # 执行主程序|Execute main program
    original_argv = sys.argv
    sys.argv = args

    try:
        genomesyn2_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
