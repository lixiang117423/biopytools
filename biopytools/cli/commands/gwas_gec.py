"""
GWAS基因组范围多重检验校正命令|GWAS Genome-wide Error Correction Command

"""

import click
import sys
import os


def _lazy_import_gwas_gec_main():
    """延迟加载gwas_gec主函数|Lazy load gwas_gec main function"""
    try:
        from ...gwas_gec.main import main as gwas_gec_main
        return gwas_gec_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_file_exists(file_path):
    """验证文件存在(仅在非帮助模式)|Validate file existence (only in non-help mode)"""
    if not _is_help_request() and not os.path.exists(file_path):
        raise click.BadParameter(f"文件不存在|File does not exist: {file_path}")
    return file_path


@click.command(
    short_help='GWAS基因组范围多重检验校正|GWAS genome-wide error correction',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--pfile', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='GWAS P值汇总统计文件|GWAS P-value summary statistics file')
@click.option('--reference', '-r',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='参考文件（VCF或PLINK binary前缀）|Reference file (VCF or PLINK binary prefix)')
@click.option('--output-dir', '-o',
              default='./gec_output',
              show_default=True,
              type=click.Path(),
              help='输出目录|Output directory')
@click.option('--threads', '-t',
              default=12,
              show_default=True,
              type=int,
              help='线程数|Number of threads')
@click.option('--memory', '-m',
              default='100G',
              show_default=True,
              help='Java内存分配|Java memory allocation')
@click.option('--maf-filter',
              default=0.05,
              show_default=True,
              type=float,
              help='MAF过滤阈值|MAF filter threshold')
@click.option('--p-cutoff',
              default=0.05,
              show_default=True,
              type=float,
              help='P值阈值|P-value threshold')
@click.option('--no-keep-ref',
              is_flag=True,
              help='不保留参考文件缓存|Do not keep reference file cache')
@click.option('--no-convert-chrom',
              is_flag=True,
              help='禁用自动染色体格式转换|Disable automatic chromosome format conversion')
@click.option('--chrom',
              help='染色体范围(如1-22或1-10,22)|Chromosome range (e.g., 1-22 or 1-10,22)')
@click.option('--chrom-col',
              default='CHR',
              show_default=True,
              help='染色体列名|Chromosome column name')
@click.option('--pos-col',
              default='BP',
              show_default=True,
              help='位置列名|Position column name')
@click.option('--p-col',
              default='P',
              show_default=True,
              help='P值列名|P-value column name')
@click.option('--alpha',
              default=0.05,
              show_default=True,
              type=float,
              help='显著性水平(FWER)|Significance level')
@click.option('--kggsee-jar',
              default='~/software/kmmsee/kggsee.jar',
              show_default=True,
              help='KGGSee JAR文件路径|KGGSee JAR file path')
def gwas_gec(pfile, reference, output_dir, threads, memory, maf_filter, p_cutoff,
             no_keep_ref, no_convert_chrom, chrom, chrom_col, pos_col, p_col, alpha, kggsee_jar):
    """
    GWAS基因组范围多重检验校正工具|GWAS Genome-wide Error Correction Tool

    使用GEC算法基于LD结构计算GWAS显著性阈值|Calculate GWAS significance thresholds based on LD structure using GEC algorithm
    支持VCF格式的参考文件|Supports VCF format reference file
    自动转换染色体格式(Chr01 -> 1)|Auto-convert chromosome format (Chr01 -> 1)

    示例|Examples: biopytools gwas-gec -i gwas.txt -r input.vcf.gz    # 使用VCF文件|Use VCF file
    """

    # 延迟加载|Lazy loading
    gwas_gec_main = _lazy_import_gwas_gec_main()

    # 构建参数列表|Build argument list
    args = ['gwas_gec.py']

    # 必需参数|Required parameters
    args.extend(['--pfile', pfile])
    args.extend(['--reference', reference])

    # 可选参数|Optional parameters
    if output_dir != './gec_output':
        args.extend(['--output-dir', output_dir])

    if threads != 16:
        args.extend(['--threads', str(threads)])

    if memory != '8g':
        args.extend(['--memory', memory])

    if maf_filter != 0.05:
        args.extend(['--maf-filter', str(maf_filter)])

    if p_cutoff != 0.05:
        args.extend(['--p-cutoff', str(p_cutoff)])

    if no_keep_ref:
        args.append('--no-keep-ref')

    if no_convert_chrom:
        args.append('--no-convert-chrom')

    if chrom:
        args.extend(['--chrom', chrom])

    if chrom_col != 'CHR':
        args.extend(['--chrom-col', chrom_col])

    if pos_col != 'BP':
        args.extend(['--pos-col', pos_col])

    if p_col != 'P':
        args.extend(['--p-col', p_col])

    if alpha != 0.05:
        args.extend(['--alpha', str(alpha)])

    if kggsee_jar != '~/software/kmmsee/kggsee.jar':
        args.extend(['--kggsee-jar', kggsee_jar])

    # 执行主程序|Execute main program
    original_argv = sys.argv
    sys.argv = args

    try:
        gwas_gec_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
