"""SNP区域基因提取命令|SNP Region Gene Extractor Command"""

import click
import sys
import os


def _lazy_import_snp_region_gene_main():
    """延迟加载snp_region_gene主函数|Lazy load snp_region_gene main function"""
    try:
        from ...snp_region_gene.main import main as snp_region_gene_main
        return snp_region_gene_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_file_exists(file_path):
    """验证文件存在(仅在非帮助模式)|Validate file existence (only in non-help mode)"""
    if not _is_help_request() and file_path and not os.path.exists(file_path):
        raise click.BadParameter(f"文件不存在|File does not exist: {file_path}")
    return file_path


@click.command(
    short_help='SNP区域基因提取工具|SNP region gene extractor',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--snp', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='SNP位置文件（格式：Chr01:24770）|SNP position file (format: Chr01:24770)')
@click.option('--gff', '-g',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='GFF3注释文件|GFF3 annotation file')
@click.option('--genome', '-G',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='基因组FASTA文件|Genome FASTA file')
@click.option('--left', '-l',
              type=int,
              default=0,
              show_default=True,
              help='SNP上游距离（bp）|Upstream distance from SNP (bp)')
@click.option('--right', '-r',
              type=int,
              default=0,
              show_default=True,
              help='SNP下游距离（bp）|Downstream distance from SNP (bp)')
@click.option('--promoter', '-p',
              type=int,
              default=2000,
              show_default=True,
              help='启动子区域距离（bp）|Promoter region distance (bp)')
@click.option('--output', '-o',
              default='./snp_region_output',
              show_default=True,
              help='输出文件前缀|Output file prefix')
@click.option('--gffread-path',
              default='gffread',
              show_default=True,
              help='gffread程序路径|gffread program path')
@click.option('--seqkit-path',
              default='seqkit',
              show_default=True,
              help='seqkit程序路径|seqkit program path')
@click.option('--keep-temp',
              is_flag=True,
              help='保留临时文件|Keep temporary files')
def snp_region_gene(snp, gff, genome, left, right, promoter, output,
                    gffread_path, seqkit_path, keep_temp):
    """
    SNP区域基因提取工具|SNP Region Gene Extractor Tool

    根据SNP位置和指定区域，提取相关基因的CDS和蛋白序列|Extract CDS and protein sequences of genes in specified regions around SNPs

    示例|Examples: biopytools snp-region-gene -i snp.txt -g annotation.gff3 -G genome.fa -l 100000 -r 100000 -o output
    """

    # 延迟加载|Lazy loading
    snp_region_gene_main = _lazy_import_snp_region_gene_main()

    # 构建参数列表|Build argument list
    args = ['snp_region_gene.py']

    # 必需参数|Required parameters
    args.extend(['-i', snp])
    args.extend(['-g', gff])
    args.extend(['-G', genome])

    # 可选参数|Optional parameters
    if left != 0:
        args.extend(['-l', str(left)])

    if right != 0:
        args.extend(['-r', str(right)])

    if promoter != 2000:
        args.extend(['-p', str(promoter)])

    if output != './snp_region_output':
        args.extend(['-o', output])

    if gffread_path != 'gffread':
        args.extend(['--gffread-path', gffread_path])

    if seqkit_path != 'seqkit':
        args.extend(['--seqkit-path', seqkit_path])

    if keep_temp:
        args.append('--keep-temp')

    # 执行主程序|Execute main program
    original_argv = sys.argv
    sys.argv = args

    try:
        snp_region_gene_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
