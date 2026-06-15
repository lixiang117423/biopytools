"""
PGGB泛基因组图构建命令|PGGB Pangenome Graph Builder Command
"""

import click
import sys
import os


def _lazy_import_pggb_main():
    """延迟加载pggb主函数|Lazy load pggb main function"""
    try:
        from ...pggb.main import main as pggb_main
        return pggb_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_file_exists(file_path):
    """验证文件存在(仅在非帮助模式)|Validate file exists (only in non-help mode)"""
    if not _is_help_request() and file_path and not os.path.exists(file_path):
        raise click.BadParameter(f"文件不存在|File does not exist: {file_path}")
    return file_path


@click.command(short_help="PGGB泛基因组图构建|PGGB Pangenome Graph Builder",
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('-i', '--input', required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='输入FASTA文件|Input FASTA file')
@click.option('-o', '--output', required=True,
              help='输出目录|Output directory')
@click.option('-t', '--threads', type=int, default=24,
              help='线程数|Number of threads (default: 24)')
@click.option('--conda-env', default='pggb_v.0.7.4',
              help='conda环境名|Conda environment name (default: pggb_v.0.7.4)')
@click.option('-s', '--segment-length', type=int, default=5000,
              help='比对分段长度|Segment length (default: 5000)')
@click.option('-p', '--map-pct-id', type=int, default=90,
              help='比对一致度|Map percent identity (default: 90)')
@click.option('-n', '--n-haplotypes', type=int, default=0,
              help='单倍型数(0=自动)|N haplotypes (0=auto)')
@click.option('--vcf-spec', default='',
              help='VCF输出参考规范|VCF output reference spec')
@click.option('--resume', is_flag=True,
              help='断点续传|Resume from existing outputs')
@click.option('--compress', is_flag=True,
              help='压缩输出|Compress output files')
@click.option('--stats', is_flag=True,
              help='生成统计信息|Generate statistics')
def pggb(**kwargs):
    """
    PGGB泛基因组图构建工具|PGGB Pangenome Graph Builder

    使用wfmash+seqwish+smoothxg构建泛基因组变异图
    Build pangenome variation graph using wfmash+seqwish+smoothxg

    示例|Examples: biopytools pggb -i genomes.fa -o output/
    """
    pggb_main = _lazy_import_pggb_main()

    # 构建参数列表|Build arguments list
    argv = ['pggb.py']

    argv.extend(['-i', kwargs['input']])
    argv.extend(['-o', kwargs['output']])
    argv.extend(['-t', str(kwargs['threads'])])
    argv.extend(['--conda-env', kwargs['conda_env']])
    argv.extend(['-s', str(kwargs['segment_length'])])
    argv.extend(['-p', str(kwargs['map_pct_id'])])

    if kwargs['n_haplotypes'] > 0:
        argv.extend(['-n', str(kwargs['n_haplotypes'])])

    if kwargs['vcf_spec']:
        argv.extend(['--vcf-spec', kwargs['vcf_spec']])

    if kwargs['resume']:
        argv.append('--resume')

    if kwargs['compress']:
        argv.append('--compress')

    if kwargs['stats']:
        argv.append('--stats')

    original_argv = sys.argv
    sys.argv = argv

    try:
        pggb_main()
    except SystemExit as e:
        sys.exit(e.code)
    finally:
        sys.argv = original_argv
