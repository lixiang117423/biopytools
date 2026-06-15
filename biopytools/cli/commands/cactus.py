"""
Cactus泛基因组分析命令|Cactus Pangenome Analysis Command
"""

import click
import sys
import os


def _lazy_import_cactus_main():
    """延迟加载cactus主函数|Lazy load cactus main function"""
    try:
        from ...cactus.cli import main as cactus_main
        return cactus_main
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


def _validate_directory_exists(dir_path):
    """验证目录存在(如不存在则创建)|Validate directory exists (create if not exists)"""
    if not _is_help_request() and dir_path:
        if not os.path.exists(dir_path):
            try:
                os.makedirs(dir_path, exist_ok=True)
            except OSError:
                raise click.BadParameter(f"无法创建输出目录|Cannot create output directory: {dir_path}")
    return dir_path


@click.command(short_help="Cactus泛基因组批量分析工具|Cactus Pangenome Batch Analysis Tool",
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--seqfile', '-s',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='序列文件|Sequence file (两列格式：样本名 + 路径|Two columns: sample_name + path)')
@click.option('--output', '-o',
              required=True,
              callback=lambda ctx, param, value: _validate_directory_exists(value) if value else None,
              help='输出目录|Output directory')
@click.option('--reference', '-r',
              required=True,
              help='参考基因组名称|Reference genome name (必须与seqfile第一个基因组名称匹配|Must match first genome in seqfile)')
@click.option('--singularity',
              default='~/miniforge3/envs/singularity_v.3.8.7/bin/singularity',
              show_default=True,
              help='Singularity可执行文件路径|Singularity executable path')
@click.option('--cactus-sif',
              default='~/software/singularity/cactus_v3.1.4.sif',
              show_default=True,
              help='Cactus SIF镜像路径|Cactus SIF image path')
@click.option('--jobstore',
              default='cactus-jobstore',
              show_default=True,
              help='Toil jobstore目录名称|Toil jobstore directory name')
@click.option('--out-name',
              default='cactus_output',
              show_default=True,
              help='输出文件前缀|Output file prefix')
@click.option('--no-cleanup',
              is_flag=True,
              help='保留jobstore不删除|Keep jobstore without deleting')
@click.option('--formats',
              multiple=True,
              # 注意：HAL和PSA是默认输出，不需要在列表中指定|Note: HAL and PSA are default outputs, no need to specify
              type=click.Choice(['gfa', 'gbz', 'odgi', 'vg', 'vcf', 'xg']),
              default=['gfa', 'gbz', 'odgi'],
              show_default=True,
              help='输出格式|Output formats (可多选|can select multiple)')
@click.option('--threads', '-t',
              type=int,
              default=12,
              show_default=True,
              help='CPU核心数|Number of CPU cores')
@click.option('--max-memory', '-m',
              default='100G',
              show_default=True,
              help='最大内存|Maximum memory')
@click.option('--bind',
              multiple=True,
              help='绑定目录到容器|Bind directory to container (可多次使用|can be used multiple times)')
@click.option('--log-level',
              type=click.Choice(['DEBUG', 'INFO', 'WARN', 'ERROR']),
              default='INFO',
              show_default=True,
              help='日志级别|Log level')
@click.option('--quiet',
              is_flag=True,
              help='静默模式|Quiet mode')
def cactus(**kwargs):
    """
    Cactus泛基因组批量分析工具|Cactus Pangenome Batch Analysis Tool

    使用Cactus构建泛基因组，支持GFA、HAL、GBZ、ODGI等多种输出格式
    Build pangenome using Cactus, supports GFA, HAL, GBZ, ODGI and other output formats

    示例|Examples: biopytools cactus -s seqfile.txt -o output/ -r ref_genome
    """
    # 获取主函数|Get main function
    cactus_main = _lazy_import_cactus_main()

    # 构建参数列表|Build arguments list
    argv = []

    # 必需参数|Required arguments
    argv.extend(['--seqfile', kwargs['seqfile']])
    argv.extend(['--output', kwargs['output']])
    argv.extend(['--reference', kwargs['reference']])

    # Singularity相关|Singularity related
    if kwargs['singularity'] != '~/miniforge3/envs/singularity_v.3.8.7/bin/singularity':
        argv.extend(['--singularity', kwargs['singularity']])

    if kwargs['cactus_sif'] != '~/software/singularity/cactus_v3.1.4.sif':
        argv.extend(['--cactus-sif', kwargs['cactus_sif']])

    # 流程参数|Pipeline parameters
    if kwargs['jobstore'] != 'cactus-jobstore':
        argv.extend(['--jobstore', kwargs['jobstore']])

    if kwargs['out_name'] != 'cactus_output':
        argv.extend(['--out-name', kwargs['out_name']])

    if kwargs['no_cleanup']:
        argv.append('--no-cleanup')

    # 输出格式|Output formats
    if kwargs['formats']:
        argv.extend(['--formats'] + list(kwargs['formats']))

    # 性能参数|Performance parameters
    argv.extend(['--threads', str(kwargs['threads'])])
    argv.extend(['--max-memory', kwargs['max_memory']])

    # 绑定目录|Bind directories
    if kwargs['bind']:
        for bind_path in kwargs['bind']:
            argv.extend(['--bind', bind_path])

    # 日志选项|Logging options
    argv.extend(['--log-level', kwargs['log_level']])

    if kwargs['quiet']:
        argv.append('--quiet')

    # 修改sys.argv并调用主函数|Modify sys.argv and call main function
    sys.argv = ['cactus'] + argv
    cactus_main()
