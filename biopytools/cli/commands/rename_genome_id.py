"""
基因组ID重命名命令|Genome ID Renamer Command
"""

import click
import sys
import os


def _lazy_import_rename_genome_id_main():
    """延迟加载rename_genome_id主函数|Lazy load rename_genome_id main function"""
    try:
        from ...fasta_id_renamer.main import main as rename_genome_id_main
        return rename_genome_id_main
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
    short_help='基因组ID重命名工具 - 按顺序重命名|Genome ID Renamer - Rename sequences in order',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--input', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='输入FASTA文件|Input FASTA file')
@click.option('--output', '-o',
              required=True,
              help='输出FASTA文件|Output FASTA file')
@click.option('--prefix', '-p',
              default='Chr',
              show_default=True,
              help='序列前缀|Sequence prefix')
@click.option('--no-zero-padding',
              is_flag=True,
              help='不使用零填充(如Chr1而非Chr01)|Do not use zero padding')
@click.option('--padding-width', '-w',
              type=int,
              default=2,
              show_default=True,
              help='填充宽度|Padding width')
@click.option('--chr-count', '-n',
              type=int,
              default=0,
              show_default=True,
              help='染色体数量，提取前N条作为染色体|Chromosome count, extract first N as chromosomes')
@click.option('--no-mapping',
              is_flag=True,
              help='不保存ID映射文件|Do not save ID mapping file')
@click.option('--mapping-file',
              default=None,
              help='ID映射文件路径|ID mapping file path')
@click.option('--log-level',
              type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'ERROR']),
              default='INFO',
              show_default=True,
              help='日志级别|Log level')
def rename_genome_id(input, output, prefix, no_zero_padding, padding_width,
                     chr_count, no_mapping, mapping_file, log_level):
    """
    基因组ID重命名工具|Genome ID Renamer

    按顺序重命名所有序列为Chr01, Chr02, ..., 并可选择性提取前N条作为染色体
    Rename all sequences in order as Chr01, Chr02, ..., optionally extract first N as chromosomes

    示例|Examples: biopytools rename-genome-id -i input.fa -o output.fa -n 20
    """

    # 延迟加载|Lazy loading
    rename_genome_id_main = _lazy_import_rename_genome_id_main()

    # 构建参数列表|Build argument list
    args = ['fasta_id_renamer.py']

    # 必需参数|Required parameters
    args.extend(['-i', input])
    args.extend(['-o', output])

    # 重命名规则|Renaming rules
    if prefix != 'Chr':
        args.extend(['-p', prefix])

    if no_zero_padding:
        args.append('--no-zero-padding')

    if padding_width != 2:
        args.extend(['-w', str(padding_width)])

    # 染色体提取|Chromosome extraction
    if chr_count > 0:
        args.extend(['-n', str(chr_count)])

    # ID映射|ID mapping
    if no_mapping:
        args.append('--no-mapping')

    if mapping_file:
        args.extend(['--mapping-file', mapping_file])

    # 日志级别|Log level
    if log_level != 'INFO':
        args.extend(['--log-level', log_level])

    # 执行主程序|Execute main program
    original_argv = sys.argv
    sys.argv = args

    try:
        rename_genome_id_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
