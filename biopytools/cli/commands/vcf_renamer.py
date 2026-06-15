"""
VCF样本名称重命名命令|VCF Sample Name Renamer Command
"""

import click
import sys
import os


def _lazy_import_vcf_renamer_main():
    """延迟加载vcf_renamer主函数|Lazy load vcf_renamer main function"""
    try:
        from ...vcf_renamer.main import main as vcf_renamer_main
        return vcf_renamer_main
    except ImportError as e:
        click.echo(f"导入错误|Import error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否是帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_input_vcf(file_path):
    """验证输入VCF文件|Validate input VCF file (only in non-help mode)"""
    if _is_help_request():
        return file_path
    if not file_path:
        raise click.BadParameter("输入VCF文件路径不能为空|Input VCF path cannot be empty")
    if not os.path.exists(file_path):
        raise click.BadParameter(f"输入VCF文件不存在|Input VCF not found: {file_path}")
    if not file_path.endswith(('.vcf.gz', '.vcf')):
        raise click.BadParameter(f"输入文件必须是VCF格式|Input must be VCF format: {file_path}")
    return file_path


@click.command(
    short_help='VCF样本名称重命名工具|VCF sample name renamer',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--input', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_input_vcf(value) if value else None,
              help='输入VCF文件路径|Input VCF file path')
@click.option('--output', '-o',
              required=True,
              type=click.Path(),
              help='输出VCF文件路径|Output VCF file path')
@click.option('--prefix', '-p',
              default='S',
              show_default=True,
              help='新样本名前缀|New sample name prefix')
@click.option('--mapping', '-m',
              type=click.Path(),
              help='样本映射文件路径|Mapping file path')
@click.option('--no-mapping',
              is_flag=True,
              help='不保留映射文件|Do not keep mapping file')
def vcf_renamer(input, output, prefix, mapping, no_mapping):
    """
    VCF样本名称重命名工具|VCF Sample Name Renamer Tool

    使用bcftools将VCF样本名重命名为简单序号名称(S1, S2, S3...)|Rename VCF sample names to sequential names using bcftools

    示例|Examples: biopytools vcf-renamer -i input.vcf.gz -o output.vcf.gz
    """

    # 延迟加载|Lazy loading
    vcf_renamer_main = _lazy_import_vcf_renamer_main()

    # 构建参数列表|Build argument list
    args = ['vcf_renamer.py']

    # 必需参数|Required parameters
    args.extend(['-i', input])
    args.extend(['-o', output])

    # 可选参数（只在非默认值时添加）|Optional parameters (add only when non-default)
    if prefix != 'S':
        args.extend(['-p', prefix])

    if mapping:
        args.extend(['-m', mapping])

    if no_mapping:
        args.append('--no-mapping')

    # 保存并恢复sys.argv|Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        # 调用原始main函数|Call original main function
        vcf_renamer_main()
    except SystemExit as e:
        # 处理程序正常退出|Handle normal program exit
        sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n用户中断|User interrupted", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"运行时错误|Runtime error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
