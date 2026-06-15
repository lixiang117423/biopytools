"""
GTX Joint Calling命令生成器|GTX Joint Calling Command Generator
"""

import click
import sys
import os


def _lazy_import_gtx_joint_main():
    """延迟加载gtx_joint主函数|Lazy load gtx_joint main function"""
    try:
        from ...gtx_joint.main import main as gtx_joint_main
        return gtx_joint_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_gtx_exec(file_path):
    """验证GTX可执行文件(仅在非帮助模式)|Validate GTX executable (only in non-help mode)"""
    if _is_help_request():
        return file_path
    if not file_path:
        file_path = '~/software/gtx/bin/gtx'
    file_path = os.path.expanduser(file_path)
    if not os.path.exists(file_path):
        raise click.BadParameter(f"GTX可执行文件不存在|GTX executable not found: {file_path}")
    if not os.access(file_path, os.X_OK):
        raise click.BadParameter(f"GTX文件不可执行|GTX file is not executable: {file_path}")
    return file_path


def _validate_reference_file(file_path):
    """验证参考基因组文件(仅在非帮助模式)|Validate reference file (only in non-help mode)"""
    if _is_help_request():
        return file_path
    if not file_path:
        raise click.BadParameter("参考基因组文件路径不能为空|Reference file path cannot be empty")
    if not os.path.exists(file_path):
        raise click.BadParameter(f"参考基因组文件不存在|Reference file not found: {file_path}")
    if not os.path.exists(f"{file_path}.fai"):
        raise click.BadParameter(
            f"参考基因组索引不存在|Reference index not found: {file_path}.fai\n"
            f"请运行|Please run: samtools faidx {file_path}"
        )
    return file_path


def _validate_gvcf_dir(dir_path):
    """验证GVCF目录(仅在非帮助模式)|Validate GVCF directory (only in non-help mode)"""
    if _is_help_request():
        return dir_path
    if not dir_path:
        raise click.BadParameter("GVCF目录路径不能为空|GVCF directory path cannot be empty")
    if not os.path.exists(dir_path):
        raise click.BadParameter(f"GVCF目录不存在|GVCF directory not found: {dir_path}")
    if not os.path.isdir(dir_path):
        raise click.BadParameter(f"GVCF路径不是目录|GVCF path is not a directory: {dir_path}")
    return dir_path


@click.command(short_help='GTX Joint Calling命令生成器|GTX Joint Calling Command Generator',
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--ref', '-r',
              required=True,
              callback=lambda ctx, param, value: _validate_reference_file(value) if value else None,
              help='参考基因组文件路径|Reference genome file path')
@click.option('--input', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_gvcf_dir(value) if value else None,
              help='GVCF文件目录|GVCF files directory')
@click.option('--output', '-o',
              required=True,
              type=click.Path(),
              help='输出目录|Output directory')
@click.option('--threads', '-t',
              default=12,
              type=int,
              show_default=True,
              help='线程数|Number of threads')
@click.option('--tmp-dir', '-T',
              default='./tmp',
              type=click.Path(),
              show_default=True,
              help='临时目录|Temporary directory')
@click.option('--script', '-s',
              default='run_gtx_joint.sh',
              show_default=True,
              help='输出脚本文件名|Output script filename')
@click.option('--faketime', '-f',
              default='2020-10-20 00:00:00',
              show_default=True,
              help='faketime时间|faketime time')
@click.option('--pattern', '-p',
              help='染色体过滤正则表达式|Chromosome filter pattern (e.g., "^Chr[0-9]+$")')
@click.option('--window', '-w',
              type=int,
              help='区间大小(bp)|Window size in bp (e.g., 10000000 for 10M)')
@click.option('--gtx', '-g',
              default='~/software/gtx/bin/gtx',
              show_default=True,
              callback=lambda ctx, param, value: _validate_gtx_exec(value) if value else None,
              help='GTX可执行文件路径|GTX executable path')
def gtx_joint(ref, input, output, threads, tmp_dir, script, faketime, pattern, window, gtx):
    """
    GTX Joint Calling命令生成工具|GTX Joint Calling Command Generator Tool

    生成GTX Joint Calling的批量执行脚本，支持按染色体或区间拆分
    Generate batch execution scripts for GTX Joint Calling, supports splitting by chromosome or windows

    示例|Examples: biopytools gtx-joint -r genome.fa -i ./gvcf -o ./output
    """

    # 延迟加载|Lazy load: import only when actually called
    gtx_joint_main = _lazy_import_gtx_joint_main()

    # 构建主函数的参数列表|Build argument list for original main function
    args = ['gtx_joint.py']

    # 必需参数|Required parameters
    if gtx != '~/software/gtx/bin/gtx':
        args.extend(['-g', gtx])
    args.extend(['-r', ref])
    args.extend(['-i', input])
    args.extend(['-o', output])

    # 可选参数|Optional parameters (add only when non-default)
    if threads != 88:
        args.extend(['-t', str(threads)])

    if tmp_dir != './tmp':
        args.extend(['-T', tmp_dir])

    if script != 'run_gtx_joint.sh':
        args.extend(['-s', script])

    if faketime != '2020-10-20 00:00:00':
        args.extend(['-f', faketime])

    if pattern:
        args.extend(['-p', pattern])

    if window:
        args.extend(['-w', str(window)])

    # 保存并恢复sys.argv|Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        # 调用原始主函数|Call original main function
        gtx_joint_main()
    except SystemExit as e:
        # 处理正常程序退出|Handle normal program exit
        sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n用户中断|User interrupted", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"执行失败|Execution failed: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
